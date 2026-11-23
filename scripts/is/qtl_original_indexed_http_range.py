#!/usr/bin/env python3
"""Read one original eQTL Catalogue Tabix .tbi index and one contiguous BGZF byte range.

Read-only small-source fallback when remote HTSlib FTP/HTTPS transport cannot
open a 4 GB BGZF file on a GitHub hosted runner. No full download, no tabix
process, no automatic retries. The two original publisher HTTP requests are
spaced apart. Returns only streaming-compatible decoded QTL TSV lines and
source integrity hashes; callers MUST never persist raw variant-effect rows.
"""
from __future__ import annotations

import gzip
import hashlib
import io
import os
from pathlib import Path
import re
import struct
import time
import urllib.parse
import urllib.request
import zlib

QTL_PATH=("/pub/databases/spot/eQTL/sumstats/"
          "QTS000025/QTD000434/QTD000434.all.tsv.gz")
HOST="ftp.ebi.ac.uk"
MAX_TBI_COMPRESSED=8_000_000
MAX_BGZF_COMPRESSED_SPAN=125_000_000
MAX_BGZF_UNCOMPRESSED=500_000_000
MIN_RETRY_INTERVAL_SECONDS=3
MAX_TBI_REFS=80
MAX_TBI_BINS_PER_REF=500000
MAX_TBI_OFFSETS_PER_REF=300000


def require(ok,why):
    if not ok:
        raise ValueError(why)


class CheckedReader:
    def __init__(self, raw):
        self.raw=memoryview(raw)
        self.offset=0

    def take(self,amount):
        require(isinstance(amount,int) and 0<=amount<=len(self.raw)-self.offset,
                "Malformed/truncated original tabix index")
        start=self.offset
        self.offset+=amount
        return self.raw[start:self.offset]

    def unpack(self, fmt):
        width=struct.calcsize(fmt)
        return struct.unpack(fmt,self.take(width))


def decode_original_tbi(raw):
    require(raw.startswith(b"\x1f\x8b"),"Original publisher .tbi is not gzip")
    payload=gzip.decompress(raw)
    print("TBI_DECODE_STAGE_1_UNCOMPRESSED_BYTES="+str(len(payload))+
          " MAGIC_HEX="+payload[:4].hex(),flush=True)
    rd=CheckedReader(payload)
    require(rd.take(4).tobytes()==b"TBI\x01","Original index is not TBI v1")
    n_ref,fmt,col_seq,col_beg,col_end,meta,skip,l_nm=rd.unpack("<8i")
    print("TBI_DECODE_STAGE_2_REFCOUNT="+str(n_ref)+
          " FORMAT="+str(fmt)+" SEQ_COL="+str(col_seq)+
          " BEGIN_COL="+str(col_beg)+" END_COL="+str(col_end)+
          " NAME_BYTES="+str(l_nm),flush=True)
    require(1<=n_ref<=MAX_TBI_REFS and 0<l_nm<10000 and
            1<=col_seq<=30 and 1<=col_beg<=30 and 1<=col_end<=30,
            "Original TBI header count/column information invalid")
    names=rd.take(l_nm).tobytes()
    require(names.endswith(b"\x00"),"TBI reference names lack terminating NUL")
    names=[x.decode("ascii") for x in names[:-1].split(b"\x00")]
    require(len(names)==n_ref and len(set(names))==len(names),
            "TBI chromosome name count or identity invalid")
    print("TBI_DECODE_STAGE_3_REFS_PARSED_COUNT="+str(len(names))+
          " HAS_10="+str("10" in names),flush=True)
    refs={}
    for name in names:
        (n_bin,)=rd.unpack("<i")
        require(0<=n_bin<=MAX_TBI_BINS_PER_REF,"Implausible original tabix bin count")
        bins={}
        for _ in range(n_bin):
            (b,)=rd.unpack("<I")
            (n_chunk,)=rd.unpack("<i")
            require(b not in bins and 0<=n_chunk<=1_000_000,
                    "Duplicated TBI bin or implausible original chunk count")
            arr=[]
            for _ in range(n_chunk):
                x,y=rd.unpack("<QQ")
                require(y>=x,"Tabix chunk ends before it starts")
                arr.append((x,y))
            bins[b]=arr
        (n_intv,)=rd.unpack("<i")
        require(0<=n_intv<=MAX_TBI_OFFSETS_PER_REF,
                "Implausible original TBI linear index length")
        linear=[rd.unpack("<Q")[0] for _ in range(n_intv)]
        refs[name]={"bins":bins,"linear_index":linear}
    print("TBI_DECODE_STAGE_4_FINAL_TRAILING_BYTES="+str(len(payload)-rd.offset),
          flush=True)
    require(len(payload)-rd.offset in {0,8},
            "Unexpected trailing original TBI data")
    if len(payload)-rd.offset==8:
        rd.unpack("<Q")  # optional n_no_coor
    return {"chromosomes":refs,"format":fmt,"col_seq":col_seq,
            "col_beg":col_beg,"col_end":col_end,"meta":meta,
            "source_tbi_uncompressed_bytes":len(payload),
            "source_tbi_compressed_SHA256":hashlib.sha256(raw).hexdigest()}


def tbi_bins_for_region(start,end):
    # Tabix reg2bins: 1-based inclusive position to 0-based half-open.
    require(1<=start<=end<2**29,"Original QTL query outside TBI positional capacity")
    beg=start-1
    stop=end-1
    out={0}
    for offset,shift in [(1,26),(9,23),(73,20),(585,17),(4681,14)]:
        out.update(range(offset+(beg>>shift),offset+(stop>>shift)+1))
    return out


def merged_one_range_from_index(tbi,region):
    name=("10" if "10" in tbi["chromosomes"] else
          "chr10" if "chr10" in tbi["chromosomes"] else None)
    require(name is not None,"No chr10 in original TBI")
    ref=tbi["chromosomes"][name]
    require(tbi["col_seq"]==2 and tbi["col_beg"]==3 and
            tbi["col_end"]==3, "Official QTL positional tabix schema drift")
    beg=region["start"]-1
    bins=tbi_bins_for_region(region["start"],region["end"])
    linear=ref["linear_index"]
    min_off=linear[beg>>14] if (beg>>14)<len(linear) else 0
    chunks=[(x,y) for k in bins for x,y in ref["bins"].get(k,[])
            if y>min_off]
    require(chunks,"No indexed chunks for original chr10 region")
    first_virtual=max(min(x for x,y in chunks),min_off)
    last_virtual=max(y for x,y in chunks)
    require(last_virtual>=first_virtual,"Index virtual offsets not monotone")
    start_byte=first_virtual>>16
    end_byte=(last_virtual>>16)+65535
    size=end_byte-start_byte+1
    require(0<size<=MAX_BGZF_COMPRESSED_SPAN,
            "Original 15-gene indexed region exceeds one-request budget")
    return {
        "indexed_chromosome_name":name,
        "index_candidate_bins":len(bins),
        "index_selected_chunks":len(chunks),
        "compressed_range_start":start_byte,
        "compressed_range_end":end_byte,
        "compressed_range_bytes":size,
        "first_virtual_offset_within_bgzf_block":first_virtual & 0xffff,
        "one_contiguous_range_spans_all_candidate_tabix_chunks":True,
    }


def parse_bgzf_range(compressed,first_block_inner_offset):
    # BGZF gzip blocks have a BC extra subfield encoding compressed block size-1.
    require(compressed[:2]==b"\x1f\x8b","Source compressed interval not BGZF gzip")
    p=0
    chunks=[]
    raw_size=0
    n_blocks=0
    while p+18<=len(compressed):
        require(compressed[p:p+3]==b"\x1f\x8b\x08" and
                compressed[p+3]&4, "BGZF source has wrong gzip FEXTRA header")
        xlen=struct.unpack_from("<H",compressed,p+10)[0]
        require(6<=xlen<=128 and p+12+xlen<=len(compressed),
                "BGZF gzip extra field malformed")
        extra_start=p+12
        q=extra_start
        bsize=None
        while q+4<=extra_start+xlen:
            sub_id=compressed[q:q+2]
            (sub_len,)=struct.unpack_from("<H",compressed,q+2)
            q+=4
            require(q+sub_len<=extra_start+xlen,
                    "BGZF extra field extends beyond header")
            if sub_id==b"BC":
                require(sub_len==2,"BGZF BC has invalid subfield size")
                (bsize,)=struct.unpack_from("<H",compressed,q)
            q+=sub_len
        require(bsize is not None,"BGZF compressed block missing BC")
        length=bsize+1
        require(20<=length<=65536,"BGZF declared block size out of range")
        require(p+length<=len(compressed),
                "Original BGZF compressed byte range ended inside last block")
        block=zlib.decompress(compressed[p:p+length],wbits=31)
        require(len(block)<=65536,"BGZF decompressed block oversized")
        if n_blocks==0:
            require(first_block_inner_offset<=len(block),
                    "Tabix virtual within-block position exceeds original block")
            block=block[first_block_inner_offset:]
        chunks.append(block)
        raw_size+=len(block)
        require(raw_size<=MAX_BGZF_UNCOMPRESSED,
                "Original indexed source decompressed output exceeds budget")
        p+=length
        n_blocks+=1
    # Trailing extra bytes from conservative range extension may end mid-BGZF;
    # allow only short partial BGZF after ≥1 complete block, and do not parse it.
    require(n_blocks>=1, "No complete BGZF blocks from original publisher")
    body=b"".join(chunks)
    text=body.decode("utf-8",errors="strict")
    if text and not text.endswith("\n"):
        text=text.rsplit("\n",1)[0]+"\n"
    return {
        "data_lines":io.StringIO(text),
        "source_BGZF_blocks_decompressed":n_blocks,
        "source_BGZF_uncompressed_bytes":raw_size,
        "source_total_QTL_full_file_downloaded":False,
    }


def fetch_original_index(url, opener=urllib.request.urlopen):
    p=urllib.parse.urlsplit(url)
    require(p.scheme=="https" and p.hostname==HOST and p.path==QTL_PATH,
            "Not original publisher source")
    req=urllib.request.Request(url+".tbi",headers={
        "Accept-Encoding":"identity",
        "User-Agent":"Ischemic-Stroke-ROSMAP-15-gene-original-tabix-index/1.0",
    })
    with opener(req,timeout=75) as res:
        require(res.status==200,"Original publisher index HTTP response !=200")
        resolved=urllib.parse.urlsplit(res.geturl())
        require(resolved.scheme=="https" and resolved.hostname==HOST,
                "Index redirected outside original publisher")
        raw=res.read(MAX_TBI_COMPRESSED+1)
    require(0<len(raw)<=MAX_TBI_COMPRESSED,
            "Original publisher .tbi index outside size budget")
    return raw


def fetch_contiguous_range(url,limit,opener=urllib.request.urlopen):
    p=urllib.parse.urlsplit(url)
    require(p.scheme=="https" and p.hostname==HOST and p.path==QTL_PATH,
            "Not exact original eQTL Catalogue QTD000434")
    start,end=limit["compressed_range_start"],limit["compressed_range_end"]
    require(0<=start<=end and end-start+1<=MAX_BGZF_COMPRESSED_SPAN,
            "Publisher source compressed range exceeds request budget")
    req=urllib.request.Request(url,headers={
        "Accept-Encoding":"identity",
        "Range":f"bytes={start}-{end}",
        "User-Agent":"Ischemic-Stroke-ROSMAP-15-gene-one-indexed-range/1.0",
    })
    with opener(req,timeout=150) as res:
        require(res.status==206,
                "Publisher did not honor bounded original source range request")
        resolved=urllib.parse.urlsplit(res.geturl())
        require(resolved.scheme=="https" and resolved.hostname==HOST,
                "Original QTL source redirected to unapproved host")
        content_range=res.headers.get("Content-Range","")
        match=re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)",content_range)
        require(match is not None and
                int(match.group(1))==start and int(match.group(2))==end and
                int(match.group(3))>end,
                "EBI Content-Range does not match tabix index offsets")
        body=res.read(end-start+2)
    require(len(body)==end-start+1,"Original QTL byte range truncated")
    return body,{"publisher_compressed_range_SHA256":hashlib.sha256(body).hexdigest(),
                 "publisher_server_content_range":content_range,
                 "original_full_4GB_source_downloaded":False}


def source_indexed_range(url,region,interval_seconds=MIN_RETRY_INTERVAL_SECONDS):
    # Deterministic stage-only diagnostics: never print original SNP lines.
    print("TBI_STAGE_1_START_ORIGINAL_INDEX_HTTP",flush=True)
    raw_index=fetch_original_index(url)
    if os.environ.get("IS_P3A_WRITE_PUBLIC_INDEX_DEBUG") == "1":
        public_index_path=Path("artifacts/is/QTD000434.published_public.tbi")
        public_index_path.parent.mkdir(parents=True,exist_ok=True)
        public_index_path.write_bytes(raw_index)
        print("TBI_STAGE_ORIGINAL_PUBLIC_INDEX_SAVED_FOR_OFFLINE_PARSER_AUDIT",flush=True)
    print("TBI_STAGE_2_INDEX_BYTES="+str(len(raw_index))+
          " MAGIC_HEX="+raw_index[:4].hex(),flush=True)
    tbi=decode_original_tbi(raw_index)
    print("TBI_STAGE_3_INDEX_DECODE_CHRCOUNT="+str(len(tbi["chromosomes"]))+
          " HAS_10="+str("10" in tbi["chromosomes"]),flush=True)
    limit=merged_one_range_from_index(tbi,region)
    print("TBI_STAGE_4_CONTIGUOUS_BGZF_RANGE_BYTES="+
          str(limit["compressed_range_bytes"]),flush=True)
    # EBI officially warns excessive partial requests can trigger DoS blocks.
    # Do not send more than the original index GET and one source byte-range GET.
    if interval_seconds>0:
        time.sleep(interval_seconds)
    print("TBI_STAGE_5_START_ONLY_QTL_RANGE_HTTP",flush=True)
    body,source=fetch_contiguous_range(url,limit)
    print("TBI_STAGE_6_QTL_RANGE_COMPRESSED_BYTES="+str(len(body)),flush=True)
    stream=parse_bgzf_range(body,limit["first_virtual_offset_within_bgzf_block"])
    print("TBI_STAGE_7_BGZF_PARSE_BLOCKS="+str(stream["source_BGZF_blocks_decompressed"]),flush=True)
    return stream["data_lines"],{
        "method":"ONE_TBI_GET_PLUS_ONE_CONTIGUOUS_HTTPS_RANGE",
        "tabix_index_SHA256":tbi["source_tbi_compressed_SHA256"],
        "tabix_index_uncompressed_bytes":tbi["source_tbi_uncompressed_bytes"],
        "index_metadata_chromosome_count":len(tbi["chromosomes"]),
        "index_and_one_range_query":limit,
        "raw_source_metadata":source,
        "original_BGZF_blocks_decompressed":stream["source_BGZF_blocks_decompressed"],
        "original_BGZF_uncompressed_bytes":stream["source_BGZF_uncompressed_bytes"],
        "individual_QTL_rows_saved_in_artifact":False,
        "whole_QTD000434_source_downloaded":False,
        "requests_to_original_EBI_source":2,
        "automatic_retry_count":0,
    }
