#!/usr/bin/env python3
"""Stream only the published BBJ 2019 CI ZIP's IS autosomal member over HTTP ranges.
Verify ZIP central directory, member CRC and size, gzip integrity and exact GRCh37
position/ref/alt. Preserve scientific HOLD even when numeric association is found.
"""
from __future__ import annotations
import argparse,gzip,hashlib,io,json,os,re,struct,tempfile,urllib.request,zlib
from pathlib import Path

URL="https://humandbs-production.ddbj.nig.ac.jp/files/hum0014/hum0014.v17.CI.v1.zip"
MEMBER="hum0014.v17.CI.v1/IS.auto.rsq07.mac10.txt.gz"
EXPECTED_SIZE=2252105310
CHR="19";POS=15573379;REF="T";ALT="C";TARGET_RSID="rs187759247"
CHUNK=32*1024*1024
def get_range(start,end):
    req=urllib.request.Request(URL,headers={"Range":f"bytes={start}-{end}","User-Agent":"BBJ2019-one-variant-verified/1.0"})
    with urllib.request.urlopen(req,timeout=180) as r:
        if r.status!=206 or r.geturl()!=URL: raise ValueError("unapproved URL or range response")
        if r.headers.get("Content-Range")!=f"bytes {start}-{end}/{EXPECTED_SIZE}":
            raise ValueError("Content-Range mismatch")
        b=r.read(end-start+2)
    if len(b)!=end-start+1:raise ValueError("truncated/mismatched range")
    return b
def zip_metadata():
    req=urllib.request.Request(URL,method="HEAD",headers={"User-Agent":"BBJ2019-one-variant-verified/1.0"})
    with urllib.request.urlopen(req,timeout=45) as r:
        if r.status!=200 or r.geturl()!=URL or int(r.headers.get("Content-Length","0"))!=EXPECTED_SIZE:
            raise ValueError("source ZIP size/URL changed")
    tail=get_range(EXPECTED_SIZE-262144,EXPECTED_SIZE-1)
    ix=tail.rfind(bytes([80,75,5,6]))
    if ix<0:raise ValueError("ZIP end missing")
    disk,other,nd,total,sz,offset,comment=struct.unpack_from("<4H2LH",tail,ix+4)
    if disk or other or total!=nd or sz>5000000 or offset+sz>EXPECTED_SIZE:
        raise ValueError("ZIP layout unexpected")
    cd=get_range(offset,offset+sz-1)
    off=0;found=None
    for _ in range(total):
        if cd[off:off+4]!=b"PK\\x01\\x02".replace(b"\\x01",bytes([1])).replace(b"\\x02",bytes([2])):
            raise ValueError("ZIP central entry invalid")
        v=struct.unpack_from("<4s6H3L5H2L",cd,off)
        compression=v[4];crc,compressed,uncompressed=v[7:10]
        nlen,elen,clen=v[10:13];local=v[16]
        name=cd[off+46:off+46+nlen].decode("utf-8")
        if name==MEMBER:
            if found:raise ValueError("duplicate target ZIP member")
            found={"name":name,"compression_method":compression,"crc32":crc,
                   "compressed_bytes":compressed,"uncompressed_bytes":uncompressed,
                   "local_header_offset":local}
        off+=46+nlen+elen+clen
    if off!=sz or not found:raise ValueError("target ZIP member not found")
    return found
def scan():
    m=zip_metadata()
    local=get_range(m["local_header_offset"],m["local_header_offset"]+29)
    if local[:4]!=bytes([80,75,3,4]):raise ValueError("invalid local ZIP header")
    namelen,extralen=struct.unpack_from("<HH",local,26)
    start=m["local_header_offset"]+30+namelen+extralen
    if m["compressed_bytes"]>1300000000 or start+m["compressed_bytes"]>EXPECTED_SIZE:
        raise ValueError("source member size limit")
    if m["compression_method"] not in (0,8):raise ValueError("unsupported ZIP compression")
    decompressor=zlib.decompressobj(-zlib.MAX_WBITS) if m["compression_method"]==8 else None
    sha=hashlib.sha256();crc=0;count=0
    with tempfile.TemporaryDirectory(prefix="bbj2019_") as work:
        fp=Path(work)/"IS.auto.rsq07.mac10.txt.gz"
        with fp.open("wb") as out:
            for n in range(0,m["compressed_bytes"],CHUNK):
                b=get_range(start+n,start+min(n+CHUNK,m["compressed_bytes"])-1)
                raw=decompressor.decompress(b) if decompressor else b
                crc=zlib.crc32(raw,crc);sha.update(raw);out.write(raw);count+=len(raw)
                if count>1300000000:raise ValueError("decompressed member cap exceeded")
            if decompressor:
                raw=decompressor.flush()
                crc=zlib.crc32(raw,crc);sha.update(raw);out.write(raw);count+=len(raw)
                if not decompressor.eof or decompressor.unused_data:raise ValueError("incomplete ZIP member DEFLATE")
        if count!=m["uncompressed_bytes"] or crc!=m["crc32"]:
            raise ValueError("original ZIP member CRC32 or size mismatches")
        hits=[];rows=0;header=[]
        with gzip.open(fp,"rt",encoding="utf-8") as handle:
            for line in handle:
                if not header:
                    header=line.split()
                    print("SOURCE_HEADER",json.dumps(header))
                    if not {"CHR","POS","Allele1","Allele2","BETA","SE","p.value"}.issubset(header):
                        raise ValueError("unexpected BBJ source schema")
                    index={name:i for i,name in enumerate(header)}
                    continue
                fields=line.split()
                rows+=1
                if len(fields)!=len(header):raise ValueError("source record field count drift")
                if fields[index["CHR"]]!=CHR or fields[index["POS"]]!=str(POS):continue
                hits.append(dict(zip(header,fields)))
                if len(hits)>1:raise ValueError("duplicate source coordinate")
        result={"technical_status":"PASS_VERIFIED_MEMBER_AND_SINGLE_COORDINATE" if hits else
                "HOLD_NO_EXACT_BBJ_2019_CI_COORDINATE","source_zip_bytes":EXPECTED_SIZE,
                "source_member":m["name"],"source_member_sha256":sha.hexdigest(),
                "member_crc32_verified":True,"gzip_full_stream_verified":True,
                "original_rows_scanned":rows,"exact_coordinate_records":len(hits),
                "target":"19:15573379:T:C","scientific_status":"HOLD",
                "sample_independence_confirmed":False,
                "OpenGWAS_bbj_a_129_data_identical_not_assumed":True}
        if hits:
            hit=hits[0]
            result["original_association"]=hit
            result["allele_alignment_status"]=("PASS_REF_ALT" if
                hit.get("Allele1")==REF and hit.get("Allele2")==ALT
                else "HOLD_ALLELES_NOT_DIRECTLY_ALIGNED")
        return result
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path);a=p.parse_args()
    try:result=scan()
    except Exception as exc:
        result={"technical_status":"HOLD_SOURCE_OR_PARSER","error_type":type(exc).__name__,
                "error":str(exc)[:300],"scientific_status":"HOLD"}
    result.update({"schema":"ndbc-BBJ2019-IS-single-variant-full-source-v1",
                   "dataset":"NHA000096","source":URL})
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({"status":result["technical_status"],
                      "rows":result.get("original_rows_scanned"),
                      "matches":result.get("exact_coordinate_records"),
                      "association":result.get("original_association"),
                      "error":result.get("error"),"scientific_status":"HOLD"}))
if __name__=="__main__":main()
