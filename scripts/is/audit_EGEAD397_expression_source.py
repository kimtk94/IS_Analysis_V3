#!/usr/bin/env python3
"""Read original ImmuNexUT E-GEAD-397 expression ZIP via short HTTPS ranges.

E-GEAD-397 is study expression metadata/counts; it may offer sample-by-celltype
evidence. This does NOT license equating expression sample count to matched
per-cell-type regression N. Never export participant IDs or row-level metadata.
Original 286MB ZIP MD5 not claimed checked by a source byte-prefix probe.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
from datetime import datetime,timezone
import hashlib
import io
import json
import os
from pathlib import Path
import re
import struct
from urllib.request import Request,urlopen
import urllib.error
import zlib

BASE="https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-000/E-GEAD-397/"
MANIFEST=BASE+"E-GEAD-397.filelist.txt"
ARCHIVE=BASE+"E-GEAD-397.processed.zip"
EXPECTED_BYTES=286026411
EXPECTED_MD5="8255d3aa3c0ccb6f703bd98d8593ffe7"
CAP=2*1024*1024


def read_public(url,byte_range=None,max_bytes=CAP):
    if not url.startswith(BASE):
        raise ValueError("Unexpected public source URL")
    h={"User-Agent":"IS-Analysis-V3-ImmuNexUT-expression-metadata-source/1",
       "Accept-Encoding":"identity"}
    if byte_range:h["Range"]="bytes="+byte_range
    with urlopen(Request(url,headers=h),timeout=65) as response:
        content=response.read(max_bytes+1)
        response_status=response.status
        cr=response.headers.get("Content-Range","")
    if len(content)>max_bytes:
        raise ValueError("Source response exceeds hard byte cap")
    if byte_range:
        if response_status!=206:
            raise ValueError("Original ZIP source did not honor byte range")
        m=re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)",cr)
        if not m:
            raise ValueError("Source has no valid Content-Range")
        start,end,total=map(int,m.groups())
        if total!=EXPECTED_BYTES or len(content)!=end-start+1:
            raise ValueError("Original ZIP length mismatch")
        if byte_range.startswith("-"):
            if end!=total-1:
                raise ValueError("Incorrect ZIP suffix range")
        else:
            a,b=map(int,byte_range.split("-"))
            if (start,end)!=(a,b):
                raise ValueError("Incorrect ZIP explicit offset")
        return content,{"start":start,"end":end,"total":total}
    if response_status!=200:
        raise ValueError("Source metadata not HTTP 200")
    return content,None


def filelist(data):
    d={}
    for line in data.decode("utf-8","strict").splitlines():
        if not line or line.startswith("#"):continue
        columns=line.split("\t")
        if len(columns)!=5:raise ValueError("Bad original filelist")
        kind,name,when,length,md5=columns
        if not length.isdecimal() or not re.fullmatch(r"[0-9a-f]{32}",md5):
            raise ValueError("Bad source manifest MD5/length")
        d[name]={"publisher_md5":md5,"source_length":int(length),"type":kind}
    zip_info=d.get("E-GEAD-397.processed.zip",{})
    if zip_info.get("publisher_md5")!=EXPECTED_MD5 or zip_info.get("source_length")!=EXPECTED_BYTES:
        raise ValueError("Published expression archive metadata drift")
    for name in ("clinical_diagnosis_age_sex_v2.txt","count.tar"):
        if name not in d:raise ValueError("Missing official expression source member")
    return d


def parse_zip_index(tail,meta):
    ix=tail.rfind(b"PK\x05\x06")
    if ix<0 or ix+22>len(tail):
        raise ValueError("Original ZIP EOCD missing")
    f=struct.unpack_from("<IHHHHIIH",tail,ix)
    _,disk,disk_cd,n_disk,n_all,cd_bytes,cd_offset,comment=f
    if disk or disk_cd or n_disk!=n_all or n_all<2 or n_all>100 or cd_bytes>CAP:
        raise ValueError("Malformed or unsupported ZIP directory")
    if (n_all==0xffff or cd_bytes==0xffffffff or cd_offset==0xffffffff):
        raise ValueError("Unexpected ZIP64 for 286MB expression archive")
    if cd_offset+cd_bytes>meta["total"] or ix+22+comment>len(tail):
        raise ValueError("Bad central directory offsets")
    return {"start":cd_offset,"end":cd_offset+cd_bytes-1,"entries":n_all,"bytes":cd_bytes}


def parse_zip_members(directory,count):
    files=[]
    offset=0
    for _ in range(count):
        if offset+46>len(directory):raise ValueError("Truncated original ZIP index")
        cols=struct.unpack_from("<IHHHHHHIIIHHHHHII",directory,offset)
        if cols[0]!=0x02014b50:raise ValueError("Bad ZIP central directory signature")
        namelen,extralen,commentlen=cols[10],cols[11],cols[12]
        end=offset+46+namelen+extralen+commentlen
        if end>len(directory):raise ValueError("Broken original ZIP member length")
        name=directory[offset+46:offset+46+namelen].decode("utf-8","strict")
        files.append({"name":name,"method":cols[4],"compressed_bytes":cols[8],
                      "uncompressed_bytes":cols[9],"crc32":cols[7],"offset":cols[16]})
        offset=end
    if offset!=len(directory):raise ValueError("Unparsed ZIP directory bytes")
    return files


def get_member_prefix(member,max_bytes=CAP):
    at=member["offset"]
    hdr,_=read_public(ARCHIVE,f"{at}-{at+29}",30)
    raw=struct.unpack("<IHHHHHIIIHH",hdr)
    if raw[0]!=0x04034b50 or raw[3]!=member["method"]:
        raise ValueError("ZIP local header differs from central index")
    name_len,extra_len=raw[-2],raw[-1]
    loc,_=read_public(ARCHIVE,f"{at+30}-{at+29+name_len+extra_len}",65536)
    if loc[:name_len].decode("utf-8","strict")!=member["name"]:
        raise ValueError("ZIP member name mismatch")
    start=at+30+name_len+extra_len
    n=min(max_bytes,member["compressed_bytes"])
    data,_=read_public(ARCHIVE,f"{start}-{start+n-1}",max_bytes)
    if member["method"]==0:
        plaintext=data
    elif member["method"]==8:
        decompressor=zlib.decompressobj(-15)
        plaintext=decompressor.decompress(data,4*CAP)
        if len(data)==member["compressed_bytes"] and (
            not decompressor.eof or zlib.crc32(plaintext)&0xffffffff!=member["crc32"]):
            raise ValueError("Entire small ZIP member CRC/deflate mismatch")
    else:raise ValueError("Unsupported ZIP member encoding")
    return plaintext,{"compressed_source_bytes_read":n,
                      "fully_fetched":n==member["compressed_bytes"],
                      "whole_zip_publisher_md5_verified":False}


def summarize_clinical(data):
    text=data.decode("utf-8-sig","strict")
    rows=list(csv.DictReader(io.StringIO(text),delimiter="\t"))
    if not rows:raise ValueError("Empty expression clinical source")
    fields=list(rows[0])
    if len(set(fields))!=len(fields):raise ValueError("Duplicate clinical metadata header")
    # Nonidentifying metadata only. Never output row values or specimen names.
    cell_cols=[k for k in fields if "cell" in k.lower() or "subset" in k.lower()]
    return {
       "public_metadata_rows":len(rows),"header_fields":fields,
       "celltype_columns_detected":cell_cols,
       "public_clinical_member_contains_direct_eQTL_per_celltype_N":False,
       "no_subject_identifiers_or_rows_exported":True,
    }


def tar_first_member(data):
    if len(data)<512:raise ValueError("No count TAR first member")
    h=data[:512]
    name=h[:100].split(b"\0",1)[0].decode("utf-8","replace")
    size=h[124:136].strip(b" \0")
    if not re.fullmatch(rb"[0-7]+",size):
        raise ValueError("Malformed count TAR member header")
    return {"first_member_name":name,
            "declared_first_member_uncompressed_bytes":int(size,8),
            "full_expression_count_matrix_rows_read":False,
            "per_cell_type_N_not_inferred_from_archive_name":True}


def run(out):
    result={"source":"E-GEAD-397 public ImmuNexUT expression archive",
            "checked_utc":datetime.now(timezone.utc).isoformat(),
            "status":"SOURCE_OR_QC_HOLD",
            "true_per_celltype_eQTL_N":"NOT_VERIFIED",
            "original_expression_ZIP_full_MD5":"NOT_VERIFIED",
            "original_expression_matrix_complete":"NOT_VERIFIED"}
    try:
        official=filelist(read_public(MANIFEST,max_bytes=6000)[0])
        tail,metadata=read_public(ARCHIVE,byte_range="-262144",max_bytes=262144)
        cd=parse_zip_index(tail,metadata)
        if cd["start"]>=metadata["start"] and cd["end"]<=metadata["end"]:
            index=tail[cd["start"]-metadata["start"]:cd["start"]-metadata["start"]+cd["bytes"]]
        else:
            index,_=read_public(ARCHIVE,byte_range=f"{cd['start']}-{cd['end']}",max_bytes=CAP)
        members=parse_zip_members(index,cd["entries"])
        for original in ("clinical_diagnosis_age_sex_v2.txt","count.tar"):
            match=[x for x in members if x["name"].split("/")[-1]==original]
            if len(match)!=1 or match[0]["uncompressed_bytes"]!=official[original]["source_length"]:
                raise ValueError("Original expression ZIP internal inventory or member size mismatch")
        clinical=next(x for x in members if x["name"].split("/")[-1]=="clinical_diagnosis_age_sex_v2.txt")
        cbytes,cmeta=get_member_prefix(clinical)
        if not cmeta["fully_fetched"] or hashlib.md5(cbytes).hexdigest()!=official["clinical_diagnosis_age_sex_v2.txt"]["publisher_md5"]:
            raise ValueError("Complete small clinical metadata publisher MD5 mismatch")
        original_counts=next(x for x in members if x["name"].split("/")[-1]=="count.tar")
        header,count_meta=get_member_prefix(original_counts)
        result.update({"status":"EXPRESSION_ARCHIVE_SOURCE_INDEX_AND_CLINICAL_MD5_PASS",
                       "original_expression_archive_members":[x["name"] for x in members],
                       "original_clinical_member_MD5":"PASS",
                       "clinical_structure":summarize_clinical(cbytes),
                       "count_tar_preview":tar_first_member(header),
                       "source_filelist_publisher_count_tar_MD5_ONLY":official["count.tar"]["publisher_md5"]})
    except (OSError,ValueError,TimeoutError,UnicodeDecodeError,urllib.error.URLError,zlib.error):
        result["status"]="SOURCE_OR_QC_HOLD"
        result["issue"]="Public expression archive unavailable, unsupported range or metadata QC mismatch"
    out.mkdir(parents=True,exist_ok=True)
    (out/"E_GEAD397_expression_source_audit.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    lines=["# ImmuNexUT E-GEAD-397 actual public expression source preflight","",
           "**"+result["status"]+"**","",
           "Never substitute 416 cohort WGS donors or expression sample numbers for per-celltype eQTL regression N.",
           ""]
    if "clinical_structure" in result:
        c=result["clinical_structure"]
        lines.extend([f"- Official clinical metadata rows: {c['public_metadata_rows']}",
                      f"- Public clinical header fields: {', '.join(c['header_fields'])}",
                      f"- Actual count TAR first member: {result['count_tar_preview']['first_member_name']}"])
    lines+=["","No original participant IDs emitted, original 286MB archive MD5 NOT independently verified.",
            "True per-cell eQTL N, variant MAF and stroke loci colocalization remain HOLD."]
    summary="\n".join(lines)+"\n"
    (out/"summary.md").write_text(summary)
    print(summary)
    if f:=os.getenv("GITHUB_STEP_SUMMARY"):
        with open(f,"a") as stream:stream.write(summary)
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--out-dir",default="artifacts/is_EGEAD397_expression_source")
    args=parser.parse_args()
    result=run(Path(args.out_dir))
    if result["status"]!="EXPRESSION_ARCHIVE_SOURCE_INDEX_AND_CLINICAL_MD5_PASS":
        raise SystemExit(1)
