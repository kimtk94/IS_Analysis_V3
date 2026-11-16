#!/usr/bin/env python3
"""Inspect ZIP central directory of public Japanese immune cis-eQTL resource.

Avoid a full ~1.8 GB ZIP download. Request only a bounded tail range and
(optional) central directory slice. No genotype data or participant IDs.
A ZIP filename cannot establish all nominal cis variant effects / eQTL
colocalization eligibility. Keep all scientific gates HOLD.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import struct
from urllib.parse import urlparse
from urllib.request import Request,urlopen
import urllib.error

SOURCE_URL="https://humandbs-production.ddbj.nig.ac.jp/files/hum0099/hum0099.v1.eqtl.v1.zip"
MAX_TAIL_BYTES=262144
MAX_CENTRAL_DIRECTORY=1048576
EOCD_SIGNATURE=b"PK\x05\x06"
CD_SIGNATURE=b"PK\x01\x02"


def source_url_ok(uri):
    p=urlparse(uri)
    if p.scheme!="https" or p.hostname!="humandbs-production.ddbj.nig.ac.jp":
        raise ValueError("Unexpected Japanese published data host")


def ranged_bytes(uri,bytes_range,limit):
    source_url_ok(uri)
    req=Request(uri,headers={
        "Range":"bytes="+bytes_range,
        "User-Agent":"IS_Analysis_V3_NHA000051_ZIP_manifest/1.0",
        "Accept-Encoding":"identity",
    })
    with urlopen(req,timeout=65) as result:
        content=result.read(limit+1)
        if len(content)>limit:
            raise ValueError("Remote ranged response exceeded declared limit")
        # Only 206 validates a byte-offset relative to the full original source.
        if result.status!=206:
            raise ValueError("HTTP Range unsupported: "+str(result.status))
        hdr=result.headers.get("Content-Range","")
        match=re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)",hdr)
        if not match:
            raise ValueError("Missing valid Content-Range for original source")
        first,last,total=[int(x) for x in match.groups()]
        if last-first+1!=len(content) or total<1 or not first<=last<total:
            raise ValueError("Inconsistent published ZIP source byte range")
        return content,{"start":first,"end":last,"total":total,"bytes_read":len(content)}


def parse_eocd(tail,absolute_start,total):
    offset=tail.rfind(EOCD_SIGNATURE)
    if offset<0 or offset+22>len(tail):
        raise ValueError("ZIP EOCD not found in last bytes")
    (signature,disk_no,disk_cd,disk_entries,total_entries,
     cd_size,cd_start,comment_len)=struct.unpack_from("<IHHHHIIH",tail,offset)
    if (disk_no!=0 or disk_cd!=0 or disk_entries!=total_entries
            or offset+22+comment_len>len(tail)):
        raise ValueError("Multivolume ZIP or malformed EOCD unsupported")
    if cd_size==0xffffffff or cd_start==0xffffffff or total_entries==0xffff:
        raise ValueError("ZIP64 directory requires additional source-specific implementation")
    if not 0<cd_size<=MAX_CENTRAL_DIRECTORY or cd_start+cd_size>total:
        raise ValueError("Implausible ZIP central directory")
    return {"central_directory_offset":cd_start,
            "central_directory_bytes":cd_size,
            "central_directory_end":cd_start+cd_size-1,
            "central_directory_entries":total_entries,
            "zip_comment_bytes":comment_len,
            "no_zip_members_downloaded":True}


def list_zip_members(central_bytes,entry_limit):
    offset=0
    names=[]
    total=0
    for _ in range(entry_limit):
        if offset+46>len(central_bytes):
            raise ValueError("Truncated ZIP directory")
        fields=struct.unpack_from("<IHHHHHHIIIHHHHHII",central_bytes,offset)
        if fields[0]!=0x02014B50:
            raise ValueError("Malformed ZIP central directory signature")
        nlen,elen,clen=fields[10],fields[11],fields[12]
        compressed,uncompressed=fields[8],fields[9]
        p=offset+46
        name=central_bytes[p:p+nlen].decode("utf-8","replace")
        if not name or len(name)>400 or not p+nlen+elen+clen<=len(central_bytes):
            raise ValueError("Malformed ZIP filename/extension")
        names.append({"member":name,
                      "compressed_bytes":compressed,
                      "uncompressed_bytes":uncompressed,
                      "member_data_not_fetched":True})
        total+=compressed
        offset=p+nlen+elen+clen
    if offset!=len(central_bytes):
        raise ValueError("Unexpected data after described ZIP central directory")
    return {"files":names,"sum_file_compressed_bytes":total,
            "all_zip_member_names_covered":True}


def run(out):
    report={
      "audited_at_utc":datetime.now(timezone.utc).isoformat(),
      "dataset":"NHA000051 Japanese healthy donor immune cell eQTL",
      "original_archive_url":SOURCE_URL,
      "archive_access":"SOURCE_UNAVAILABLE",
      "actual_nominal_assoc_rows_examined":False,
      "all_cis_pairs_confirmed":False,
      "effect_SE_and_EAF_confirmed":False,
      "genome_build_independently_reconciled_to_GWAS":False,
      "compatible_with_coloc":"NOT_YET_VERIFIED",
      "research_use":"PUBLIC_SOURCE_STRUCTURE_ONLY",
    }
    try:
        tail,range_status=ranged_bytes(SOURCE_URL,f"-{MAX_TAIL_BYTES}",MAX_TAIL_BYTES)
        report["ZIP_footer_range"]=range_status
        d=parse_eocd(tail,range_status["start"],range_status["total"])
        report["ZIP_central_directory"]=d
        start=d["central_directory_offset"]
        end=d["central_directory_end"]
        rel=start-range_status["start"]
        if rel>=0 and end<=range_status["end"]:
            directory=tail[rel:rel+d["central_directory_bytes"]]
        else:
            directory,meta=ranged_bytes(SOURCE_URL,f"{start}-{end}",MAX_CENTRAL_DIRECTORY)
            if meta["start"]!=start or meta["end"]!=end:
                raise ValueError("Central directory range mismatches EOCD")
        listing=list_zip_members(directory,d["central_directory_entries"])
        report["ZIP_member_manifest"]=listing
        report["archive_access"]="ACTUAL_PUBLIC_ZIP_CENTRAL_DIRECTORY_VERIFIED"
    except (urllib.error.URLError,OSError,TimeoutError,ValueError,struct.error) as exc:
        report["archive_access"]="SOURCE_OR_UNSUPPORTED_ZIP_HOLD"
        report["error_short"]=f"{type(exc).__name__}: {str(exc)[:250]}"
    out.mkdir(parents=True,exist_ok=True)
    (out/"NHA000051_original_zip_directory.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    lines=[
      "# Japanese NHA000051 full public archive central-directory inventory",
      "",
      "**Archive source status:** "+report["archive_access"],
      "",
      "No complete eQTL associations have been fetched. Alleles/beta/SE, N, coverage and statistical completeness remain HOLD.",
    ]
    if directory:=report.get("ZIP_central_directory"):
        lines.extend(["",f"Original archive members: {directory['central_directory_entries']}","",
                      "| ZIP archive member | Compressed bytes |",
                      "|---|---:|"])
        for item in report["ZIP_member_manifest"]["files"][:35]:
            lines.append(f"| {item['member']} | {item['compressed_bytes']:,} |")
    if err:=report.get("error_short"):
        lines.append(f"\nSource ZIP manifest issue: {err}")
    md="\n".join(lines)+"\n"
    (out/"summary.md").write_text(md,encoding="utf-8")
    print(md)
    if env:=os.getenv("GITHUB_STEP_SUMMARY"):
        with open(env,"a",encoding="utf-8") as f:
            f.write(md)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",default="artifacts/is_eas_japanese_zip_manifest")
    a=p.parse_args()
    run(Path(a.out_dir))


if __name__=="__main__":
    main()
