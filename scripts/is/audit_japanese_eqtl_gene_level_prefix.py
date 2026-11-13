#!/usr/bin/env python3
"""Read only bounded byte-ranges of public Japanese immune eQTL nested ZIP.

1) Validate central directory + immutable ZIP CRC / member metadata.
2) Extract only methods DOCX (~22 KiB) and public gene-level tar.gz prefix
   (max 2MiB raw) using local-header offsets, WITHOUT full 1.8GB download.
3) Inspect methods text and first nested TAR names/header to decide if full
   cis-eQTL all-variant beta/SE/allele evidence could be present.

No patient's genotypes, donor IDs or association rows are exported.
GSM/GWAS joint analysis, EAS locus harmonization and coloc stay HOLD.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import io
import json
import os
from pathlib import Path
import re
import struct
from urllib.parse import urlparse
import urllib.error
import xml.etree.ElementTree as ET
import zipfile
import zlib

from audit_japanese_eqtl_zip_central_dir import (
    SOURCE_URL, ranged_bytes, parse_eocd, MAX_TAIL_BYTES,
)

MAX_PREFIX=2*1024*1024
MAX_DOCX=400000


def exact_range(url,start,length,limit):
    if not 0<=start or not 1<=length<=limit:
        raise ValueError("Invalid bounded ZIP source range")
    data,metadata=ranged_bytes(url,f"{start}-{start+length-1}",limit)
    if metadata["start"]!=start or metadata["end"]!=start+length-1 or len(data)!=length:
        raise ValueError("Source byte-range mismatch")
    return data


def zip_directory_entries(data,n):
    cursor=0
    files=[]
    for _ in range(n):
        if cursor+46>len(data):
            raise ValueError("Truncated ZIP central directory")
        f=struct.unpack_from("<IHHHHHHIIIHHHHHII",data,cursor)
        if f[0]!=0x02014B50:
            raise ValueError("Invalid central-directory member signature")
        name_len,extra_len,comment_len=f[10],f[11],f[12]
        end=cursor+46+name_len+extra_len+comment_len
        if end>len(data):
            raise ValueError("Member name exceeds ZIP directory")
        name=data[cursor+46:cursor+46+name_len].decode("utf-8","strict")
        files.append({
            "name":name,"compression_method":f[4],"crc32":f[7],
            "compressed_size":f[8],"uncompressed_size":f[9],
            "local_header_offset":f[16],
            "general_purpose_flags":f[3],
        })
        cursor=end
    if cursor!=len(data):
        raise ValueError("Unaccounted directory bytes")
    return files


def locate_member_data(url,member):
    start=member["local_header_offset"]
    head=exact_range(url,start,30,30)
    f=struct.unpack("<IHHHHHIIIHH",head)
    if f[0]!=0x04034B50:
        raise ValueError("Invalid ZIP local member header")
    if f[3]!=member["compression_method"]:
        raise ValueError("Local/central ZIP compression method differs")
    namelen,extra=f[-2],f[-1]
    meta=exact_range(url,start+30,namelen+extra,10240)
    name=meta[:namelen].decode("utf-8","strict")
    if name!=member["name"]:
        raise ValueError("Local and central member names differ")
    return start+30+namelen+extra


def member_first_uncompressed(url,member,size_limit):
    """Byte-prefix only. zipdeflate member may lack end-of-stream in prefix."""
    if size_limit>MAX_PREFIX or size_limit<512:
        raise ValueError("Unsupported uncompressed prefix cap")
    offset=locate_member_data(url,member)
    method=member["compression_method"]
    if method not in (0,8):
        raise ValueError("Unexpected ZIP member compression method")
    read=min(member["compressed_size"],MAX_PREFIX)
    fragment=exact_range(url,offset,read,MAX_PREFIX)
    if method==8:
        decoded=zlib.decompressobj(-15).decompress(fragment,size_limit)
    else:
        decoded=fragment[:size_limit]
    if not decoded:
        raise ValueError("Unable to decode ZIP member prefix")
    return decoded,{
        "compression_method":method,
        "prefix_member_compressed_bytes_read":read,
        "prefix_member_uncompressed_bytes":len(decoded),
        "only_prefix_not_entire_member":read<member["compressed_size"],
    }


def extract_complete_small_member(url,member):
    if member["compressed_size"]>MAX_DOCX or member["uncompressed_size"]>MAX_DOCX:
        raise ValueError("Methods member exceeds safe complete content cap")
    offset=locate_member_data(url,member)
    src=exact_range(url,offset,member["compressed_size"],MAX_DOCX)
    method=member["compression_method"]
    if method==0:
        decoded=src
    elif method==8:
        z=zlib.decompressobj(-15)
        decoded=z.decompress(src,MAX_DOCX+1)
        if not z.eof:
            raise ValueError("ZIP methods member DEFLATE stream truncated")
    else:
        raise ValueError("Unknown methods ZIP compression")
    if len(decoded)!=member["uncompressed_size"]:
        raise ValueError("Uncompressed methods member size mismatch")
    if zlib.crc32(decoded)&0xFFFFFFFF!=member["crc32"]:
        raise ValueError("Methods document ZIP CRC mismatch")
    return decoded


def first_tar_records_from_nested_gzip(binary_prefix,max_records=14):
    """Decode beginning of a gzip compressed TAR only, no entire archive."""
    if not binary_prefix.startswith(b"\x1f\x8b"):
        raise ValueError("Gene-level member is not gzip bytes")
    stream=zlib.decompressobj(16+zlib.MAX_WBITS)
    inflated=stream.decompress(binary_prefix,2*1024*1024)
    if len(inflated)<512:
        raise ValueError("TAR header unavailable in bounded original prefix")
    names=[]
    off=0
    while off+512<=len(inflated) and len(names)<max_records:
        header=inflated[off:off+512]
        if header==b"\0"*512:
            break
        name=header[:100].split(b"\0",1)[0].decode("utf-8","replace")
        size_raw=header[124:136].strip(b"\0 ")
        if not re.fullmatch(rb"[0-7]+",size_raw):
            raise ValueError("Nested TAR member size invalid")
        size=int(size_raw,8)
        # Presence of truncated large content is expected for prefix-only inspection.
        names.append({
            "name":name,"uncompressed_member_size":size,
            "is_directory":header[156:157]==b"5",
            "member_header_complete":True,
            "association_contents_not_verified":True,
        })
        end=off+512+((size+511)//512)*512
        if end>len(inflated):
            break
        off=end
    return {
        "nested_tar_first_member_headers":names,
        "bytes_of_nested_tar_inflated":len(inflated),
        "full_gene_level_tar_unpacked":False,
        "full_variant_cis_pairs_examined":False,
    }


def methods_docx_text(blob):
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        if "word/document.xml" not in zf.namelist():
            raise ValueError("DOCX does not contain WordprocessingML document")
        xml=zf.read("word/document.xml")
    if len(xml)>2000000:
        raise ValueError("Unexpected methods XML size")
    root=ET.fromstring(xml)
    namespace={"w":"http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    paras=[]
    for item in root.findall(".//w:p",namespace):
        texts=[node.text for node in item.findall(".//w:t",namespace) if node.text]
        if texts:
            paras.append("".join(texts))
    if not paras:
        raise ValueError("Original method document has no text")
    # Preserve bounded material for independent methods review, no patient IDs.
    text="\n".join(paras)
    return {
        "text_length_characters":len(text),
        "paragraphs":len(paras),
        "method_text_first_12000_chars":text[:12000],
        "source_method_docx_crc_verified":True,
    }


def run(out):
    report={
        "checked_utc":datetime.now(timezone.utc).isoformat(),
        "accession":"NHA000051",
        "source":"NBDC unrestricted-access public immune eQTL ZIP",
        "url":SOURCE_URL,
        "archive_gene_level_source_access":"SOURCE_OR_QC_HOLD",
        "read_all_gene_level_associations":False,
        "complete_tested_cis_variant_distribution":"NOT_VERIFIED",
        "GWAS_variant_effect_allele_harmonization":"NOT_RUN",
        "EAS_coloc":"NOT_RUN",
    }
    try:
        tail,meta=ranged_bytes(SOURCE_URL,f"-{MAX_TAIL_BYTES}",MAX_TAIL_BYTES)
        footer=parse_eocd(tail,meta["start"],meta["total"])
        begin=footer["central_directory_offset"]
        length=footer["central_directory_bytes"]
        if begin>=meta["start"] and begin+length-1<=meta["end"]:
            directory=tail[begin-meta["start"]:begin-meta["start"]+length]
        else:
            directory=exact_range(SOURCE_URL,begin,length,1048576)
        members=zip_directory_entries(directory,footer["central_directory_entries"])
        genes=[m for m in members if m["name"].endswith("/eQTL_gene_level.tar.gz")]
        methods=[m for m in members if m["name"].endswith("/eQTL_method_riken_ishigaki.docx")]
        if len(genes)!=1 or len(methods)!=1:
            raise ValueError("Missing expected public source gene-level QTL or methods member")
        report["gene_level_ZIP_member_meta"]=genes[0]
        report["methods_ZIP_member_meta"]=methods[0]
        raw,stage=member_first_uncompressed(SOURCE_URL,genes[0],MAX_PREFIX)
        report["gene_level_prefix_download_meta"]=stage
        report["gene_level_tar_first_members"]=first_tar_records_from_nested_gzip(raw)
        doc=extract_complete_small_member(SOURCE_URL,methods[0])
        report["source_methods"]=methods_docx_text(doc)
        report["archive_gene_level_source_access"]="PUBLIC_GENE_LEVEL_NESTED_GZIP_PREFIX_AND_METHODS_CRC_VERIFIED"
    except (urllib.error.URLError,OSError,TimeoutError,ValueError,
            struct.error,zipfile.BadZipFile,zlib.error,ET.ParseError) as exc:
        report["archive_gene_level_source_access"]="SOURCE_OR_QC_HOLD"
        report["source_issue"]=f"{type(exc).__name__}: {str(exc)[:300]}"
    out.mkdir(parents=True,exist_ok=True)
    (out/"japanese_gene_level_methods_prefix.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    lines=["# Original Japanese NHA000051 gene-level eQTL format and methods pilot","",
           "**"+report["archive_gene_level_source_access"]+"**",""]
    if sub:=report.get("gene_level_tar_first_members"):
        lines.extend(["| First nested TAR member name | Declared uncompressed bytes |",
                      "|---|---:|"])
        for x in sub["nested_tar_first_member_headers"]:
            lines.append(f"| {x['name']} | {x['uncompressed_member_size']:,} |")
    if methods:=report.get("source_methods"):
        lines.append(f"\nMethods DOCX CRC PASS, {methods['paragraphs']} paragraphs available in machine-readable evidence.")
    if err:=report.get("source_issue"):
        lines.append("\nSource HOLD: "+err)
    lines+=["","Only first TAR headers examined. NO variant-gene beta/SE/allele row QC, no tested cis completeness, no GWAS mapping or coloc."]
    md="\n".join(lines)+"\n"
    (out/"summary.md").write_text(md,encoding="utf-8")
    print(md)
    if env:=os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(env,"a",encoding="utf-8") as f:
            f.write(md)
    return report


if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--out-dir",default="artifacts/is_japanese_gene_source")
    args=a.parse_args()
    run(Path(args.out_dir))
