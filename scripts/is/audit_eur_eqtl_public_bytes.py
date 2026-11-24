#!/usr/bin/env python3
"""Read-only byte-range proof of public EUR cis-eQTL candidate source formats.

Uses only ~256 KiB from each remote resource, never tries to download the
~2 GB GTEx EUR tar or complete EBI catalogue. Publisher source:
 - GTEx v8 European-American single-tissue eQTL archive, public GCS bucket
 - EBI eQTL Catalogue BLUEPRINT monocyte QTD000021 .all.tsv.gz

An object HEAD or a compressed header is NOT a full nominal cis-eQTL locus
or valid GWAS-eQTL harmonization; report source provenance and explicit holds.
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
import zlib

GTEX_EUR_URL=(
  "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/"
  "GTEx_Analysis_v8_eQTL_EUR.tar"
)
EBI_BLUEPRINT_URL=(
  "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
  "QTS000002/QTD000021/QTD000021.all.tsv.gz"
)
MAX_PREFIX_BYTES=262144
SOURCE_HOSTS={"storage.googleapis.com","ftp.ebi.ac.uk"}
REQUIRED_ALL_STATS={
  "variant","pvalue","gene_id","beta","se","chromosome","position",
  "ref","alt","maf"
}
SOURCE_BUILD={
    "GTEx_v8_EUR_eQTL":"GRCh38 (GTEx v8 QTL variant ID; verify source archive README)",
    "eQTL_Catalogue_BLUEPRINT":"GRCh38 (catalogue official column specification)",
    "GIGASTROKE_EUR_EAS":"GRCh37 (source names)"
}


def safe_url(url):
    obj=urlparse(url)
    if obj.scheme!="https" or obj.hostname not in SOURCE_HOSTS:
        raise ValueError("Unapproved public source host or insecure URL")


def remote_head(url):
    safe_url(url)
    result={"source_url":url,"http_ok":False,
            "real_content_examined":False,"statistically_ready":False}
    try:
        req=Request(url,method="HEAD",
                    headers={"User-Agent":"IS-Analysis-V3-EUR-eQTL-source-QC/1.0"})
        with urlopen(req,timeout=45) as response:
            length=response.headers.get("Content-Length","")
            result["http_status"]=response.status
            result["remote_advertised_bytes"]=int(length) if length.isdecimal() else None
            result["accept_ranges"]=response.headers.get("Accept-Ranges")
            result["http_ok"]=response.status==200
            result["final_host"]=urlparse(response.geturl()).hostname
    except (urllib.error.URLError,OSError,TimeoutError,ValueError) as exc:
        result["access_issue"]=f"{type(exc).__name__}: {str(exc)[:180]}"
    return result


def remote_prefix(url,limit=MAX_PREFIX_BYTES):
    safe_url(url)
    if not 64<=limit<=MAX_PREFIX_BYTES:
        raise ValueError("Byte request outside declared cap")
    req=Request(url,headers={
        "User-Agent":"IS-Analysis-V3-EUR-eQTL-source-QC/1.0",
        "Range":f"bytes=0-{limit-1}",
        "Accept-Encoding":"identity"
    })
    with urlopen(req,timeout=70) as response:
        actual=response.status
        if actual not in (200,206):
            raise ValueError("Unexpected source GET response "+str(actual))
        body=response.read(limit)
        if not body:
            raise ValueError("No prefix data")
        if actual==200 and response.headers.get("Content-Length","").isdigit():
            # If the server ignores range requests we still never read past limit.
            advertised=int(response.headers["Content-Length"])
        else:
            advertised=None
        result={
            "http_status":actual,"bytes_inspected":len(body),
            "http_range_respected":actual==206,
            "server_advertised_length_on_200":advertised,
        }
        return body,result


def first_tar_member(prefix):
    """USTAR header only, no tar extraction or attempts to read entire archive."""
    if len(prefix)<512:
        raise ValueError("Insufficient TAR prefix")
    h=prefix[:512]
    name=h[:100].split(b"\0",1)[0].decode("utf-8","replace")
    mode=h[100:108].split(b"\0",1)[0].decode("ascii","replace")
    size_raw=h[124:136].strip(b" \0").decode("ascii","replace")
    magic=h[257:263]
    if not name or any(ord(c)<32 for c in name):
        raise ValueError("Unrecognized TAR first member name")
    if not re.fullmatch(r"[0-7]+",size_raw):
        raise ValueError("Unparseable TAR first-member size")
    size=int(size_raw,8)
    first_member_type=chr(h[156]) if h[156] else "0"
    return {
        "first_member_name":name,
        "first_member_size_bytes":size,
        "first_member_type":first_member_type,
        "ustar_magic_present":magic.startswith(b"ustar"),
        "does_not_enumerate_entire_tar":True,
    }


def gzip_prefix_header(prefix):
    if len(prefix)<4 or not prefix.startswith(b"\x1f\x8b"):
        raise ValueError("Remote eQTL candidate is not gzip/BGZF")
    obj=zlib.decompressobj(16+zlib.MAX_WBITS)
    decoded=obj.decompress(prefix,1_000_000)
    if not decoded:
        raise ValueError("No compressed TSV data in prefix")
    txt=decoded.decode("utf-8","replace")
    header=txt.split("\n",1)[0].strip("\r")
    columns=[s.lstrip("#") for s in header.split("\t")]
    if len(columns)<5:
        raise ValueError("Unrecognized nominal association header")
    if len(set(columns))!=len(columns):
        raise ValueError("Duplicate eQTL field name in source")
    required_present=sorted(REQUIRED_ALL_STATS & set(columns))
    missing=sorted(REQUIRED_ALL_STATS - set(columns))
    return {
        "source_gzip_magic_verified":True,
        "source_header_columns":columns,
        "full_nominal_required_fields_present":not missing,
        "present_required_fields":required_present,
        "missing_required_fields":missing,
        "first_header_bytes_inflated":len(decoded),
        "has_ebiqtl_expected_GRCh38_variant_key":"variant" in columns,
        "still_need_association_rows_and_variant_coverage":True,
        "effect_orientation_from_publisher":"ALT is effect allele per EBI Catalogue documentation (not yet oriented to GWAS)",
    }


def run(out):
    report={
        "checked_at_utc":datetime.now(timezone.utc).isoformat(),
        "task":"GTEx v8 EUR eQTL TAR and eQTL Catalogue BLUEPRINT nominal association source byte-prefix audit",
        "source_genome_builds":SOURCE_BUILD,
        "GIGASTROKE_genome_build_liftover":"NOT_RUN",
        "ancestry_specific_effect_alignment":"NOT_RUN",
        "nominal_cis_region_complete":"NOT_VERIFIED",
        "full_eqtl_file_downloaded":False,
        "coloc_ready":False,
        "MAGMA_ready":False,
        "LD_SuSiE_ready":False,
        "candidate_sources":{},
    }
    for key,url,check in (
        ("GTEx_v8_EUR_archive",GTEX_EUR_URL,"tar"),
        ("EBI_BLUEPRINT_monocyte_all_associations",EBI_BLUEPRINT_URL,"gzip")
    ):
        status=remote_head(url)
        result={"source":key,"head":status,
                "data_access":"SOURCE_NOT_VERIFIED",
                "source_full_cis_eQTL_QC":"HOLD"}
        try:
            chunk,details=remote_prefix(url)
            result["prefix"]=details
            if check=="tar":
                result["tar_first_member"]=first_tar_member(chunk)
            else:
                result["nominal_tsv_schema"]=gzip_prefix_header(chunk)
            result["data_access"]="PUBLIC_SOURCE_BYTE_PREFIX_CONFIRMED"
        except (urllib.error.URLError,OSError,TimeoutError,ValueError,zlib.error) as exc:
            result["data_access"]="SOURCE_PREFIX_UNAVAILABLE"
            result["source_issue"]=f"{type(exc).__name__}: {str(exc)[:200]}"
        report["candidate_sources"][key]=result

    out.mkdir(parents=True,exist_ok=True)
    (out/"eur_eqtl_public_source_byte_qc.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    lines=[
        "# EUR cis-eQTL source preflight: real remote bytes",
        "",
        "Neither source is approved for GWAS colocalization at this stage.",
        "",
        "| Source | Actual public bytes | Source structure | Coloc |",
        "|---|---|---|---|"
    ]
    for key,item in report["candidate_sources"].items():
        structure=("TAR first member" if key.startswith("GTEx") else "Compressed nominal field names")
        lines.append(f"| {key} | {item['data_access']} | {structure} | HOLD |")
    lines.extend(["",
      "EBI eQTL Catalogue positions are GRCh38 and ALT is effect allele. Original IS GWAS is GRCh37: proper liftover and strand/orientation checks mandatory.",
      "GTEx EUR eQTL TAR header confirms archive structure only; no assumption about contained full tested cis associations.",
      "Successful HTTP metadata/header query does not prove tested-variant completeness, relevant tissue/ancestry or downstream causal association.",
    ])
    doc="\n".join(lines)+"\n"
    (out/"summary.md").write_text(doc,encoding="utf-8")
    print(doc)
    if env:=os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(env,"a",encoding="utf-8") as f:
            f.write(doc)
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",default="artifacts/is_eur_cis_eqtl_byte_qc")
    args=p.parse_args()
    run(Path(args.out_dir))


if __name__=="__main__":
    main()
