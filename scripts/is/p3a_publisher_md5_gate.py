#!/usr/bin/env python3
"""Read-only publisher-MD5 provenance check for two original GIGASTROKE GWAS files.

Download exactly the official GRCh37 compressed files (streaming, no file
retention), compare official live YAML MD5, previously pinned publisher MD5,
previously computed full-file SHA256 and HTTP byte count. Does not inspect
individual variants, infer per-SNP N, or certify independent stroke replication.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
import urllib.request
from pathlib import Path

BASE = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
    "GCST90104001-GCST90105000"
)
SOURCES = {
    "EUR": {
        "accession": "GCST90104540",
        "bytes": 283_482_604,
        "publisher_md5_pinned": "1bfe4ae8a5042fb24cb0a562b05b0d2f",
        "full_sha256_pinned": "e5d9809b608650e36a2808ca27bbb9de5b883ff9661d6d6056d02b3f792eacc0",
        "original_rows_in_prior_source_qc": 7_482_032,
    },
    "EAS": {
        "accession": "GCST90104545",
        "bytes": 260_689_208,
        "publisher_md5_pinned": "904e12e33ec3c2eeefe4dc93bec4995c",
        "full_sha256_pinned": "c521d5695414bdedcc84d2ea21d588e330c75a672f36a87c3fd092b4c55687a3",
        "original_rows_in_prior_source_qc": 6_725_216,
    },
}
HEADERS = {"User-Agent": "IS-Academic-Research-MD5-Verify/1.0"}

def ensure(value, reason):
    if not value:
        raise ValueError(reason)

def build_urls(spec):
    acc=spec["accession"]
    ensure(re.fullmatch(r"GCST[0-9]{8}",acc) is not None,"unexpected accession")
    file=f"{acc}_buildGRCh37.tsv.gz"
    return f"{BASE}/{acc}/{file}", f"{BASE}/{acc}/{file}-meta.yaml"

def publisher_metadata(text, spec):
    acc=spec["accession"]
    fields={}
    for key in ("gwas_id","genome_assembly","data_file_name","data_file_md5sum"):
        match=re.search(r"(?m)^\s*"+re.escape(key)+r":\s*[^a-fA-F0-9_-]*?([^\r\n#]+)", text)
        ensure(match is not None, "publisher YAML missing "+key)
        fields[key]=match.group(1).strip()
    ensure(fields["gwas_id"]==acc and fields["genome_assembly"]=="GRCh37"
           and fields["data_file_name"]==f"{acc}_buildGRCh37.tsv.gz",
           "publisher YAML accession/assembly/filename mismatch")
    ensure(fields["data_file_md5sum"]==spec["publisher_md5_pinned"],
           "publisher release MD5 drift")
    return fields

def validate_result(result, ancestry):
    spec=SOURCES[ancestry]
    ensure(result.get("schema")=="is-p3a-publisher-original-byte-integrity-v1"
           and result.get("ancestry")==ancestry
           and result.get("accession")==spec["accession"]
           and result.get("file_size_bytes")==spec["bytes"]
           and result.get("actual_file_md5")==spec["publisher_md5_pinned"]
           and result.get("actual_file_sha256")==spec["full_sha256_pinned"]
           and result.get("publisher_yaml_md5")==spec["publisher_md5_pinned"]
           and result.get("full_gwas_file_verified_with_publisher_md5") is True
           and result.get("full_gwas_sample_size_approved") is False
           and result.get("archived_P3A_GWAS_allele_effect_verified") is False
           and result.get("scientific_status")=="HOLD",
           "not a validated matching original publisher file or scientific status")
    return True

def verify_remote(ancestry):
    spec=SOURCES[ancestry]
    url, meta=build_urls(spec)
    req=urllib.request.Request(meta,headers=HEADERS)
    with urllib.request.urlopen(req,timeout=45) as response:
        ensure(response.status==200,"publisher YAML HTTP status")
        yaml=response.read(120_001)
        ensure(len(yaml)<=120_000, "publisher metadata too long")
    fields=publisher_metadata(yaml.decode("utf-8-sig"),spec)
    req=urllib.request.Request(url,headers=HEADERS)
    digest_md5=hashlib.md5(usedforsecurity=False)
    digest_sha=hashlib.sha256()
    total=0
    with urllib.request.urlopen(req,timeout=70) as response:
        ensure(response.status==200,"GWAS source must be full HTTP 200")
        declared=response.headers.get("Content-Length")
        if declared:
            ensure(int(declared)==spec["bytes"],"source Content-Length changed")
        while True:
            chunk=response.read(2*1024*1024)
            if not chunk:break
            total+=len(chunk)
            ensure(total<=spec["bytes"],"original GWAS exceeded bounded source size")
            digest_sha.update(chunk)
            digest_md5.update(chunk)
    result={
        "schema":"is-p3a-publisher-original-byte-integrity-v1",
        "ancestry":ancestry,"accession":spec["accession"],
        "reference_accession_from_live_publisher_yaml":fields["gwas_id"],
        "reference_assembly_from_live_publisher_yaml":fields["genome_assembly"],
        "file_size_bytes":total,
        "actual_file_md5":digest_md5.hexdigest(),
        "actual_file_sha256":digest_sha.hexdigest(),
        "publisher_yaml_md5":fields["data_file_md5sum"],
        "earlier_original_variant_rows_reported":spec["original_rows_in_prior_source_qc"],
        "variant_rows_reparsed_this_run":False,
        "full_gwas_file_verified_with_publisher_md5": total==spec["bytes"] and
            digest_sha.hexdigest()==spec["full_sha256_pinned"] and
            digest_md5.hexdigest()==fields["data_file_md5sum"],
        "full_gwas_sample_size_approved":False,
        "archived_P3A_GWAS_allele_effect_verified":False,
        "independent_BBJ_excluded_replication":False,
        "scientific_status":"HOLD",
    }
    validate_result(result,ancestry)
    return result

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--ancestry",choices=list(SOURCES),required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    try:
        result=verify_remote(args.ancestry)
    except Exception as exc:
        result={"schema":"is-p3a-publisher-original-byte-integrity-v1",
                "ancestry":args.ancestry,"scientific_status":"HOLD",
                "technical_status":"HOLD_SOURCE_UNAVAILABLE_OR_DIGEST_MISMATCH",
                "error_type":type(exc).__name__,"error":str(exc)[:240],
                "full_gwas_file_verified_with_publisher_md5":False}
    else:
        result["technical_status"]="PASS_PUBLISHER_MD5_SHA256_SIZE"
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:result.get(k) for k in
                      ("ancestry","technical_status","scientific_status",
                       "file_size_bytes","full_gwas_file_verified_with_publisher_md5","error")},
                     ensure_ascii=False))
    return 0 if result["technical_status"]=="PASS_PUBLISHER_MD5_SHA256_SIZE" else 1

if __name__=="__main__":
    raise SystemExit(main())
