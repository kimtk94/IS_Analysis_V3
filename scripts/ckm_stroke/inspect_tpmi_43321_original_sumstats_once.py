#!/usr/bin/env python3
"""One-shot original TPMI 433.21 PLINK2 GWAS file integrity/header preflight.

Only this small (58.3MB) official source is requested ONCE; no retry or
throttle bypass. Full compressed bytes hashed in ephemeral tmpdir and deleted.
No per-variant beta/SE/P or original whole-file GWAS exported in artifacts.
Do not infer effective case-control N from PLINK2 OBS_CT alone.
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re
import tempfile
import urllib.request
from urllib.parse import urlsplit

URL="https://pheweb.ibms.sinica.edu.tw/download/gwas_plink2/433.21"
HOST="pheweb.ibms.sinica.edu.tw"
EXPECTED_BYTES=58_302_393
SAFETY_CAP=120_000_000
P_ALLOWED=("$","NA","N/A","nan","NaN",".")
EXPECTED_COHORT_CASES=9249
EXPECTED_COHORT_CONTROLS=304660

def read_full_original(source_url=URL,opener=urllib.request.urlopen,
                       expected_bytes=EXPECTED_BYTES):
    u=urlsplit(source_url)
    if u.scheme!="https" or u.hostname!=HOST or u.path!="/download/gwas_plink2/433.21":
        raise ValueError("not the sole authorized original TPMI 433.21 public GWAS link")
    req=urllib.request.Request(source_url,headers={
        "User-Agent":"CKM-Stroke-TPMI-43321-Single-Original-Source-QC/1.0",
        "Accept":"application/octet-stream,application/gzip,text/plain"})
    with tempfile.TemporaryDirectory(prefix="tpmi43321_") as td:
        filename=Path(td)/"original_download_ephemeral"
        downloaded=0
        hashfile=hashlib.sha256()
        with opener(req,timeout=90) as resp:
            loc=urlsplit(resp.geturl())
            if resp.status!=200 or loc.scheme!="https" or loc.hostname!=HOST:
                raise ValueError("original TPMI GWAS download source HTTP status or host changed")
            contentlen=resp.headers.get("Content-Length","")
            if contentlen and contentlen.isdecimal() and int(contentlen)!=expected_bytes:
                raise ValueError("official original download content length changed")
            with filename.open("wb") as dest:
                while True:
                    chunk=resp.read(1024*1024)
                    if not chunk:
                        break
                    downloaded+=len(chunk)
                    if downloaded>SAFETY_CAP:
                        raise ValueError("single-request original GWAS download exceeds cap")
                    hashfile.update(chunk)
                    dest.write(chunk)
        if downloaded!=expected_bytes:
            raise ValueError("full source bytes differ from previous HEAD metadata")
        result=inspect_header_and_sample(filename)
        result["original_URL"]=source_url
        result["exact_original_file_sha256_from_current_download"]=hashfile.hexdigest()
        result["original_file_compressed_bytes"]=downloaded
        result["file_body_downloaded_once_without_retry"]=True
        result["file_deleted_after_source_audit"]=True
        result["download_rate_limit_not_bypassed"]=True
        result["original_effect_beta_SE_P_values_not_published_in_artifact"]=True
        return result

def inspect_header_and_sample(path,sample_rows=200):
    with path.open("rb") as fd:
        magic=fd.read(8)
    if magic.startswith(b"\x1f\x8b"):
        source_type="gzip_text"
        opener=lambda: gzip.open(path,"rt",encoding="utf-8-sig",errors="strict",newline="")
    elif magic.startswith(b"PK") or magic.startswith(b"\x28\xb5\x2f\xfd"):
        raise ValueError("original summary file is an unsupported archive/compression, do not guess parser")
    else:
        source_type="plain_text"
        opener=lambda: path.open("rt",encoding="utf-8-sig",errors="strict",newline="")
    with opener() as handle:
        header=None
        for n,line in enumerate(handle):
            if n>30 or len(line)>16000:
                raise ValueError("source GWAS header exceeds safety cap")
            parts=line.rstrip("\r\n").split("\t")
            if len(parts)<5:
                parts=line.strip().split()
            if len(parts)>=5 and ("P" in parts or "P_VALUE" in parts or "PVAL" in parts or "p" in parts):
                header=parts
                break
        if header is None or len(set(header))!=len(header):
            raise ValueError("TPMI original SNP-level association header missing or duplicated")
        canonical={x.strip().upper().lstrip("#"):i for i,x in enumerate(header)}
        def has(*names):
            return next((canonical[n] for n in names if n in canonical),None)
        chrom=has("CHROM","CHR")
        pos=has("POS","BP","POSITION")
        var=has("ID","SNP","RSID")
        ref=has("REF","OTHER_ALLELE","A2")
        alt=has("ALT","EFFECT_ALLELE","A1")
        pidx=has("P","P_VALUE","PVAL")
        beta=has("BETA","LOG(OR)","LOG_OR","LOGOR")
        se=has("SE","LOG(OR)_SE","BETA_SE","SE_BETA","LOGOR_SE")
        odds=has("OR","ODDS_RATIO")
        nidx=has("OBS_CT","N","N_SAMPLES","SAMPLESIZE","N_TOTAL")
        af=has("AF","A1_FREQ","ALT_FREQS","EAF","FREQ","AF_ALT")
        if chrom is None or pos is None or pidx is None or alt is None:
            raise ValueError("original association file missing variant CHROM POS P and allele")
        nonempty=0
        numeric_p=0
        error_rows=0
        observed_n=[]
        samples=[]
        for line in handle:
            if not line.strip():
                continue
            fields=line.rstrip("\r\n").split("\t")
            if len(fields)<5:
                fields=line.strip().split()
            if len(fields)!=len(header):
                raise ValueError("source original variant row different from header width")
            nonempty+=1
            if nonempty>sample_rows:
                break
            if not fields[chrom] or not fields[pos].isdigit():
                error_rows+=1
                continue
            ptxt=fields[pidx].strip()
            if ptxt not in P_ALLOWED:
                try:
                    p=float(ptxt)
                    if math.isfinite(p) and 0<=p<=1:
                        numeric_p+=1
                    else:
                        error_rows+=1
                except ValueError:
                    error_rows+=1
            if nidx is not None and fields[nidx].isdigit():
                observed_n.append(int(fields[nidx]))
            if len(samples)<2:
                samples.append({"chromosome":fields[chrom],
                                "position_valid":bool(fields[pos].isdigit()),
                                "effect_allele_present":bool(fields[alt]),
                                "variant_source_ID_present":var is not None and bool(fields[var])})
        if nonempty<20 or error_rows>nonempty*.05:
            raise ValueError("insufficient valid original SNP rows for source header QC")
    return {
        "technical_status":"PASS_TPMI_43321_ORIGINAL_GWAS_FULL_BYTE_HASH_AND_BOUNDED_HEADER",
        "file_format_detected":source_type,
        "published_page_endpoint":"433.21",
        "source_phenotype":"Cerebral artery occlusion, with cerebral infarction",
        "reference_genome":"GRCh38_PHEWEB_METHODS",
        "published_phenotype_cases":EXPECTED_COHORT_CASES,
        "published_phenotype_controls":EXPECTED_COHORT_CONTROLS,
        "original_source_columns":header,
        "source_field_coverage":{
           "has_variant_position":chrom is not None and pos is not None,
           "has_original_variant_id":var is not None,
           "has_original_REF":ref is not None,
           "has_original_ALT_or_A1":alt is not None,
           "has_original_P":pidx is not None,
           "has_BETA":beta is not None,
           "has_standard_error":se is not None,
           "has_OR":odds is not None,
           "has_variantwise_observation_count_OBS_CT_or_N":nidx is not None,
           "has_allele_frequency":af is not None,
        },
        "first_bounded_source_records_checked":min(nonempty,sample_rows),
        "first_bounded_source_records_valid_P":numeric_p,
        "invalid_source_QC_records_in_bounded_sample":error_rows,
        "bounded_source_no_effect_statistics_export":samples,
        "OBS_CT_is_NOT_a_case_control_effective_N":True,
        "whole_original_GWAS_effect_rows_not_validated_by_bounded_QC":True,
        "original_GWAS_source_not_identical_to_2026_Neurology_discovery_proven":True,
        "independence_from_GIGASTROKE_participants_established":False,
        "scientific_status":"HOLD",
        "ready_for_ldsc":False,
        "MR_coloc_LDSC_executed":False
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    try:
        r=read_full_original()
    except Exception as exc:
        r={"technical_status":"HOLD_TPMI_43321_SINGLE_DOWNLOAD_OR_GWAS_SCHEMA_INVALID",
           "error_type":type(exc).__name__,"error_detail":str(exc)[:300],
           "scientific_status":"HOLD","ready_for_ldsc":False,
           "MR_coloc_LDSC_executed":False}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n","utf-8")
    print(json.dumps({"status":r["technical_status"],
                      "bytes":r.get("original_file_compressed_bytes"),
                      "format":r.get("file_format_detected"),
                      "fields":r.get("source_field_coverage"),
                      "error":r.get("error_detail"),
                      "scientific_status":r["scientific_status"]}))
    if r["technical_status"]!="PASS_TPMI_43321_ORIGINAL_GWAS_FULL_BYTE_HASH_AND_BOUNDED_HEADER":
        raise SystemExit(1)
if __name__=="__main__":
    main()
