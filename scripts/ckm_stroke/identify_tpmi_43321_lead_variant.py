#!/usr/bin/env python3
"""TPMI 433.21 pinned original: identify the ONE genomewide significant variant.

Single unthrottled 58MB source download, no retry.  Validate SHA256 from
separate source QC run, stream all 1,739,048 records, publish only the
significant *public association* and aggregate positive controls.
No PGS weights; A1 is tested allele. Do not infer Ncase or Ncontrol from OBS_CT.
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import tempfile
import urllib.request
from urllib.parse import urlsplit

from scripts.ckm_stroke.audit_tpmi_43321_full_gwas_qc import (
    EXPECTED_SHA256, EXPECTED_HEADER, URL, HOST, EXPECTED_BYTES, SAFETY_CAP,
)

EXPECTED_ROWS=1_739_048
EXPECTED_OBS=204_099
P_THRESHOLD=5e-8

def parse_significant(columns):
    if len(columns)!=17:
        raise ValueError("original TPMI significant row wrong header width")
    (chrom,pos,vid,ref,alt,prov,a1,omitted,frequency,model,obs,or_text,
     se_text,z_text,p_text,err,rsid)=columns
    if not (chrom.isdecimal() and 1<=int(chrom)<=22
            and pos.isdecimal() and int(pos)>0
            and set(ref).issubset(set("ACGT")) and ref
            and set(alt).issubset(set("ACGT")) and alt
            and a1 in (ref,alt) and model=="ADD" and err=="." and obs.isdecimal()):
        raise ValueError("significant GWAS source association unexpected genotype or error convention")
    try:
        raw_p=float(p_text);or_value=float(or_text)
        se=float(se_text);z=float(z_text);freq=float(frequency)
    except ValueError as exc:
        raise ValueError("non-numeric significant original TPMI association") from exc
    if (not all(math.isfinite(x) for x in (raw_p,or_value,se,z,freq))
        or not 0<=raw_p<=P_THRESHOLD or or_value<=0 or se<=0
        or not 0<=freq<=1 or int(obs)!=EXPECTED_OBS):
        raise ValueError("genomewide significant variant value range or N mismatch")
    beta=math.log(or_value)
    z_from_or=beta/se
    p_from_z=math.erfc(abs(z)/math.sqrt(2))
    delta_z=abs(z_from_or-z)
    if delta_z>0.05:
        raise ValueError("effect direction/logOR-Z mismatch")
    if raw_p>0 and p_from_z>0 and abs(math.log10(raw_p)-math.log10(p_from_z))>0.1:
        raise ValueError("original source P/Z test relation discrepancy")
    return {
       "chrom_GRCh38":chrom,"position_GRCh38":int(pos),
       "source_variant_ID":vid,"reported_rsid":rsid,
       "REF_GRCh38":ref,"ALT_GRCh38":alt,"A1_effect_allele":a1,
       "A1_frequency_original":freq,"source_test":model,
       "original_OBS_CT":int(obs),
       "OR_original_A1":or_value,
       "logOR_original_A1":beta,"LOGOR_SE_original":se,
       "Z_STAT_original":z,"P_original":raw_p,
       "OR_SE_vs_Z_abs_delta":delta_z,
       "P_vs_Z_abs_log10_delta":abs(math.log10(raw_p)-math.log10(p_from_z)) if raw_p>0 and p_from_z>0 else None,
       "allele_tested_A1_NOT_assume_ALT":True,
       "origin":"pinned-original-TPMI433.21-PLINK2-source",
       "variant_id_and_alleles_are_public_GWAS_not_individual_genotype":True,
    }

def scan_original(path):
    found=[];rows=0;max_line=0;count_p0=0
    with gzip.open(path,"rt",encoding="utf-8-sig",newline="") as source:
        head=source.readline().rstrip("\r\n").split("\t")
        if head!=EXPECTED_HEADER:
            raise ValueError("original TPMI file header changed")
        for rawline in source:
            if len(rawline)>100000:
                raise ValueError("unexpected long original GWAS record")
            if not rawline.strip():continue
            fields=rawline.rstrip("\r\n").split("\t")
            if len(fields)!=17:raise ValueError("source GWAS variant schema drift")
            rows+=1
            ptxt=fields[14]
            try:p=float(ptxt)
            except ValueError:raise ValueError("source original P not numeric across full file")
            if not math.isfinite(p) or p<0 or p>1:
                raise ValueError("invalid source original P outside [0,1]")
            if p==0:count_p0+=1
            if p<=P_THRESHOLD:
                found.append(parse_significant(fields))
            if rows>EXPECTED_ROWS:raise ValueError("original row count exceeds original QC")
    if rows!=EXPECTED_ROWS or len(found)!=1 or count_p0!=0:
        raise ValueError("original all-row count, single genomewide hit or P=0 differs from independent QC")
    return {
        "schema":"ckm-TPMI43321-pinned-source-genomewide-significant-lead-v1",
        "technical_status":"PASS_PINNED_1739048_SOURCE_ROWS_ONE_SIGNIFICANT_VARIANT",
        "all_original_rows_verified":rows,"genomewide_threshold":P_THRESHOLD,
        "source_genomewide_significant_count":len(found),
        "source_P_eq_zero_count":count_p0,
        "original_source_genome_build":"GRCh38",
        "original_significant_associations":found,
        "source_GWAS_N_effective_case_control_verified":False,
        "GIGASTROKE_EUR_GRCh37_same_variant_harmonization_verified":False,
        "TPMI_independent_replication_of_GIGASTROKE_claim_approved":False,
        "scientific_status":"HOLD","ready_for_ldsc":False,
        "new_MR_coloc_LDSC_analysis_performed":False,
    }

def download_once(opener=urllib.request.urlopen,expected_bytes=EXPECTED_BYTES):
    with tempfile.TemporaryDirectory(prefix="tpmi43321_lead_") as tmp:
        path=Path(tmp)/"source.gz"
        size=0;h=hashlib.sha256()
        req=urllib.request.Request(URL,headers={"User-Agent":"CKM-Stroke-One-Source-Pinned-GWAS-Lead/1.0"})
        with opener(req,timeout=90) as response:
            u=urlsplit(response.geturl())
            if response.status!=200 or u.scheme!="https" or u.hostname!=HOST:
                raise ValueError("source host or HTTP status changed")
            length=response.headers.get("Content-Length","")
            if length and length.isdecimal() and int(length)!=expected_bytes:
                raise ValueError("source content length differs from previous verified source")
            with path.open("wb") as out:
                while True:
                    data=response.read(1048576)
                    if not data:break
                    size+=len(data)
                    if size>SAFETY_CAP:raise ValueError("source download safety cap exceeded")
                    h.update(data);out.write(data)
        if size!=expected_bytes or h.hexdigest()!=EXPECTED_SHA256:
            raise ValueError("source full compressed file SHA256/length mismatch")
        result=scan_original(path)
        result["downloaded_original_source_bytes"]=size
        result["pinned_source_sha256_verified"]=h.hexdigest()
        result["source_download_count"]=1
        result["temporary_full_file_deleted"]=True
        result["entire_GWAS_summary_statistics_NOT_redistributed"]=True
        return result

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True,type=Path)
    args=ap.parse_args()
    try:
        out=download_once()
    except Exception as e:
        out={"technical_status":"HOLD_ORIGINAL_TPMI_SOURCE_VARIANT_IDENTITY",
             "error_type":type(e).__name__,"error_detail":str(e)[:250],
             "scientific_status":"HOLD","ready_for_ldsc":False,
             "new_MR_coloc_LDSC_analysis_performed":False}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n","utf-8")
    print(json.dumps({"status":out["technical_status"],
                      "source_rows":out.get("all_original_rows_verified"),
                      "lead":out.get("original_significant_associations"),
                      "scientific_status":out["scientific_status"],
                      "error":out.get("error_detail")}))
    if out["technical_status"]!="PASS_PINNED_1739048_SOURCE_ROWS_ONE_SIGNIFICANT_VARIANT":
        raise SystemExit(1)

if __name__=="__main__":
    main()
