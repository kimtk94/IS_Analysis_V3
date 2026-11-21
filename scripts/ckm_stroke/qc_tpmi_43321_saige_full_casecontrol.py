#!/usr/bin/env python3
"""Full public TPMI 433.21 SAIGE GWAS source audit (aggregate output only).

This is a *separate* SAIGE model, not the previous TPMI PLINK2 model.
Stream the 816,992,833 byte original gzip, check decompression CRC and raw
download SHA256/size. Count EVERY complete association's N_case/N_ctrl, valid
alleles/beta/SE/p and coordinate build evidence without releasing row data.
No binary effective-N approval, harmonized allele pairs, or MR/LDSC/coloc.
"""
from __future__ import annotations
import argparse
from collections import Counter
import gzip
import hashlib
import json
import math
from pathlib import Path
import urllib.request
from urllib.parse import urlsplit

URL="https://pheweb.ibms.sinica.edu.tw/download/saige_gwas/433.21"
SOURCE_BYTES=816_992_833
EXPECTED_HEADER=("CHR","POS","MarkerID","Allele1","Allele2","AC_Allele2",
 "AF_Allele2","imputationInfo","BETA","SE","Tstat","var",
 "p.value","p.value.NA","Is.SPA","AF_case","AF_ctrl","N_case","N_ctrl","QC")
REFERENCE_PHENO_CASES=9249
REFERENCE_PHENO_CONTROLS=304660
MAX_GZIP_BYTES=1_000_000_000
MAX_ROWS=50_000_000
MAX_UNPACKED_BYTES=15_000_000_000
MAX_LINE=8192
DNA=set("ACGT")
class BoundedHashingReader:
    def __init__(self,source):
        self.source=source
        self.count=0
        self.sha=hashlib.sha256()
    def read(self,n=-1):
        b=self.source.read(min(n if n>=0 else 1048576,1048576))
        self.count+=len(b)
        if self.count>MAX_GZIP_BYTES:raise ValueError("source compressed size cap exceeded")
        self.sha.update(b)
        return b

def parse_n(value):
    # Never round non-integer N_case / N_ctrl.
    try:
        n=int(value)
        if str(n)==value and n>0:
            return n
    except (ValueError,TypeError,OverflowError):
        pass
    return None

def record_fields(fields, counters, pairs, chr_counts):
    """Validate an original SAIGE association row without exporting effects."""
    if len(fields)!=len(EXPECTED_HEADER):
        raise ValueError("original SAIGE column-count drift")
    chrom,pos,ident,ref,alt=fields[:5]
    chr_counts[chrom]+=1
    if not pos.isdecimal() or int(pos)<=0:raise ValueError("invalid source coordinate")
    if not ref or not alt or ref==alt or any(c not in DNA for c in ref+alt):
        counters["non_SNV_or_invalid_alleles"]+=1
    elif len(ref)!=1 or len(alt)!=1:
        counters["indel_not_SNV"]+=1
    try:
        n_case=parse_n(fields[17])
        n_ctrl=parse_n(fields[18])
        if n_case is None or n_ctrl is None:
            counters["missing_or_noninteger_N_case_ctrl"]+=1
        else:
            pairs[(n_case,n_ctrl)]+=1
            counters["positive_integer_case_control_N"]+=1
            if (n_case,n_ctrl)==(REFERENCE_PHENO_CASES,REFERENCE_PHENO_CONTROLS):
                counters["matches_phenotype_page_N"]+=1
            if n_case+n_ctrl==204099:
                counters["sum_equals_PLINK2_OBS_CT"]+=1
    except (TypeError,ValueError,OverflowError):
        counters["missing_or_noninteger_N_case_ctrl"]+=1
    nums=[]
    for i in (8,9,12):
        try:nums.append(float(fields[i]))
        except (ValueError,OverflowError):nums.append(float("nan"))
    beta,se,p=nums
    if all(math.isfinite(z) for z in nums) and se>0 and 0<=p<=1:
        counters["finite_beta_positive_se_valid_p"]+=1
    else:
        counters["invalid_beta_se_p"]+=1
    qc=fields[19]
    if qc not in ("", "NA",".","PASS","Pass","pass"):
        counters["nontrivial_QC_flag"]+=1
    counters["records"]+=1

def full_audit():
    urlsplit_result=urlsplit(URL)
    if urlsplit_result.scheme!="https" or urlsplit_result.hostname!="pheweb.ibms.sinica.edu.tw":
        raise ValueError("not approved original provider")
    req=urllib.request.Request(URL,headers={
        "User-Agent":"IS_Analysis_V3-SAIGE-full-N_case-N_ctrl-QC/1.0",
        "Accept":"application/octet-stream"})
    counts=Counter();pairs=Counter();chr_counts=Counter()
    uncomp=0
    with urllib.request.urlopen(req,timeout=180) as response:
        final=urlsplit(response.geturl())
        if response.status!=200 or final.scheme!="https" or final.hostname!="pheweb.ibms.sinica.edu.tw":
            raise ValueError("source URL/HTTP drift")
        length=response.headers.get("Content-Length")
        if length is None or int(length)!=SOURCE_BYTES:
            raise ValueError("expected original compressed source length differs")
        reader=BoundedHashingReader(response)
        with gzip.GzipFile(fileobj=reader,mode="rb") as gz:
            headerline=gz.readline(MAX_LINE+1)
            uncomp+=len(headerline)
            if len(headerline)>MAX_LINE:
                raise ValueError("original header exceeds limit")
            header=tuple(headerline.decode("utf-8-sig").rstrip("\r\n").split("\t"))
            if header!=EXPECTED_HEADER:
                raise ValueError("source header schema changed")
            while True:
                line=gz.readline(MAX_LINE+1)
                if not line:break
                if len(line)>MAX_LINE:raise ValueError("source row line over limit")
                uncomp+=len(line)
                if uncomp>MAX_UNPACKED_BYTES:raise ValueError("inflated source exceeds cap")
                if counts["records"]>=MAX_ROWS:raise ValueError("source exceeds row count limit")
                if not line.strip():continue
                fields=line.rstrip(b"\r\n").decode("utf-8").split("\t")
                record_fields(fields,counts,pairs,chr_counts)
        # Consume all source bytes and confirm no trailer was left unread.
        while reader.read(1048576):pass
    if reader.count!=SOURCE_BYTES:raise ValueError("original source compressed bytes incomplete")
    if counts["records"]<100000:raise ValueError("source too few GWAS associations")
    out={"schema":"TPMI43321-SAIGE-full-original-casecontrol-QC-v1",
         "technical_status":"PASS_FULL_GZIP_AND_SOURCE_CASE_CONTROL_SCAN",
         "source":URL,"source_compressed_bytes":reader.count,
         "computed_original_sha256":reader.sha.hexdigest(),
         "sha256_prior_official_reference_approved":False,
         "full_gzip_crc_verified":True,
         "uncompressed_bytes":uncomp,
         "source_header_verified":True,
         "records":counts["records"],
         "by_chr":dict(sorted(chr_counts.items())),
         "QC_counts":dict(counts),
         "distinct_case_control_pairs":len(pairs),
         "most_frequent_case_control_pairs":[
             {"N_case":case,"N_ctrl":ctrl,"count":n,
              "sum":case+ctrl} for (case,ctrl),n in pairs.most_common(25)
         ],
         "all_N_pairs_identical_to_phenotype_page":(
             len(pairs)==1 and pairs.get((REFERENCE_PHENO_CASES,REFERENCE_PHENO_CONTROLS),0)==counts["records"]),
         "missing_QC_pair_suppresses_full_N_approval":(
             counts["positive_integer_case_control_N"]!=counts["records"]),
         "PLINK2_model_OBS_CT":204099,
         "SAIGE_case_control_N_are_NOT_PLINK2_model_case_control_N":True,
         "SAIGE_GRCh38_build_from_source_crosscheck_not_yet_verified":True,
         "genomewide_effects_not_exported":True,
         "SAIGE_GIGASTROKE_harmonization_executed":False,
         "SAIGE_LDSC_neff_approved":False,
         "scientific_status":"HOLD"}
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True,type=Path)
    args=ap.parse_args()
    try:result=full_audit()
    except Exception as exc:
        result={"schema":"TPMI43321-SAIGE-full-original-casecontrol-QC-v1",
                "technical_status":"HOLD_FULL_SOURCE_OR_PARSER",
                "error_type":type(exc).__name__,"error":str(exc)[:360],
                "scientific_status":"HOLD"}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["technical_status"],
        "rows":result.get("records"),
        "sample_pairs":result.get("most_frequent_case_control_pairs"),
        "all_rows_have_phenotype_N":result.get("all_N_pairs_identical_to_phenotype_page"),
        "qc":result.get("QC_counts"),
        "source_sha256":result.get("computed_original_sha256"),
        "error":result.get("error"),"scientific_status":result["scientific_status"]},ensure_ascii=False))
    if result["technical_status"]!="PASS_FULL_GZIP_AND_SOURCE_CASE_CONTROL_SCAN":
        raise SystemExit(1)
if __name__=="__main__":
    main()
