#!/usr/bin/env python3
"""One pinned TPMI 433.21 PLINK2 gzip: full-genome row QC, never GWAS export.

Only official 58,302,393 byte public source; no retries or downloads of 817MB
SAIGE file; original compressed SHA256 MUST equal prior independently observed
hash. Validate ALL decompressed source rows and report aggregate-only JSON.
This does not demonstrate case-level independence or validate effective N.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import json
import math
from pathlib import Path
import tempfile
import urllib.request
from urllib.parse import urlsplit

from scripts.ckm_stroke.inspect_tpmi_43321_original_sumstats_once import (
    URL, HOST, EXPECTED_BYTES, SAFETY_CAP
)

EXPECTED_SHA256="ce14daef23f9fadf921030d9ed2429383d927f6da391c620d9abf870045b2e62"
EXPECTED_HEADER=[
    "#CHROM","POS","ID","REF","ALT","PROVISIONAL_REF?","A1",
    "OMITTED","A1_FREQ","TEST","OBS_CT","OR","LOG(OR)_SE",
    "Z_STAT","P","ERRCODE","rsID"
]
CHROMS={str(n) for n in range(1,23)}
VALID_DNA=set("ACGT")
MAX_UNCOMPRESSED_BYTES=8_000_000_000
MAX_VARIANTS=25_000_000
LINE_LIMIT=100_000
MISSING={"NA","N/A",".","","nan","NaN","-nan"}

def parse_float(value):
    if value in MISSING:
        return None
    try:
        x=float(value)
    except (ValueError,OverflowError):
        return None
    return x if math.isfinite(x) else None

def original_single_get(opener=urllib.request.urlopen,expected_size=EXPECTED_BYTES):
    target=urlsplit(URL)
    if target.scheme!="https" or target.hostname!=HOST:
        raise ValueError("official GWAS URL host not approved")
    request=urllib.request.Request(URL,headers={
        "User-Agent":"CKM-Stroke-TPMI-43321-Full-Original-QC-One-Pass/1.0",
        "Accept":"application/gzip,application/octet-stream"})
    with tempfile.TemporaryDirectory(prefix="ckm_tpmi_full_qc_") as td:
        path=Path(td)/"source_original_ephemeral.gz"
        h=hashlib.sha256()
        size=0
        with opener(request,timeout=90) as response:
            actual=urlsplit(response.geturl())
            if response.status!=200 or actual.scheme!="https" or actual.hostname!=HOST:
                raise ValueError("unexpected source response or off-site redirect")
            length=response.headers.get("Content-Length","")
            if length and length.isdecimal() and int(length)!=expected_size:
                raise ValueError("official source content length differs from pinned HEAD")
            with path.open("wb") as out:
                while True:
                    chunk=response.read(1<<20)
                    if not chunk:
                        break
                    size+=len(chunk)
                    if size>SAFETY_CAP:
                        raise ValueError("source download size cap exceeded")
                    out.write(chunk)
                    h.update(chunk)
        if size!=expected_size:
            raise ValueError("original compressed bytes incomplete")
        sha=h.hexdigest()
        if sha!=EXPECTED_SHA256:
            raise ValueError("original source content SHA256 drift; DO NOT analyze changed file")
        result=scan_all_rows(path)
        result["full_original_compressed_sha256"]=sha
        result["pinned_prior_original_sha256_match"]=True
        result["original_compressed_bytes"]=size
        result["source_url"]=URL
        result["original_full_file_deleted_after_scan"]=True
        result["original_source_download_attempts"]=1
        result["source_file_or_variant_effects_uploaded"]=False
        return result

def round_bucket(n):
    if n<0:
        return "<0"
    if n<1000:return "<1k"
    if n<10000:return "1k-10k"
    if n<100000:return "10k-100k"
    if n<300000:return "100k-300k"
    if n<313909:return "300k-313909"
    if n==313909:return "=313909"
    return ">313909"

def scan_all_rows(path):
    if not path.is_file():
        raise ValueError("full original compressed source absent")
    stats=Counter()
    bychr=Counter()
    testtypes=Counter()
    errcodes=Counter()
    obs_buckets=Counter()
    freq_buckets=Counter()
    repeated_chrom_pos=0
    positions_out_of_order=0
    seen_variant_keys=set()  # 128-bit hash of exact source variant+test; collision negligible, not mathematical proof
    extrema={"min_OBS_CT":None,"max_OBS_CT":None,
             "min_A1_FREQ":None,"max_A1_FREQ":None,
             "min_OR":None,"max_OR":None}
    def ext(key,value,fn):
        current=extrema[key]
        if current is None or fn(value,current):
            extrema[key]=value
    previous_position=None
    prev_chrom=None
    uncompressed_bytes=0
    with gzip.open(path,"rt",encoding="utf-8-sig",errors="strict",newline="") as f:
        header=f.readline()
        uncompressed_bytes+=len(header.encode("utf-8"))
        if header.rstrip("\r\n").split("\t")!=EXPECTED_HEADER:
            raise ValueError("original TPMI 433.21 PLINK2 full source header differs from pinned first-run header")
        for number,line in enumerate(f,start=2):
            if len(line)>LINE_LIMIT:
                raise ValueError("oversized source TSV association row")
            uncompressed_bytes+=len(line.encode("utf-8"))
            if uncompressed_bytes>MAX_UNCOMPRESSED_BYTES:
                raise ValueError("uncompressed source byte cap exceeded")
            if not line.strip():
                stats["empty_lines"]+=1
                continue
            if stats["records"]>=MAX_VARIANTS:
                raise ValueError("original source record count cap exceeded")
            col=line.rstrip("\r\n").split("\t")
            if len(col)!=len(EXPECTED_HEADER):
                raise ValueError("source original association TSV truncated/extra cells on line "+str(number))
            stats["records"]+=1
            chrom,pos,ident,ref,alt,provisional,a1,omitted,af,model,obs,odds,se,z,p,err,rsid=col
            bychr[chrom]+=1
            testtypes[model]+=1
            errcodes[err]+=1
            if chrom not in CHROMS:
                stats["non_autosomal_rows"]+=1
            try:
                bp=int(pos)
                if bp<=0:raise ValueError("invalid")
            except ValueError:
                raise ValueError("nonpositive or nonnumeric GRCh38 position at source line "+str(number))
            if chrom==prev_chrom:
                if bp<previous_position:
                    positions_out_of_order+=1
                if bp==previous_position:
                    repeated_chrom_pos+=1
            prev_chrom=chrom
            previous_position=bp
            if not ident or not ref or not alt or not a1:
                stats["missing_identifier_or_allele"]+=1
            else:
                alleles_valid=all(c in VALID_DNA for c in ref+alt) and a1 in (ref,alt)
                if not alleles_valid:
                    stats["A1_not_REF_ALT_or_invalid_DNA"]+=1
                if ref==alt:stats["REF_equal_ALT"]+=1
            if model!="ADD":
                stats["non_ADD_test"]+=1
            identhash=hashlib.blake2b(
                ("\t".join([chrom,pos,ref,alt,a1,model])).encode("ascii","backslashreplace"),
                digest_size=16).digest()
            if identhash in seen_variant_keys:
                stats["repeated_variant_test_128bit_hash"]+=1
            else:
                seen_variant_keys.add(identhash)
            freq=parse_float(af)
            if freq is None or not 0<=freq<=1:
                stats["invalid_or_missing_A1_FREQ"]+=1
            else:
                stats["freq_numeric"]+=1
                ext("min_A1_FREQ",freq,lambda a,b:a<b)
                ext("max_A1_FREQ",freq,lambda a,b:a>b)
                maf=min(freq,1-freq)
                if maf<0.0001:freq_buckets["MAF_<0.0001"]+=1
                elif maf<0.001:freq_buckets["MAF_0.0001_0.001"]+=1
                elif maf<0.01:freq_buckets["MAF_0.001_0.01"]+=1
                elif maf<0.05:freq_buckets["MAF_0.01_0.05"]+=1
                else:freq_buckets["MAF_ge_0.05"]+=1
            if obs.isdecimal():
                n=int(obs)
                if n>0:
                    stats["OBS_CT_valid"]+=1
                    obs_buckets[round_bucket(n)]+=1
                    ext("min_OBS_CT",n,lambda a,b:a<b)
                    ext("max_OBS_CT",n,lambda a,b:a>b)
                    if n>313909:stats["OBS_CT_exceeds_phenotype_page_total"]+=1
                else:
                    stats["OBS_CT_nonpositive_or_missing"]+=1
            else:
                stats["OBS_CT_nonpositive_or_missing"]+=1
            val_p=parse_float(p)
            val_or=parse_float(odds)
            val_se=parse_float(se)
            val_z=parse_float(z)
            valid_p=val_p is not None and 0<=val_p<=1
            if valid_p:
                stats["P_valid"]+=1
                if val_p==0:stats["P_reported_zero"]+=1
                if val_p<=5e-8:stats["P_le_5e_8"]+=1
            else:stats["P_invalid_or_missing"]+=1
            if val_or is None or val_or<=0:stats["OR_invalid_or_missing"]+=1
            else:
                stats["OR_positive"]+=1
                ext("min_OR",val_or,lambda a,b:a<b)
                ext("max_OR",val_or,lambda a,b:a>b)
            if val_se is None or val_se<=0:stats["LOGOR_SE_invalid_or_missing"]+=1
            else:stats["LOGOR_SE_positive"]+=1
            if val_z is None:stats["Z_STAT_invalid_or_missing"]+=1
            else:stats["Z_STAT_finite"]+=1
            if (valid_p and val_or is not None and val_or>0
                and val_se is not None and val_se>0 and val_z is not None
                and freq is not None and 0<=freq<=1 and obs.isdecimal() and int(obs)>0
                and a1 in (ref,alt) and model=="ADD" and err in ("NA",".","")):
                stats["fully_numeric_ADD_no_error"]+=1
                beta=math.log(val_or)
                inferred_z=beta/val_se
                if abs(inferred_z-val_z)>=0.05:
                    stats["abs_OR_SE_z_diff_ge_0.05"]+=1
                # A large P-v-Wald discrepancy may reflect logistic Firth/SPA or printing precision:
                # inventory only, never silently "correct" original P values.
                if val_p>0 and abs(val_z)<37:
                    estimated=math.erfc(abs(val_z)/math.sqrt(2))
                    if estimated>0:
                        stats["P_vs_Z_comparable"]+=1
                        drift=abs(math.log10(val_p)-math.log10(estimated))
                        if drift>=0.1:stats["P_vs_Z_minuslog10_diff_ge_0.1"]+=1
                        if drift>=1:stats["P_vs_Z_minuslog10_diff_ge_1"]+=1
            # avoid outputting raw per-SNP effects or identifiers
    if stats["records"]==0 or uncompressed_bytes==0 or stats["records"]>MAX_VARIANTS:
        raise ValueError("missing genomewide original association rows")
    return {
        "schema":"TPMI43321-original-gwas-complete-row-QC-v1",
        "technical_status":"PASS_FULL_SOURCE_SHA256_AND_ALL_ROWS_PARSED_WITH_QC_FLAGS",
        "source":"TPMI_433.21_PheWeb_GRCh38_logistic_PLINK2",
        "pinned_original_columns":EXPECTED_HEADER,
        "rows_total_in_original":stats["records"],
        "uncompressed_source_bytes":uncompressed_bytes,
        "by_chr_all_original_rows":dict(sorted(bychr.items(),key=lambda x:(not x[0].isdigit(),int(x[0]) if x[0].isdigit() else 0))),
        "original_additive_test_categories":dict(testtypes.most_common(20)),
        "original_ERRORCODE_categories":dict(errcodes.most_common(20)),
        "row_level_QC_counts":dict(stats),
        "observed_sample_N_bins":dict(obs_buckets),
        "original_MAF_bins_from_A1_FREQ":dict(freq_buckets),
        "value_min_max":extrema,
        "source_rows_out_of_coordinate_order":positions_out_of_order,
        "adjacent_identical_chromosome_position_rows":repeated_chrom_pos,
        "duplicate_variant_test_128bit_hashes":stats["repeated_variant_test_128bit_hash"],
        "duplicate_hash_is_not_mathematical_proof_of_duplicate_identity":True,
        "OR_to_logOR_conversion_needed_before_MR_coloc":True,
        "A1_is_tested_allele_not_assume_ALT":True,
        "P_vs_Wald_differences_are_descriptive_not_genomic_control":True,
        "OBS_CT_is_not_SNP_case_control_effective_N":True,
        "source_cumulative_publisher_expected_full_SHA_not_available":True,
        "GIGASTROKE_and_TPMI_individual_overlap_resolved":False,
        "TPMI_2026_study_phenotype_equivalence_resolved":False,
        "all_SNPs_eligible_for_harmonization":False,
        "downstream_recommended_gate":"CHECK_QC_COUNTS_AND_LICENSE_BEFORE_FILTERED_VCF_HARMONIZATION",
        "scientific_status":"HOLD",
        "ready_for_ldsc":False,
        "MR_coloc_LDSC_executed":False,
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    arg=ap.parse_args()
    try:
        out=original_single_get()
    except Exception as e:
        out={"technical_status":"HOLD_FULL_TPMI_SOURCE_OR_QC_FAILED",
             "error_type":type(e).__name__,"error_detail":str(e)[:300],
             "scientific_status":"HOLD","ready_for_ldsc":False,
             "MR_coloc_LDSC_executed":False}
    arg.out.parent.mkdir(parents=True,exist_ok=True)
    arg.out.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":out["technical_status"],"rows":out.get("rows_total_in_original"),
                      "fully_numeric_ADD_no_error":out.get("row_level_QC_counts",{}).get("fully_numeric_ADD_no_error"),
                      "source_SHA_verified":out.get("pinned_prior_original_sha256_match",False),
                      "source_QC_errors":{k:v for k,v in out.get("row_level_QC_counts",{}).items()
                                       if "invalid" in k or "missing" in k or "diff_ge" in k},
                      "error":out.get("error_detail"),"scientific_HOLD":out["scientific_status"]}))
    if out["technical_status"]!="PASS_FULL_SOURCE_SHA256_AND_ALL_ROWS_PARSED_WITH_QC_FLAGS":
        raise SystemExit(1)

if __name__=="__main__":
    main()
