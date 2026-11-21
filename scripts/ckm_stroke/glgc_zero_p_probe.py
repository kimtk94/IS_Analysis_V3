#!/usr/bin/env python3
"""Bounded original GLGC pvalue=0 rescue-candidate audit.

No causal inference, no MR-ready exports, no invented p values.
The output only records observed original negative-log10 P values for
canonical A/C/G/T SNP rows with numeric pvalue exactly zero.
"""
import argparse
import csv
import datetime as dt
import gzip
import json
import math
import os
from pathlib import Path

from scripts.ckm_stroke.glgc_pysam_probe import fetch_official_header
from scripts.ckm_stroke.glgc_region_probe import BASE, FILE_NAMES, REQUIRED

DNA={"A","C","G","T"}
COLS=["CHROM","POS_b37","rsID","REF","ALT","N","POOLED_ALT_AF",
      "EFFECT_SIZE","SE","pvalue","pvalue_neg_log10","z_abs"]


def audit_pzero(trait, chrom, factory=None, header_fetcher=None, max_rows=5000000):
    chrom=str(chrom)
    if trait not in FILE_NAMES or chrom not in {str(i) for i in range(1,23)}:
        raise ValueError("unsupported trait or chromosome")
    if type(max_rows) is not int or max_rows<100:
        raise ValueError("invalid full-source maximum")
    url=BASE+FILE_NAMES[trait]
    observed=(header_fetcher or fetch_official_header)(url)
    if observed.get("status")!="PASS":
        return {"status":"HOLD","reason":"official source header inaccessible",
                "header_evidence":observed,"scientific_status":"HOLD"}, []
    header=observed["header"].split("\t")
    required=REQUIRED | {"N","POOLED_ALT_AF","pvalue_neg_log10","rsID"}
    if not required.issubset(set(header)) or len(set(header))!=len(header):
        return {"status":"HOLD","reason":"required original GLGC columns missing/duplicate",
                "scientific_status":"HOLD"}, []
    ix={x:header.index(x) for x in required}
    if factory is None:
        import pysam
        factory=pysam.TabixFile
    result_rows=[]
    source_rows=0
    eligible_snp=0
    invalid_pzero_metadata=0
    overflow=False
    reported_logps=[]
    other_errors=0
    seen_zero_p_keys=set()
    duplicate_zero_p_keys=0
    with factory(url) as source:
        for raw in source.fetch(chrom):
            source_rows+=1
            if source_rows>max_rows:
                overflow=True
                break
            fields=raw.rstrip("\n").split("\t")
            if len(fields)!=len(header):
                other_errors+=1
                continue
            try:
                source_chr=fields[ix["CHROM"]].strip().removeprefix("chr")
                pos=int(fields[ix["POS_b37"]])
                nea=fields[ix["REF"]].strip().upper()
                ea=fields[ix["ALT"]].strip().upper()
            except ValueError:
                other_errors+=1
                continue
            if (source_chr!=chrom or pos<=0 or
                len(ea)!=1 or len(nea)!=1 or
                ea not in DNA or nea not in DNA or ea==nea):
                continue
            try:
                beta=float(fields[ix["EFFECT_SIZE"]])
                se=float(fields[ix["SE"]])
                p=float(fields[ix["pvalue"]])
            except ValueError:
                continue
            if not all(map(math.isfinite,(beta,se,p))) or se<=0:
                continue
            eligible_snp+=1
            if p!=0:
                continue
            key=(chrom,pos,tuple(sorted((ea,nea))))
            if key in seen_zero_p_keys:
                duplicate_zero_p_keys+=1
            seen_zero_p_keys.add(key)
            raw_n=fields[ix["N"]].strip()
            raw_af=fields[ix["POOLED_ALT_AF"]].strip()
            raw_logp=fields[ix["pvalue_neg_log10"]].strip()
            valid=True
            try:
                logp=float(raw_logp)
                n=float(raw_n)
                af=float(raw_af)
                if (not all(map(math.isfinite,(logp,n,af))) or
                    logp<=0 or n<=0 or af<0 or af>1):
                    valid=False
            except ValueError:
                valid=False
            if not valid:
                invalid_pzero_metadata+=1
            else:
                reported_logps.append(logp)
            result_rows.append({
                "CHROM":chrom, "POS_b37":pos, "rsID":fields[ix["rsID"]].strip(),
                "REF":nea,"ALT":ea,"N":raw_n,"POOLED_ALT_AF":raw_af,
                "EFFECT_SIZE":fields[ix["EFFECT_SIZE"]].strip(),
                "SE":fields[ix["SE"]].strip(),
                "pvalue":fields[ix["pvalue"]].strip(),
                "pvalue_neg_log10":raw_logp,
                "z_abs":abs(beta/se),
            })
    status="PASS" if not overflow and not invalid_pzero_metadata and not other_errors and not duplicate_zero_p_keys else "HOLD"
    return {
        "status":status, "scientific_status":"HOLD",
        "scope":"ONE_ORIGINAL_GLGCEUR_CHROMOSOME_ZERO_P_SNP_EVIDENCE_NOT_MR",
        "source":"GLGC_Graham2021","ancestry":"EUR","source_url":url,
        "chromosome":chrom,"trait":trait,
        "source_header_status":observed["status"],
        "source_rows":source_rows,"eligible_snp_rows_including_pzero":eligible_snp,
        "source_pzero_canonical_snp_count":len(result_rows),
        "pzero_with_valid_original_logp":len(reported_logps),
        "pzero_with_invalid_auxiliary_metadata":invalid_pzero_metadata,
        "source_parse_errors":other_errors,
        "duplicate_zero_p_variant_keys":duplicate_zero_p_keys,
        "row_cap_hit":overflow,
        "original_logp_min":min(reported_logps) if reported_logps else None,
        "original_logp_max":max(reported_logps) if reported_logps else None,
        "original_logp_ge_300":sum(v>=300 for v in reported_logps),
        "source_full_genome_checksum_verified":False,
        "instrument_significance_or_harmonization":"NOT_RUN",
        "note":"P==0 is an observed source numeric encoding, not proof of non-association. Original -log10P must be validated separately against beta/SE and source transformation before any downstream association inference."
    },result_rows


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--trait",default="LDL_C")
    p.add_argument("--chromosome",default="19")
    p.add_argument("--max-rows",type=int,default=5000000)
    p.add_argument("--expected-original-rows",type=int,default=None)
    p.add_argument("--expected-pzero-rows",type=int,default=None)
    p.add_argument("--json-output",required=True)
    p.add_argument("--tsv-output",required=True)
    a=p.parse_args()
    summary,records=audit_pzero(a.trait,a.chromosome,max_rows=a.max_rows)
    if a.expected_original_rows is not None and summary.get("source_rows")!=a.expected_original_rows:
        summary["status"]="HOLD"
        summary["expected_original_rows_mismatch"]=a.expected_original_rows
    if a.expected_pzero_rows is not None and summary.get("source_pzero_canonical_snp_count")!=a.expected_pzero_rows:
        summary["status"]="HOLD"
        summary["expected_pzero_rows_mismatch"]=a.expected_pzero_rows
    summary["timestamp_utc"]=dt.datetime.now(dt.timezone.utc).isoformat()
    summary["source_audit_commit"]=os.getenv("GITHUB_SHA","")
    output=Path(a.json_output)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    if summary["status"]=="PASS":
        tsv=Path(a.tsv_output)
        tsv.parent.mkdir(parents=True,exist_ok=True)
        with gzip.open(tsv,"wt",encoding="utf-8",newline="") as out:
            writer=csv.DictWriter(out,fieldnames=COLS,delimiter="\t")
            writer.writeheader()
            writer.writerows(records)
    print(json.dumps(summary,ensure_ascii=False))
    return int(summary["status"]!="PASS")


if __name__=="__main__":
    raise SystemExit(main())
