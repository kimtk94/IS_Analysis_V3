#!/usr/bin/env python3
"""Numerically stable SOURCE-P-vs-Wald diagnostic for GLGC zero-P rows.

This audit never turns zero into a conventional floating point P-value;
it compares the original -log10(P) with a two-sided normal approximation
computed in the logarithmic domain from original beta and SE.
No LD, covariate scale, cohort overlap, instrument selection or MR GO.
"""
import argparse
import csv
import gzip
import json
import math
import statistics
from pathlib import Path

from scipy.special import log_ndtr

LOG10=math.log(10.0)
LOG2=math.log(2.0)


def normal_two_sided_neglog10p(z):
    if not math.isfinite(z):
        raise ValueError("finite z required")
    return -(LOG2+float(log_ndtr(-abs(z))))/LOG10


def compare(rows,expected_count=None):
    observed=[]
    keys=set()
    errors=[]
    for i,row in enumerate(rows,1):
        try:
            if float(row["pvalue"])!=0:
                raise ValueError("source P must equal numeric 0")
            beta=float(row["EFFECT_SIZE"])
            se=float(row["SE"])
            logp=float(row["pvalue_neg_log10"])
            af=float(row["POOLED_ALT_AF"])
            sample_n=float(row["N"])
            pos=int(row["POS_b37"])
            chrom=str(row["CHROM"])
            ea=row["ALT"].strip().upper()
            nea=row["REF"].strip().upper()
            if not all(map(math.isfinite,(beta,se,logp,af,sample_n))) or se<=0 or logp<=0 or sample_n<=0 or not 0<=af<=1 or pos<=0 or ea==nea:
                raise ValueError("invalid source association field")
            if len(ea)!=1 or len(nea)!=1 or ea not in "ACGT" or nea not in "ACGT":
                raise ValueError("not a canonical biallelic SNP")
            key=(chrom,pos,tuple(sorted((ea,nea))))
            if key in keys:
                raise ValueError("duplicate variant key")
            keys.add(key)
            z=beta/se
            wald=normal_two_sided_neglog10p(z)
            observed.append({
                "CHROM":chrom,"POS_b37":pos,"rsID":row["rsID"],
                "REF":nea,"ALT":ea,
                "source_negative_log10_p":logp,
                "wald_negative_log10_p":wald,
                "delta_log10_source_minus_wald":logp-wald,
                "z_abs":abs(z),"N":sample_n,"EAF":af
            })
        except (ValueError,TypeError,KeyError,OverflowError) as exc:
            errors.append({"row":i,"error":str(exc)[:200]})
    if expected_count is not None and len(observed)!=expected_count:
        errors.append({"row":"total","error":"unexpected number of original source zero-P SNPs"})
    diffs=sorted(abs(x["delta_log10_source_minus_wald"]) for x in observed)
    quantile=lambda q: diffs[min(len(diffs)-1,int((len(diffs)-1)*q))] if diffs else None
    return ({
        "status":"PASS" if not errors else "HOLD",
        "scientific_status":"HOLD",
        "scope":"GLGC_CHR19_LDL_SOURCE_PZERO_WALD_LOGDOMAIN_DIAGNOSTIC_ONLY",
        "number_original_pzero_variants":len(observed),
        "source_vs_wald_abs_log10_delta_median":statistics.median(diffs) if diffs else None,
        "source_vs_wald_abs_log10_delta_p95":quantile(0.95),
        "source_vs_wald_abs_log10_delta_max":max(diffs) if diffs else None,
        "abs_delta_gt_1":sum(v>1 for v in diffs),
        "abs_delta_gt_5":sum(v>5 for v in diffs),
        "abs_delta_gt_50":sum(v>50 for v in diffs),
        "source_negative_log10_p_range":[min(x["source_negative_log10_p"] for x in observed),max(x["source_negative_log10_p"] for x in observed)] if observed else [],
        "wald_negative_log10_p_range":[min(x["wald_negative_log10_p"] for x in observed),max(x["wald_negative_log10_p"] for x in observed)] if observed else [],
        "errors":errors,
        "source_file_whole_genome_hash_verified":False,
        "causal_or_mr_readiness":"NOT_ASSESSED",
        "interpretation":"The discrepancy is descriptive, not an independently calibrated significance or biology test. Beta/SE normal approximation may differ from published meta-analysis exact P and phenotype transformation; manual review and harmonization still required."
    },observed)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input",required=True)
    p.add_argument("--json-output",required=True)
    p.add_argument("--tsv-output",required=True)
    p.add_argument("--expected-count",type=int,default=332)
    args=p.parse_args()
    with gzip.open(args.input,"rt",encoding="utf-8") as f:
        input_rows=list(csv.DictReader(f,delimiter="\t"))
    summary,records=compare(input_rows,args.expected_count)
    out=Path(args.json_output)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    if summary["status"]=="PASS":
        dest=Path(args.tsv_output)
        dest.parent.mkdir(parents=True,exist_ok=True)
        with gzip.open(dest,"wt",encoding="utf-8",newline="") as result:
            cols=list(records[0].keys())
            writer=csv.DictWriter(result,fieldnames=cols,delimiter="\t")
            writer.writeheader()
            writer.writerows(records)
    print(json.dumps(summary,ensure_ascii=False))
    return int(summary["status"]!="PASS")


if __name__=="__main__":
    raise SystemExit(main())
