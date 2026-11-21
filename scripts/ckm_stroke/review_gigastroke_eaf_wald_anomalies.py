#!/usr/bin/env python3
"""Audit one pinned, whole-source EAF/Wald result artifact without re-downloading GWAS.

Strictly revalidates source publisher MD5, Locke/Yengo allele orders and 22
chromosome frequency/P-vs-Wald counts via upstream independent reviewer.
Only aggregate categories enter the report. Different phenotypes having
outliers on the same chromosome does NOT establish the same variant.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

from scripts.ckm_stroke.audit_gigastroke_eaf_gc import review_diagnostics

TRAITS=("AIS","LAA","CE","SVD")
P_CATEGORIES=("reported_p_is_zero","abs_log10p_delta_lt_001",
              "abs_log10p_delta_001_to_010","abs_log10p_delta_010_to_1",
              "abs_log10p_delta_ge_1")
HRS_CATEGORIES=("hrs_absdiff_lt_002","hrs_absdiff_002_to_005",
                "hrs_absdiff_005_to_010","hrs_absdiff_010_to_020",
                "hrs_absdiff_ge_020","hrs_reference_not_available")
MAF_CATEGORIES=("maf_lt_001","maf_001_to_005","maf_ge_005")

def percentage(n,den):
    if type(n) is not int or type(den) is not int or den<=0 or n<0 or n>den:
        raise ValueError("invalid source evidence percentage")
    return round(100*n/den,6)

def analyze(original,checked):
    expected=review_diagnostics(original)
    if expected!=checked:
        raise ValueError("independent submitted EAF/Wald review disagrees with original publisher-MD5 source evidence")
    rows=checked.get("outcomes",[])
    if len(rows)!=4 or [x.get("phenotype") for x in rows]!=list(TRAITS):
        raise ValueError("wrong complete four subtype order")
    profiles=[]
    excessive_hrs_chrs=[]
    p_moderate_chrs=[]
    for item in rows:
        phenotype=item["phenotype"]
        n=item["matched_valid_SNPs"]
        p=item["P_vs_BETA_SE_Wald"]
        hrs=item["HRS_reference_freq_absdiff"]
        maf=item["MAF_distribution"]
        if any(type(v) is not int or v<0 for c in (p,hrs,maf) for v in c.values()):
            raise ValueError("noninteger original source summary counts")
        if (sum(p.values())!=n or sum(hrs.values())!=n or sum(maf.values())!=n):
            raise ValueError("source phenotype histogram does not reconcile")
        extreme_hrs=[]
        p_delta_mid=[]
        for chromosome in item["by_chr"]:
            ch=chromosome["chr"]
            if type(ch) is not int or not 1<=ch<=22:
                raise ValueError("bad source chr")
            if chromosome["HRS_absdiff_ge_020"]:
                extreme_hrs.append({"chr":ch,"n":chromosome["HRS_absdiff_ge_020"]})
            # Full original evidence has more detailed P bins than the compact review.
            src=next(x for x in original["outcomes"] if x["phenotype"]==phenotype)
            full_chr=src["eaf_gc_source_descriptive"]["by_chr"][ch-1]
            mid=full_chr["categories"].get("abs_log10p_delta_010_to_1",0)
            if type(mid) is not int or mid<0:
                raise ValueError("invalid full-source P/Wald bin")
            if mid:
                p_delta_mid.append({"chr":ch,"n":mid})
        if sum(x["n"] for x in extreme_hrs)!=hrs["hrs_absdiff_ge_020"]:
            raise ValueError("extreme HRS frequency count not reconciled per chromosome")
        if sum(x["n"] for x in p_delta_mid)!=p["abs_log10p_delta_010_to_1"]:
            raise ValueError("mid P/Wald discrepancy count not reconciled per chromosome")
        excessive_hrs_chrs.append({x["chr"] for x in extreme_hrs})
        p_moderate_chrs.append({x["chr"] for x in p_delta_mid})
        profiles.append({
           "phenotype":phenotype,
           "valid_source_pre_munge_SNPs":n,
           "source_MAF_below_005":maf["maf_lt_001"]+maf["maf_001_to_005"],
           "source_MAF_below_005_pct":percentage(maf["maf_lt_001"]+maf["maf_001_to_005"],n),
           "HRS_absdiff_ge_005":hrs["hrs_absdiff_005_to_010"]+
                                  hrs["hrs_absdiff_010_to_020"]+hrs["hrs_absdiff_ge_020"],
           "HRS_absdiff_ge_010":hrs["hrs_absdiff_010_to_020"]+hrs["hrs_absdiff_ge_020"],
           "HRS_absdiff_ge_020":hrs["hrs_absdiff_ge_020"],
           "HRS_absdiff_ge_020_chr":extreme_hrs,
           "P_vs_Wald_delta_ge_001":p["abs_log10p_delta_001_to_010"]+
                                    p["abs_log10p_delta_010_to_1"]+p["abs_log10p_delta_ge_1"],
           "P_vs_Wald_delta_ge_010":p["abs_log10p_delta_010_to_1"]+p["abs_log10p_delta_ge_1"],
           "P_vs_Wald_delta_ge_1":p["abs_log10p_delta_ge_1"],
           "P_vs_Wald_delta_010_to_1_chromosomes":p_delta_mid,
           "source_P_is_zero":p["reported_p_is_zero"],
           "source_BETA_not_exact_logOR":
               item["BETA_vs_logOR"]["beta_logOR_absdelta_1e10_to_1e6"]+
               item["BETA_vs_logOR"]["beta_logOR_absdelta_ge_1e6"],
        })
    common_hrs=sorted(set.intersection(*excessive_hrs_chrs))
    common_p=sorted(set.intersection(*p_moderate_chrs))
    all_original=len({x["accession"] for x in profiles if "accession" in x})==4
    return {
        "schema":"ckm-gigastroke-four-original-EAF-Wald-secondary-aggregate-review-v1",
        "technical_status":"PASS_MD5_PINNED_ANOMALY_LOCUS_CHROMOSOME_AGGREGATE_ONLY",
        "base_source_technical_status":checked["technical_status"],
        "source_outcome_frequencies_are_not_GIANT_HRS_frequencies":True,
        "per_variant_shared_identity_across_subtypes_established":False,
        "GC_correction_existence_or_factor_inferred":False,
        "original_effect_data_or_raw_SNP_ids_exported":False,
        "same_chromosome_is_not_same_SNP":True,
        "extreme_HRS_reference_delta_shared_chromosomes":common_hrs,
        "P_vs_Wald_mid_discrepancy_shared_chromosomes_all_four":common_p,
        "outcomes":profiles,
        "next_manual_review":[
          "Determine whether chr11 HRS reference-frequency outlier is same SNP in all four phenotypes; aggregate counts cannot tell",
          "Inspect chr4 P/Wald 0.1-1 log10 differences, especially AIS and CE, against source precision and published methods",
          "Verify original public GWAS GC and cohort-specific N_eff provenance before any LDSC regression"
        ],
        "scientific_status":"HOLD","ready_for_ldsc":False,
        "outcome_N_eff_approved":False,
        "LDSC_h2_intercept_rg_MR_coloc_executed":False,
    }

def to_markdown(report):
    lines=[
      "# GIGASTROKE original EUR 4-subtype EAF / P-vs-Wald anomaly review",
      "",
      "Technical gate: "+report["technical_status"],
      "Scientific status: HOLD. This review does NOT perform LDSC, correct source summary statistics, or approve meta-analysis N_eff.",
      "",
      "| Outcome | Full original pre-munge SNPs | MAF <0.05 | HRS difference >=0.05 | HRS difference >=0.20 | P/Wald delta >=0.1 log10 | P/Wald delta >=1 log10 |",
      "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in report["outcomes"]:
        lines.append("| {phenotype} | {valid_source_pre_munge_SNPs:,} | {source_MAF_below_005:,} | {HRS_absdiff_ge_005:,} | {HRS_absdiff_ge_020:,} | {P_vs_Wald_delta_ge_010:,} | {P_vs_Wald_delta_ge_1:,} |".format(**r))
    lines.extend([
       "",
       "Extreme HRS reference-frequency difference shared at the chromosome level: "+
         (", ".join("chr"+str(x) for x in report["extreme_HRS_reference_delta_shared_chromosomes"]) or "none"),
       "",
       "**Important:** Shared chromosome does not imply the same rsID or variant. Population differences do not prove error.",
       "",
       "P-vs-Wald differences may reflect rounding, model-specific statistics or processing. They do NOT establish a genomic-control lambda.",
       "",
       "The frequency reference is GIANT's HRS tested-allele frequency, not stroke-cohort allele frequency. No association effects or SNP identifiers were published.",
       "",
       "All comparisons are limited to the pre-munge Yengo2018/HM3/EUR-LDscore intersected candidate panel. Actual final LDSC regression SNP count remains unknown.",
       "",
       "## Per-trait chromosome-level review",
    ])
    for p in report["outcomes"]:
        lines.append("### "+p["phenotype"])
        lines.append("Extreme HRS discrepancy chromosome counts: "+json.dumps(p["HRS_absdiff_ge_020_chr"],ensure_ascii=False))
        lines.append("Moderate P/Wald discrepancy chromosome counts: "+json.dumps(p["P_vs_Wald_delta_010_to_1_chromosomes"],ensure_ascii=False))
        lines.append("")
    return "\n".join(lines)+"\n"

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--full-source-evidence",required=True,type=Path)
    ap.add_argument("--independent-review",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    ap.add_argument("--markdown",required=True,type=Path)
    a=ap.parse_args()
    try:
        source=a.full_source_evidence.read_bytes()
        checked=a.independent_review.read_bytes()
        report=analyze(json.loads(source),json.loads(checked))
        report["full_original_aggregate_input_sha256"]=hashlib.sha256(source).hexdigest()
        report["original_eaf_wald_review_input_sha256"]=hashlib.sha256(checked).hexdigest()
    except Exception as exc:
        report={"technical_status":"HOLD_EAF_WALD_ANOMALY_REVIEW_INVALID",
                "error_type":type(exc).__name__,"error_detail":str(exc)[:300],
                "scientific_status":"HOLD","ready_for_ldsc":False}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    if report.get("technical_status")=="PASS_MD5_PINNED_ANOMALY_LOCUS_CHROMOSOME_AGGREGATE_ONLY":
        a.markdown.parent.mkdir(parents=True,exist_ok=True)
        a.markdown.write_text(to_markdown(report),encoding="utf-8")
    print(json.dumps({"status":report["technical_status"],
                      "shared_HRS_extreme_chr":report.get("extreme_HRS_reference_delta_shared_chromosomes"),
                      "traits":{x["phenotype"]:{"maf_below_005":x["source_MAF_below_005"],
                            "HRS_delta_ge_020":x["HRS_absdiff_ge_020"],
                            "P_vs_Wald_delta_ge_010":x["P_vs_Wald_delta_ge_010"]}
                            for x in report.get("outcomes",[])},
                      "science":report["scientific_status"],
                      "error":report.get("error_detail")},ensure_ascii=False))
    if report["technical_status"]!="PASS_MD5_PINNED_ANOMALY_LOCUS_CHROMOSOME_AGGREGATE_ONLY":
        raise SystemExit(1)

if __name__=="__main__":
    main()
