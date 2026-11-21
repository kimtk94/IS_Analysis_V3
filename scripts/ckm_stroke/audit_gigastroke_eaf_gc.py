#!/usr/bin/env python3
"""Original-source EAF and P-vs-Wald diagnostic bins, not GC correction.

One source GWAS stream; only aggregate bins are emitted. GIANT's HRS
tested-allele reference frequency is NOT an equivalent cohort-specific
GIGASTROKE effect-allele frequency. Report descriptive differences only.
P-vs-BETA/SE differences cannot uniquely identify genomic control, and the
original 4 outcome GWAS lack variant-specific N_eff.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path

from scripts.ckm_stroke.audit_bmi_duplicate_zero_p import minus_log10_two_sided_p_from_z
from scripts.ckm_stroke.audit_gigastroke_independent_full_direction import direction_counts
from scripts.ckm_stroke.review_gigastroke_full_evidence import PINNED

MAF_KEYS=("maf_lt_001","maf_001_to_005","maf_ge_005")
HRS_KEYS=("hrs_absdiff_lt_002","hrs_absdiff_002_to_005",
          "hrs_absdiff_005_to_010","hrs_absdiff_010_to_020",
          "hrs_absdiff_ge_020","hrs_reference_not_available")
WALD_KEYS=("reported_p_is_zero","abs_log10p_delta_lt_001",
           "abs_log10p_delta_001_to_010","abs_log10p_delta_010_to_1",
           "abs_log10p_delta_ge_1")
LOGOR_KEYS=("beta_logOR_absdelta_lt_1e10","beta_logOR_absdelta_1e10_to_1e6",
            "beta_logOR_absdelta_ge_1e6","source_OR_nonpositive_or_nonfinite")
EAF_KEYS=("raw_eaf_boundary_zero_or_one","raw_eaf_strict_between_zero_and_one")
DIAGNOSTIC_KEYS=MAF_KEYS+HRS_KEYS+WALD_KEYS+LOGOR_KEYS+EAF_KEYS
ORIENTED={"same","swapped","complement_same","complement_swapped"}
KEEP_FREQ={"same","complement_same"}

def diagnostic(fields,reference_hrs,source_direction):
    """Variant summary categories, no SNP ID or association estimates retained."""
    if source_direction not in ORIENTED:
        raise ValueError("cannot check source EAF without verified allele orientation")
    eaf,beta,se,p=(float(fields[i]) for i in (2,3,4,5))
    if (not all(map(math.isfinite,(eaf,beta,se,p)))
        or not 0<=eaf<=1 or se<=0 or not 0<=p<=1):
        raise ValueError("invalid allele frequency or association for EAF/GC check")
    eaf_aligned=eaf if source_direction in KEEP_FREQ else 1-eaf
    maf=min(eaf,1-eaf)
    maf_cat=("maf_lt_001" if maf<0.01 else
             "maf_001_to_005" if maf<0.05 else "maf_ge_005")
    if reference_hrs is None:
        hrs_cat="hrs_reference_not_available"
    else:
        hrs=float(reference_hrs)
        if not math.isfinite(hrs) or not 0<=hrs<=1:
            raise ValueError("invalid HRS source reference-allele frequency")
        delta=abs(eaf_aligned-hrs)
        hrs_cat=("hrs_absdiff_lt_002" if delta<0.02 else
                 "hrs_absdiff_002_to_005" if delta<0.05 else
                 "hrs_absdiff_005_to_010" if delta<0.1 else
                 "hrs_absdiff_010_to_020" if delta<0.2 else
                 "hrs_absdiff_ge_020")
    z=beta/se
    if not math.isfinite(z):
        raise ValueError("nonfinite original Z from BETA/SE")
    if p==0:
        wald_cat="reported_p_is_zero"
    else:
        from_z=minus_log10_two_sided_p_from_z(z)
        delta=abs(-math.log10(p)-from_z)
        if not math.isfinite(delta):
            raise ValueError("nonfinite P versus Wald divergence")
        wald_cat=("abs_log10p_delta_lt_001" if delta<0.01 else
                  "abs_log10p_delta_001_to_010" if delta<0.1 else
                  "abs_log10p_delta_010_to_1" if delta<1 else
                  "abs_log10p_delta_ge_1")
    try:
        odds=float(fields[6])
        if not math.isfinite(odds) or odds<=0:
            raise ValueError("invalid original OR")
        diff=abs(beta-math.log(odds))
        or_cat=("beta_logOR_absdelta_lt_1e10" if diff<1e-10 else
                "beta_logOR_absdelta_1e10_to_1e6" if diff<1e-6 else
                "beta_logOR_absdelta_ge_1e6")
    except ValueError:
        or_cat="source_OR_nonpositive_or_nonfinite"
    bound="raw_eaf_boundary_zero_or_one" if eaf in (0,1) else "raw_eaf_strict_between_zero_and_one"
    return (maf_cat,hrs_cat,wald_cat,or_cat,bound)

def summarize(records_by_coord):
    """Only called after duplicate-coordinate quarantine to avoid double counts."""
    by_chr={ch:Counter() for ch in range(1,23)}
    total=Counter()
    for (ch,pos),row in records_by_coord.items():
        if ch not in by_chr or len(row)!=5:
            raise ValueError("invalid original effect-diagnostic chromosome/shape")
        for item in row:
            if item not in DIAGNOSTIC_KEYS:
                raise ValueError("invalid original effect-diagnostic classification")
            total[item]+=1
            by_chr[ch][item]+=1
    return {
        "matched_original_source_SNPs_with_valid_EAF_and_Wald_inputs":len(records_by_coord),
        "categories":dict(total),
        "by_chr":[{"chr":ch,"matched_diagnostic_SNPs":sum(by_chr[ch].get(k,0) for k in MAF_KEYS),
                   "categories":dict(by_chr[ch])} for ch in range(1,23)],
        "HRS_frequency_is_reference_only":True,
        "HRS_ancestry_and_cohort_vs_GIGASTROKE_frequency_equivalence_not_assumed":True,
        "Wald_P_discrepancies_are_NOT_proof_of_genomic_control":True,
        "missing_per_variant_N_eff":True,
        "source_original_effect_sizes_uploaded":False,
    }

def _number(counts,key):
    value=counts.get(key,0)
    if type(value) is not int or value<0:
        raise ValueError("nonnegative integer expected: "+key)
    return value

def review_diagnostics(report):
    # Independently require publisher full MD5, Locke direction and 22 chr attrition.
    directions=direction_counts(report)
    if (report.get("LDSC_h2_or_rg_executed") is not False
        or report.get("scientific_status")!="HOLD"):
        raise ValueError("scientific LDSC gate unexpectedly changed")
    results=[]
    for r in report["outcomes"]:
        trait=r["phenotype"]
        if trait not in PINNED or r.get("raw_md5_verified") is not True:
            raise ValueError("unrecognized original GWAS or missing publisher full MD5")
        evidence=r.get("eaf_gc_source_descriptive",{})
        if evidence.get("HRS_frequency_is_reference_only") is not True or \
           evidence.get("Wald_P_discrepancies_are_NOT_proof_of_genomic_control") is not True or \
           evidence.get("missing_per_variant_N_eff") is not True or \
           evidence.get("source_original_effect_sizes_uploaded") is not False:
            raise ValueError("EAF / GC source assumptions lost")
        count=evidence.get("matched_original_source_SNPs_with_valid_EAF_and_Wald_inputs")
        if count!=r["counts"]["unique_harmonizable_BMI_positions"]:
            raise ValueError("valid frequency / Wald SNP set differs from full original source")
        groups=evidence.get("categories")
        if not isinstance(groups,dict) or set(groups)-set(DIAGNOSTIC_KEYS):
            raise ValueError("unrecognized diagnostic category")
        for key in DIAGNOSTIC_KEYS:_number(groups,key)
        for keys in (MAF_KEYS,HRS_KEYS,WALD_KEYS,LOGOR_KEYS,EAF_KEYS):
            if sum(_number(groups,k) for k in keys)!=count:
                raise ValueError("whole-genome EAF/GC category partition mismatch")
        chroms=evidence.get("by_chr")
        if not isinstance(chroms,list) or len(chroms)!=22 or [x.get("chr") for x in chroms]!=list(range(1,23)):
            raise ValueError("22-chromosome EAF/GC evidence missing")
        cumulative=Counter()
        for row in chroms:
            raw=row.get("categories",{})
            if not isinstance(raw,dict) or set(raw)-set(DIAGNOSTIC_KEYS):
                raise ValueError("unrecognized chromosome QC category")
            n=row.get("matched_diagnostic_SNPs")
            if type(n) is not int or n<0:
                raise ValueError("invalid chromosome observed SNP count")
            for keys in (MAF_KEYS,HRS_KEYS,WALD_KEYS,LOGOR_KEYS,EAF_KEYS):
                if sum(_number(raw,k) for k in keys)!=n:
                    raise ValueError("chromosome EAF/GC category partition mismatch")
            if n!=r["chr_attrition_pre_munge"][row["chr"]-1]["allele_aligned_numeric_valid_pre_munge"]:
                raise ValueError("chromosome EAF/GC count differs from original MD5 SNP ledger")
            cumulative.update(raw)
        if any(_number(groups,k)!=cumulative.get(k,0) for k in DIAGNOSTIC_KEYS):
            raise ValueError("chromosome and genome-wide EAF/GC counts disagree")
        results.append({"phenotype":trait,"accession":PINNED[trait][0],
                        "full_publisher_MD5_verified":True,
                        "matched_valid_SNPs":count,
                        "MAF_distribution":{k:_number(groups,k) for k in MAF_KEYS},
                        "HRS_reference_freq_absdiff":{k:_number(groups,k) for k in HRS_KEYS},
                        "P_vs_BETA_SE_Wald":{k:_number(groups,k) for k in WALD_KEYS},
                        "BETA_vs_logOR":{k:_number(groups,k) for k in LOGOR_KEYS},
                        "EAF_boundary_frequency":_number(groups,"raw_eaf_boundary_zero_or_one"),
                        "by_chr":[{"chr":x["chr"],
                                   "n":x["matched_diagnostic_SNPs"],
                                   "P_vs_Wald_discrepancy_ge_1":
                                      _number(x["categories"],"abs_log10p_delta_ge_1"),
                                   "HRS_absdiff_ge_020":
                                      _number(x["categories"],"hrs_absdiff_ge_020")}
                                  for x in chroms]})
    if {x["phenotype"] for x in results}!=set(PINNED):
        raise ValueError("missing four original EUR outcomes")
    return {"schema":"ckm-gigastroke-four-full-eaf-hrs-wald-descriptive-v1",
            "technical_status":"PASS_ORIGINAL_MD5_FOUR_EAF_AND_WALD_DESCRIPTIVE_ONLY",
            "outcomes":results,
            "original_EAF_not_directly_equivalent_to_GIANT_HRS_freq":True,
            "GC_method_or_correction_factor_INFERRED":False,
            "N_eff_approved":False,
            "LDSC_munged_regression_SNP_count_established":False,
            "LDSC_h2_intercept_rg_run":False,
            "scientific_status":"HOLD","ready_for_ldsc":False}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--evidence",required=True,type=Path)
    parser.add_argument("--out",required=True,type=Path)
    a=parser.parse_args()
    try:
        r=review_diagnostics(json.loads(a.evidence.read_text("utf-8")))
    except Exception as exc:
        r={"technical_status":"HOLD_EAF_GC_SOURCE_OR_DIAGNOSTIC_PROVENANCE",
           "error_type":type(exc).__name__,"error_detail":str(exc)[:240],
           "scientific_status":"HOLD","ready_for_ldsc":False,
           "N_eff_approved":False,"LDSC_h2_intercept_rg_run":False}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":r["technical_status"],
                      "results":[{"phenotype":x["phenotype"],
                                  "n":x["matched_valid_SNPs"],
                                  "p_vs_wald_discrepancy_ge1":
                                    x["P_vs_BETA_SE_Wald"]["abs_log10p_delta_ge_1"],
                                  "HRS_diff_ge02":
                                    x["HRS_reference_freq_absdiff"]["hrs_absdiff_ge_020"]}
                                 for x in r.get("outcomes",[])],
                      "scientific":r["scientific_status"],"error":r.get("error_detail")}))
    if r["technical_status"]!="PASS_ORIGINAL_MD5_FOUR_EAF_AND_WALD_DESCRIPTIVE_ONLY":
        raise SystemExit(1)

if __name__=="__main__":
    main()
