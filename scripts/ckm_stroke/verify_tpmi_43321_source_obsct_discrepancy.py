#!/usr/bin/env python3
"""Independent aggregate cross-check of pinned all-SNP TPMI GWAS CI artifact.

Strictly consumes prior original-source aggregate JSON, not raw GWAS effects.
Separates publisher PheWeb phenotype N from PLINK2 regression OBS_CT;
NEITHER is permitted to stand in for variant-level binary N_eff.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

PINNED_SHA="ce14daef23f9fadf921030d9ed2429383d927f6da391c620d9abf870045b2e62"
SOURCE_RUN=44872892525
EXPECTED_ROWS=1_739_048
PHEWEB_CASES=9_249
PHEWEB_CONTROLS=304_660
OBS_CT=204_099
N_CHR=22
EXPECTED_HEADER_COLUMNS=17

def validate(report):
    if (report.get("schema")!="TPMI43321-original-gwas-complete-row-QC-v1"
        or report.get("technical_status")!="PASS_FULL_SOURCE_SHA256_AND_ALL_ROWS_PARSED_WITH_QC_FLAGS"
        or report.get("full_original_compressed_sha256")!=PINNED_SHA
        or report.get("pinned_prior_original_sha256_match") is not True
        or report.get("original_compressed_bytes")!=58_302_393
        or report.get("rows_total_in_original")!=EXPECTED_ROWS
        or len(report.get("pinned_original_columns",[]))!=EXPECTED_HEADER_COLUMNS
        or report.get("scientific_status")!="HOLD"
        or report.get("ready_for_ldsc") is not False
        or report.get("MR_coloc_LDSC_executed") is not False):
        raise ValueError("original source integrity, complete genome or scientific HOLD drift")
    bychr=report.get("by_chr_all_original_rows",{})
    if (len(bychr)!=N_CHR
        or set(bychr)!={str(x) for x in range(1,23)}
        or sum(bychr.values())!=EXPECTED_ROWS):
        raise ValueError("22-chromosome original source count ledger disagrees")
    c=report.get("row_level_QC_counts",{})
    for label in ("records","fully_numeric_ADD_no_error","P_valid",
                  "OR_positive","LOGOR_SE_positive","OBS_CT_valid",
                  "freq_numeric","P_vs_Z_comparable"):
        if c.get(label)!=EXPECTED_ROWS:
            raise ValueError("whole-original association data complete-case QC drift: "+label)
    if (c.get("P_le_5e_8")!=1
        or report.get("duplicate_variant_test_128bit_hashes")!=0
        or report.get("original_additive_test_categories")!={"ADD":EXPECTED_ROWS}
        or report.get("original_ERRORCODE_categories")!={".":EXPECTED_ROWS}
        or report.get("source_rows_out_of_coordinate_order")!=0):
        raise ValueError("original source mutation or sort/TEST condition drift")
    freq=report.get("original_MAF_bins_from_A1_FREQ",{})
    if (sum(freq.values())!=EXPECTED_ROWS
        or any(not isinstance(v,int) or v<0 for v in freq.values())):
        raise ValueError("MAF bins do not account for each genomic original row")
    extrema=report.get("value_min_max",{})
    if extrema.get("min_OBS_CT")!=OBS_CT or extrema.get("max_OBS_CT")!=OBS_CT:
        raise ValueError("PLINK2 original OBS_CT is not the previously verified 204099 constant")
    if report.get("original_full_file_deleted_after_scan") is not True:
        raise ValueError("full raw GWAS was not removed from ephemeral audit")
    total=PHEWEB_CASES+PHEWEB_CONTROLS
    delta=total-OBS_CT
    if total!=313_909 or delta!=109_810:
        raise ValueError("PheWeb phenotype and PLINK2 N mismatch math drift")
    return {
      "schema":"ckm-TPMI-43321-source-QC-sample-denominator-reconciliation-v1",
      "technical_status":"PASS_PINNED_SOURCE_FULL_QC_RECONCILED_SAMPLE_DENOMINATOR_MISMATCH",
      "original_full_source_run":SOURCE_RUN,
      "original_full_source_SHA256":PINNED_SHA,
      "original_all_SNPs_checked":EXPECTED_ROWS,
      "original_individual_chr_reference_ledger_passed":True,
      "source_logistic_tests":"ADD",
      "source_errorcode_all_rows":".",
      "original_single_variant_P_le_5e8_count":1,
      "PheWeb_phenotype_reported_ncases":PHEWEB_CASES,
      "PheWeb_phenotype_reported_ncontrols":PHEWEB_CONTROLS,
      "PheWeb_phenotype_reported_total_n":total,
      "PLINK2_per_variant_regression_OBS_CT_all_SNPs":OBS_CT,
      "PLINK2_OBS_CT_min_equals_max":True,
      "observed_total_difference":delta,
      "observed_total_difference_percent_of_page_N":round(100*delta/total,5),
      "PLINK2_sample_could_be_an_unrelated_subset_not_source_confirmed":True,
      "TPMI_2025_original_cohort_205779_independent_reported_in_Nature":True,
      "infer_exact_204099_from_205779_NOT_ALLOWED":True,
      "phenotype_N_vs_regression_N_reconciled_by_individual_records":False,
      "SNP_case_specific_N_available":False,
      "SNP_control_specific_N_available":False,
      "binary_case_control_Neff_approved":False,
      "MAF_bins":freq,
      "effect_allele":"A1_IN_PLINK2_NOT_AUTOMATICALLY_ALT",
      "OR_NEEDS_POSITIVE_LOG_TRANSFORM_BEFORE_MR":True,
      "future_mr_coloc_needs_full_case_definition_QC_and_harmonization":True,
      "study_overlap_with_GIGASTROKE_excluded_by_individual_ID":False,
      "actual_LDSC_h2_intercept_rg_MR_coloc_not_executed":True,
      "technical_full_gwas_input_check":"PASS",
      "scientific_primary_GWAS_Neff_and_independent_EAS_replication":"HOLD",
      "scientific_status":"HOLD","ready_for_ldsc":False
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--full-source-artifact",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    try:
        r=validate(json.loads(a.full_source_artifact.read_text("utf-8")))
    except Exception as e:
        r={"technical_status":"HOLD_PINNED_TPMI_SOURCE_AGGREGATE_AUDIT",
           "error_type":type(e).__name__,"error_detail":str(e)[:250],
           "scientific_status":"HOLD","ready_for_ldsc":False}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n","utf-8")
    print(json.dumps({
      "status":r["technical_status"],
      "source_SNPs":r.get("original_all_SNPs_checked"),
      "OBS_CT":r.get("PLINK2_per_variant_regression_OBS_CT_all_SNPs"),
      "PheWeb_total":r.get("PheWeb_phenotype_reported_total_n"),
      "mismatch_percent":r.get("observed_total_difference_percent_of_page_N"),
      "source_Neff_scientific_HOLD":r["scientific_status"],
      "error":r.get("error_detail")
    }))
    if r["technical_status"]!="PASS_PINNED_SOURCE_FULL_QC_RECONCILED_SAMPLE_DENOMINATOR_MISMATCH":
        raise SystemExit(1)
if __name__=="__main__":
    main()
