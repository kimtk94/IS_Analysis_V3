#!/usr/bin/env python3
"""Validate original GIGASTROKE effective-sample-size evidence, never estimate h2."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
from scripts.ckm_stroke.audit_gigastroke_supplement_cohort_neff import (
    EXPECTED_CONTROLS, EXPECTED_SHA256, EXPECTED_SUBTYPES
)
from scripts.ckm_stroke.audit_gigastroke_bmi_allele_harmonization import COLUMNS

MISSING_CASES = {"AIS": 103, "LAA": 0, "CE": 28, "SVD": 39}

def numeric(value, label):
    if isinstance(value, bool) or not isinstance(value, (int,float)):
        raise ValueError("invalid numeric " + label)
    x = float(value)
    if not math.isfinite(x):
        raise ValueError("nonfinite numeric " + label)
    return x

def mapped(records, label):
    if not isinstance(records, list) or len(records) != 4 or any(not isinstance(v,dict) for v in records):
        raise ValueError("four original " + label + " required")
    names = [v.get("phenotype") for v in records]
    if len(set(names)) != 4 or set(names) != set(EXPECTED_SUBTYPES):
        raise ValueError("duplicate or missing original " + label)
    return {v["phenotype"]: v for v in records}

def evaluate(cohort, effect):
    """The only successful scientific outcome is HOLD, with no numerical N approval."""
    if (cohort.get("schema") != "ckm-gigastroke-corrected-primary-cohort-ldsc-gate-v1"
        or cohort.get("publisher_source_sha256") != EXPECTED_SHA256
        or cohort.get("technical_status") != "PASS_PRIMARY_CORRECTED_COHORT_AND_TABLE3_ONLY"
        or cohort.get("publisher_original_table1_eur_cohort_count") != 44
        or cohort.get("original_full_sumstats_downloaded") is not False
        or cohort.get("ready_for_ldsc") is not False
        or cohort.get("LDSC_h2_rg_executed") is not False
        or cohort.get("scientific_status") != "HOLD"):
        raise ValueError("unverified corrected first-party cohort evidence")
    if (effect.get("schema") != "ckm-gigastroke-eur-bounded-original-logor-scale-v1"
        or effect.get("technical_status") != "PASS_FOUR_BOUNDED_LOGOR_EFFECT_SCALE_ONLY"
        or effect.get("full_source_verified_in_this_run") is not False
        or effect.get("LDSC_h2_intercept_rg_executed") is not False
        or effect.get("ready_for_ldsc") is not False
        or effect.get("scientific_status") != "HOLD"):
        raise ValueError("unverified original EBI effect prefix")
    ukb = cohort.get("uk_biobank",{})
    if (ukb.get("AIS_case_count") != 2211
        or ukb.get("UKB_reported_control_count") != 369419
        or ukb.get("overlap_in_BMI_actual_count_known") is not False):
        raise ValueError("UKB overlap provenance unsupported")
    outcomes = mapped(cohort.get("outcomes"), "outcomes")
    probes = mapped(effect.get("prefixes"), "effect prefixes")
    records = []
    for phenotype, spec in EXPECTED_SUBTYPES.items():
        o = outcomes[phenotype]
        p = probes[phenotype]
        cases, controls = spec["cases"], EXPECTED_CONTROLS
        if (o.get("accession") != spec["accession"]
            or p.get("accession") != spec["accession"]
            or o.get("EUR_case_count_from_original_table1") != cases
            or o.get("EUR_control_count_from_original_table1") != controls
            or o.get("EUR_total_N_from_original_table1") != spec["total"]
            or o.get("illustrative_neff_not_per_snp_or_approved") is not True
            or o.get("this_pipeline_ldsc_executed") is not False):
            raise ValueError("phenotype or original sample metadata changed " + phenotype)
        if (p.get("source_per_variant_N_available") is not False
            or p.get("source_rsid_available") is not False
            or p.get("source_header_columns") != list(COLUMNS)
            or p.get("outcome_beta_scale") != "LOG_ODDS_ratio_approx_beta"
            or p.get("sample_prefix_compressed_bytes") != 131072
            or p.get("bounded_sample_not_genomewide_proof") is not True
            or p.get("observed",{}).get("parsed_effect_rows",0) < 500
            or numeric(p.get("maximum_abs_beta_minus_lnOR"), phenotype+" effect") > 1e-10):
            raise ValueError("per-variant N missing or effect-source mismatch " + phenotype)
        ci = o.get("illustrative_cohortwise_neff",{})
        missing = MISSING_CASES[phenotype]
        rows = ci.get("missing_cohort_control_rows")
        if (ci.get("reported_case_total_from_rows") != cases
            or ci.get("cases_missing_cohort_control_count") != missing
            or ci.get("cases_with_reported_cohort_controls") != cases-missing
            or ci.get("not_per_variant_N_and_not_approved_for_LDSC") is not True
            or not isinstance(rows,list)
            or any(not isinstance(r,dict) or not isinstance(r.get("cases"),int)
                   or isinstance(r.get("cases"),bool) or r["cases"]<=0 for r in rows)
            or sum(r["cases"] for r in rows)!=missing):
            raise ValueError("missing cohort controls misrepresented " + phenotype)
        pooled = numeric(o.get("illustrative_harmonic_neff"), phenotype+" pooled")
        per_cohort = numeric(ci.get("cohortwise_sum_harmonic_neff_illustrative"), phenotype+" cohort")
        if (abs(pooled-round(4*cases*controls/(cases+controls),3))>0.001
            or not 0 < per_cohort < pooled):
            raise ValueError("illustrative N_eff scenario invalid " + phenotype)
        records.append({
            "phenotype": phenotype, "accession": spec["accession"],
            "source_EUR_cases":cases, "source_EUR_controls":controls,
            "source_total_N_is_NOT_variant_N":spec["total"],
            "missing_controls_cohort_case_rows":missing,
            "original_logOR_sample_checked":True,
            "source_per_variant_N_available":False,
            "illustrative_not_approved_for_LDSC":{
                "total_first_case_control_harmonic":pooled,
                "sum_reported_cohort_harmonics":per_cohort,
                "absolute_difference":round(pooled-per_cohort,3),
                "relative_difference_pct":round(100*(pooled-per_cohort)/pooled,3),
            },
            "effective_N_model":"UNIDENTIFIED",
            "scientific_gate":"HOLD_PER_VARIANT_META_WEIGHTED_N_UNVERIFIED",
            "allow_pooled_N_as_LDSC_input":False,
            "allow_EAC_as_LDSC_N":False,
        })
    return {
        "schema":"ckm-gigastroke-meta-neff-method-policy-v1",
        "checked_at_utc":datetime.now(timezone.utc).isoformat(),
        "technical_status":"PASS_ORIGINAL_META_NEFF_EVIDENCE_BLOCKS_UNJUSTIFIED_N",
        "publisher_corrected_supplement_sha256":EXPECTED_SHA256,
        "original_study_method_references":{
            "paper":"https://doi.org/10.1038/s41586-022-05165-3",
            "meta_analysis":"fixed_effect_inverse_variance",
            "reported_cov_LDSC":"on_meta_analysis_before_genomic_control_correction",
            "effective_allele_count_QC_is_NOT_effective_N":True,
            "paper_methods_citation_documented_not_live_text_verified_by_CI":True,
        },
        "outcomes":records,
        "next_requirements":[
            "Recover meta-analysis variant-specific participating cohorts and weighting or justifiable N model",
            "Review original Cov-LDSC without genomic control versus public harmonized file",
            "Verify final 22-chromosome reference and regression SNP attrition",
            "Review cohort overlap and report analysis sensitivity/assumptions before any h2",
        ],
        "scientific_status":"HOLD", "ready_for_ldsc":False,
        "approved_N_eff_available":False,
        "LDSC_h2_intercept_rg_run":False,
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--cohort",required=True,type=Path)
    parser.add_argument("--prefix",required=True,type=Path)
    parser.add_argument("--out",required=True,type=Path)
    args=parser.parse_args()
    try:
        r=evaluate(json.loads(args.cohort.read_text("utf-8")),
                   json.loads(args.prefix.read_text("utf-8")))
    except Exception as e:
        r={"technical_status":"HOLD_INVALID_NEFF_EVIDENCE",
           "scientific_status":"HOLD", "ready_for_ldsc":False,
           "approved_N_eff_available":False, "LDSC_h2_intercept_rg_run":False,
           "error_type":type(e).__name__,"error_detail":str(e)[:240]}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":r["technical_status"],
                      "delta_neff_pct":{x["phenotype"]:x["illustrative_not_approved_for_LDSC"]["relative_difference_pct"] for x in r.get("outcomes",[])},
                      "scientific":r["scientific_status"],"error":r.get("error_detail")}))
    if r["technical_status"]!="PASS_ORIGINAL_META_NEFF_EVIDENCE_BLOCKS_UNJUSTIFIED_N":
        raise SystemExit(1)

if __name__=="__main__":
    main()
