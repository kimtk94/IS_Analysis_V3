#!/usr/bin/env python3
"""Read-only ALDH2 BBJ ischemic-stroke conditional LD sensitivity science gate.

Input Phase8F is BBJ ischemic stroke (AIS), not Koyanagi alcohol GWAS.
Contrasts single-lead conditioned summary statistics using 1000G JPT104
vs EAS504 proxy LD, alongside separate Koyanagi alcohol SuSiE reliability.
NO claims that PIP=1 implies causal independence.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path

CUTOFF=5e-8
LOW_LD_R2=0.2
EXPECTED_ROWS=3033
EXPECTED_MODELS=4

def load(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))

def analyze(rows,ld,cutoff=CUTOFF,low_ld=LOW_LD_R2):
    fields={f"{ld}_p_cond",f"{ld}_r2_to_rs671",f"{ld}_z_cond"}
    if not rows or not fields.issubset(rows[0]):
        raise ValueError(f"Unexpected conditional schema: {ld}")
    valid=[]
    for r in rows:
        try:
            p=float(r[f"{ld}_p_cond"])
            r2=float(r[f"{ld}_r2_to_rs671"])
            z=float(r[f"{ld}_z_cond"])
        except (ValueError,TypeError) as e:
            raise ValueError("Invalid conditional statistic") from e
        if not (math.isfinite(p) and 0<=p<=1 and
                math.isfinite(r2) and 0<=r2<=1+1e-7 and math.isfinite(z)):
            raise ValueError("Out of range conditional statistic")
        valid.append((p,r2,r,z))
    if not valid: raise ValueError("Empty LD reference")
    hits=[x for x in valid if x[0] < cutoff]
    independent=[x for x in hits if x[1] <= low_ld]
    top=sorted(valid,key=lambda x:x[0])[:5]
    return {
        "variants":len(valid),"p_cutoff":cutoff,"low_ld_r2_max":low_ld,
        "gws_residual":len(hits),
        "gws_low_ld_residual":len(independent),
        "gws_high_ld_residual":len(hits)-len(independent),
        "top_five":[{
            "variant_id":x[2]["variant_id"],"rsid":x[2]["rsid"],
            "conditional_p":x[0],"r2_to_rs671":x[1],
            "conditional_z":x[3],
            "gws_and_low_ld":x[0]<cutoff and x[1]<=low_ld
        } for x in top],
        "gws_variant_ids":[x[2]["variant_id"] for x in hits],
        "gws_low_ld_variant_ids":[x[2]["variant_id"] for x in independent],
    }

def science_gate(rows, condition_gate, finite_gate, stable):
    if len(rows)!=EXPECTED_ROWS or len({r["variant_id"] for r in rows})!=len(rows):
        raise ValueError("Expected 3033 distinct BBJ IS SNPs")
    j=analyze(rows,"jpt")
    e=analyze(rows,"eas")
    if int(condition_gate["ALDH2"]["source_input_variants"])!=574:
        raise ValueError("Unexpected alcohol source variants")
    if finite_gate["total_models"]!=EXPECTED_MODELS:
        raise ValueError("Unexpected SuSiE model count")
    if not finite_gate["reliability_warning_all_4_models"] or not finite_gate["sensitivity_warning_all_4_models"]:
        raise ValueError("SuSiE reliability flags changed; re-audit required")
    alc=stable["results"]["ALDH2"]
    if alc["status"]!="BLOCKED" or not all(alc["R_reliability_flag"]):
        raise ValueError("Koyanagi ALDH2 fine-mapping science gate unexpectedly changed")
    residual_flag="HIGH_LD_PROXY_ONLY" if e["gws_residual"]>0 and e["gws_low_ld_residual"]==0 else "REVIEW"
    status=("NO_INDEPENDENT_SECONDARY_SIGNAL_SUPPORTED_WITH_EXISTING_PROXY_LD"
        if j["gws_low_ld_residual"]==0 and e["gws_low_ld_residual"]==0
        else "LOW_LD_SECONDARY_SIGNAL_REVIEW_REQUIRED")
    return {
        "status":status,
        "trait_BBJ_conditioned":"ISCHEMIC_STROKE_BBJ",
        "trait_Koyanagi_finemapping":"ALCOHOL_CONSUMPTION",
        "ld_reference_JPT":"1000G_JPT104",
        "ld_reference_EAS":"1000G_EAS504",
        "sensitivity_JPT":j,
        "sensitivity_EAS":e,
        "EAS_residual_interpretation":residual_flag,
        "single_lead_conditioning":"EXTERNAL_LD_APPROXIMATION_NOT_FULL_JOINT_MULTI_SIGNAL",
        "locus_alcohol_finemap_status":alc["status"],
        "alcohol_finemap_cs_counts":alc["credible_set_counts"],
        "alcohol_finemap_R_reliability_flag":alc["R_reliability_flag"],
        "alcohol_finemap_R_sensitivity_flag":alc["R_sensitivity_flag"],
        "alcohol_gwas_ld_mismatch_s":float(condition_gate["ALDH2"]["LD_mismatch_s"]),
        "instrument_independence_demonstrated":False,
        "second_alcohol_causal_signal_demonstrated":False,
        "secondary_BBJ_IS_signal_demonstrated":False,
        "ALDH2_BLOCKED_lifted":False,
        "external_ld_not_cohort_matched":True,
        "Koyanagi_BBJ_sample_overlap_unresolved":True,
        "large_Japanese_cohort_matched_LD_required":True,
        "minimum_next_validation":[
            "cohort_matched_adequate_Japanese_LD",
            "multi_signal_conditional_analysis_or_joint_finemapping",
            "independent_non_BBJ_IS_outcome_for_replication",
            "explicit_sample_overlap_handling_before_MR",
            "re-evaluate_alcohol_SuSiE_model_reliability_and_sensitivity"
        ],
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--conditional-tsv",type=Path,required=True)
    p.add_argument("--alcohol-condition-gate",type=Path,required=True)
    p.add_argument("--finite-gate",type=Path,required=True)
    p.add_argument("--finite-stability",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    rows=load(a.conditional_tsv)
    report=science_gate(rows,
        json.loads(a.alcohol_condition_gate.read_text()),
        json.loads(a.finite_gate.read_text()),
        json.loads(a.finite_stability.read_text()))
    if a.out_dir.exists() and any(a.out_dir.iterdir()):
        raise FileExistsError("Refuse to overwrite existing audit")
    a.out_dir.mkdir(parents=True,exist_ok=True)
    report["utc"]=datetime.now(timezone.utc).isoformat()
    report["sources"]={
        "BBJ_IS_single_lead_conditioned":str(a.conditional_tsv),
        "Koyanagi_alcohol_gate":str(a.alcohol_condition_gate),
        "finite_R_reliability_gate":str(a.finite_gate),
        "finite_R_stability":str(a.finite_stability)
    }
    report["scope"]="DIAGNOSTIC_ONLY_NO_CANONICAL_MUTATION"
    output=a.out_dir/"ALDH2_BBJ_IS_AND_ALCOHOL_CONDITIONAL_SCIENCE_GATE.json"
    output.write_text(json.dumps(report,indent=2)+"\n")
    print("ALDH2_CONDITIONAL_LD_GATE_COMPLETE",report["status"])
    for k in ("sensitivity_JPT","sensitivity_EAS"):
        x=report[k]
        print(k, "gws_residual",x["gws_residual"],
              "low_ld_gws",x["gws_low_ld_residual"],
              "top_p",x["top_five"][0]["conditional_p"],
              "top_r2",x["top_five"][0]["r2_to_rs671"])
    print("ALDH2_FINE_MAPPING",report["locus_alcohol_finemap_status"])

if __name__=="__main__":
    main()
