#!/usr/bin/env python3
"""Fail-closed descriptive patient resampling audit; not biological inference."""
import argparse
import json
import math
import pathlib

SOURCE_SHA="461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a"

def validate(stage24,stage25):
    if stage24.get("schema")!="hcc-ng-stage24-descriptive-spearman-v1":
        raise ValueError("Stage24 source result missing")
    if stage25.get("schema")!="hcc-ng-stage25-patient-influence-v1":
        raise ValueError("unknown Stage25 metrics schema")
    if stage25.get("source_RDS_sha256")!=SOURCE_SHA or stage24.get("source_RDS_sha256")!=SOURCE_SHA:
        raise ValueError("unpinned source RDS")
    if stage25.get("mapping_scope")!=stage24.get("mapping_scope"):
        raise ValueError("crosswalk scope changed")
    n=stage25.get("strict_gene_protein_pair_n")
    if type(n)!=int or n<5000 or n!=stage24.get("crosswalk_protein_gene_pair_n"):
        raise ValueError("wrong gene/protein pair denominator")
    if stage25.get("patient_n")!=48 or stage25.get("leave_one_patient_out_iterations")!=48:
        raise ValueError("invalid independent patient resampling")
    if stage25.get("patient_bootstrap_resamples")!=80 or stage25.get("patient_bootstrap_seed")!=20261031:
        raise ValueError("bootstrap not reproducible")
    b=stage24["patient_mean_48"]["median"]
    p=stage25.get("primary_median_spearman_rho")
    if type(p) not in (float,int) or abs(p-b)>0.000002:
        raise ValueError("Stage24 vs Stage25 baseline rho drift")
    nums=["leave_one_out_median_min","leave_one_out_median_max",
          "leave_one_out_median_q05","leave_one_out_median_q95",
          "leave_one_out_max_absolute_shift",
          "patient_bootstrap_median_rho_q025","patient_bootstrap_median_rho_q975"]
    for k in nums:
        v=stage25.get(k)
        if type(v) not in (float,int) or not math.isfinite(v) or abs(v)>2:
            raise ValueError("unphysical statistic: "+k)
    m=stage25
    if not (-1<=m["leave_one_out_median_min"]<=m["leave_one_out_median_q05"]
           <=m["leave_one_out_median_q95"]<=m["leave_one_out_median_max"]<=1):
        raise ValueError("invalid leave-one-out ordered quantiles")
    if not (-1<=m["patient_bootstrap_median_rho_q025"]
           <=m["patient_bootstrap_median_rho_q975"]<=1):
        raise ValueError("invalid patient bootstrap quantiles")
    shift=m["leave_one_out_max_absolute_shift"]
    if shift<0 or shift+0.000002<max(
         abs(m["leave_one_out_median_min"]-p),
         abs(m["leave_one_out_median_max"]-p)):
        raise ValueError("LOO shift smaller than observed extrema")
    if type(m.get("leave_one_out_max_shift_gt_0_05_review_flag")) is not bool or (
        m["leave_one_out_max_shift_gt_0_05_review_flag"]!=(shift>0.05)):
        raise ValueError("LOO review flag inconsistent")
    required_false=[
       "technical_preprocessing_robustness_proven",
       "original_ng_correlation_reproduced",
       "gene_level_FDR_or_per_gene_CI_calculated",
       "clinical_or_causal_validity_proven",
       "patient_gene_rho_or_expression_values_exported",
    ]
    if any(m.get(k) is not False for k in required_false):
        raise ValueError("forbidden scientific or privacy claim")
    return {
      "schema":"hcc-ng-stage25-patient-influence-review-v1",
      "source_validation":"PASS_SOURCE_AND_AGGREGATE_COUNTS",
      "patient_bootstrap_validation":"PASS_RESAMPLING_DESCRIPTIVE_ONLY",
      "clinical_or_biological_validity":"NOT_TESTED",
      "bootstrap_95_percentile_interval_not_automatically_95_coverage":True,
      "reanalysis_original_ng_source_gene_set":"HOLD",
      "Suter_s5_preprocessing_parameters_verified":"HOLD",
      "V03":"PARTIAL","V07":"PARTIAL",
      "review_flag_large_single_patient_shift":m["leave_one_out_max_shift_gt_0_05_review_flag"],
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--stage24",required=True)
    ap.add_argument("--stage25",required=True)
    ap.add_argument("--out",required=True)
    a=ap.parse_args()
    x=json.loads(pathlib.Path(a.stage24).read_text(encoding="utf-8"))
    y=json.loads(pathlib.Path(a.stage25).read_text(encoding="utf-8"))
    report=validate(x,y)
    out=pathlib.Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))

if __name__=="__main__":main()
