#!/usr/bin/env python3
"""Fail-closed Stage24 cross-modality descriptive correlation audit.

Accept aggregate-only results and enforce denominator consistency. This is
NOT a gene-level inferential, differential expression or recurrence analysis.
"""
import argparse
import json
import math
import pathlib

def verify(mapping, result):
    if mapping.get("audit_status")!="PASS_AGGREGATE_COUNTS_ONLY":
        raise ValueError("source-reviewed crosswalk did not pass Stage23")
    if result.get("schema")!="hcc-ng-stage24-descriptive-spearman-v1":
        raise ValueError("unexpected correlation schema")
    if result.get("source_RDS_sha256")!="461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a":
        raise ValueError("source RDS sha256 changed")
    if result.get("mapping_scope")!="source_RNA_x_current_UniProt_reviewed_strict_one_to_one":
        raise ValueError("annotation provenance unclear")
    pairs=mapping.get("strict_one_gene_one_protein_genes_n")
    if type(pairs)!=int or pairs<5000 or pairs>5500 or result.get("crosswalk_protein_gene_pair_n")!=pairs:
        raise ValueError("mapping denominator mismatch")
    if result.get("source_biopsy_n")!=50 or result.get("independent_patient_n")!=48:
        raise ValueError("biopsy/patient denominator mismatch")
    extra=result.get("patient_repeated_biopsy_groups_n")
    singleton=result.get("singleton_only_patient_n")
    if not type(extra)==int or not type(singleton)==int or not (1<=extra<=2 and singleton==48-extra):
        raise ValueError("duplicate-patient control mismatch")
    for label in ("patient_mean_48","singleton_only","biopsy_50_nonindependent_sensitivity"):
        stats=result.get(label)
        if not isinstance(stats,dict):raise ValueError("missing correlation stratum")
        finite=stats.get("finite_feature_n");undefined=stats.get("undefined_constant_feature_n")
        if type(finite)!=int or type(undefined)!=int or finite<5000 or finite+undefined!=pairs:
            raise ValueError("invalid finite correlation count")
        for k in ("median","q10","q90"):
            if type(stats.get(k)) not in (float,int) or not math.isfinite(stats[k]) or abs(stats[k])>1:
                raise ValueError("impossible Spearman rho")
        if not stats["q10"] <= stats["median"] <= stats["q90"]:
            raise ValueError("invalid correlation quantile order")
        for k in ("fraction_positive","fraction_abs_ge_0_5"):
            v=stats.get(k)
            if type(v) not in (float,int) or not math.isfinite(v) or not 0<=v<=1:
                raise ValueError("invalid rho proportion")
    delta=result.get("median_rho_patient_mean_minus_singleton_only")
    if type(delta) not in (float,int) or abs(delta-(result["patient_mean_48"]["median"]-result["singleton_only"]["median"]))>0.000005:
        raise ValueError("patient duplicate sensitivity summary mismatch")
    must_be_false=(
        "full_gene_universe_original_Ng_reproduced",
        "original_Ng_5481_gene_correlation_reproduced",
        "multiple_testing_correction_or_gene_level_pvalues_calculated",
        "direct_etiology_or_recurrence_tested",
        "scientific_significance_or_causality_proven",
        "feature_patient_identifiers_or_expression_exported",
    )
    for k in must_be_false:
        if result.get(k) is not False:
            raise ValueError("prohibited biological/privacy claim: "+k)
    return {
        "schema":"hcc-ng-stage24-postflight-review-v1",
        "source_stage23_count_audit":"PASS",
        "source_stage24_correlations_audited":"PASS_DESCRIPTIVE_ONLY",
        "mapped_strict_gene_protein_pairs_n":pairs,
        "biopsies":50,
        "unique_patients":48,
        "duplicate_patient_groups":extra,
        "singleton_only_patients":singleton,
        "patient_level_median_rho":result["patient_mean_48"]["median"],
        "singleton_only_median_rho":result["singleton_only"]["median"],
        "original_paper_gene_set_replication":"HOLD",
        "original_Ng_correlation_replication":"NOT_REPRODUCED",
        "V03":"PARTIAL",
        "V07":"PARTIAL_DESCRIPTIVE_NO_SIGNIFICANCE",
        "clinical_and_causal_conclusions":False,
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mapping",required=True)
    ap.add_argument("--correlations",required=True)
    ap.add_argument("--out",required=True)
    args=ap.parse_args()
    mapping=json.loads(pathlib.Path(args.mapping).read_text(encoding="utf-8"))
    corr=json.loads(pathlib.Path(args.correlations).read_text(encoding="utf-8"))
    result=verify(mapping,corr)
    out=pathlib.Path(args.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))

if __name__=="__main__":main()
