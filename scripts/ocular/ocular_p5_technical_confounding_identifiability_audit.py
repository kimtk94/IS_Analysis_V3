#!/usr/bin/env python3
"""GSE313736: exact design identifiability and pooled-cell selection decomposition.

Inputs: original public aggregate-only P5e selection JSON, P5b 9-pool LOO
summary and P5d pooled-library model summary. This does NOT fit a new causal
drug model, assume independent donor-libraries, or access genotypes/UMIs.
"""
import argparse
import hashlib
import json
import math
import tempfile
from pathlib import Path
from datetime import datetime, timezone
import numpy as np

POOLS = tuple("ABCDEFGHI")


def full_rank_audit():
    # 18 distinct original 10x libraries: 9 paired source pools x 2 arms.
    treatment = np.array([arm for pool in POOLS for arm in (0., 1.)])
    source_pool = np.repeat(np.eye(9), 2, axis=0)
    library = np.eye(18)
    pool_treat = np.column_stack((source_pool, treatment))
    library_treat = np.column_stack((library, treatment))
    r_pool = int(np.linalg.matrix_rank(pool_treat))
    r_library = int(np.linalg.matrix_rank(library_treat))
    # The treatment column is exactly the sum of all 9 DEX-library dummies.
    reconstruction = library[:, 1::2].sum(axis=1)
    deviation = float(np.max(np.abs(reconstruction - treatment)))
    if (r_pool, r_library, deviation) != (10, 18, 0.0):
        raise ValueError("DESIGN_ALIAS_RANK_OR_IDENTITY_DRIFT")
    if float(np.linalg.norm(treatment-library@np.linalg.lstsq(library,treatment,rcond=None)[0]))>1e-12:
        raise ValueError("DEX_NOT_IN_LIBRARY_SPAN")
    return {
        "n_pool_pairs": 9, "n_original_10x_libraries":18,
        "model_pool_plus_treatment_design_columns": 10,
        "model_pool_plus_treatment_rank":r_pool,
        "model_library_fixed_effect_plus_treatment_design_columns":19,
        "model_library_fixed_effect_plus_treatment_rank":r_library,
        "new_rank_from_treatment_after_unrestricted_library_effects":r_library-18,
        "treatment_equals_sum_of_9_DEX_library_indicator_columns":True,
        "unrestricted_library_fixed_effect_model_drug_effect_identifiable":False,
        "note": ("Donor pairing and 9-pool sensitivity models can estimate a DEX-arm "
                 "contrast only by assuming no unrestricted condition-linked "
                 "library-preparation nuisance; they cannot identify drug versus "
                 "arbitrary original-library effects from this source alone."),
    }


def selection_pool_decomposition(sel):
    if sel.get("status") != "PASS_CELL_STATE_SELECTION_DIAGNOSTIC":
        raise ValueError("SOURCE_P5E_SELECTION_QC_NOT_PASS")
    if sel.get("biological_donors_with_both_arms") != 59 or sel.get("cells") != 169398:
        raise ValueError("ORIGINAL_P5E_SOURCE_DENOMINATOR_DRIFT")
    source_arms=sel["arms"]
    pools=sel["per_pool_aggregate"]
    if len(pools)!=9 or {x["pool"] for x in pools}!=set(POOLS):
        raise ValueError("MISSING_OR_DUPLICATE_SOURCE_POOLS")
    rows=[]
    for x in sorted(pools,key=lambda z:z["pool"]):
        src0=x["no-dex"];src1=x["dex"]
        n0=int(src0["total_sng"]);n1=int(src1["total_sng"])
        if min(n0,n1)<=0:
            raise ValueError("EMPTY_POOL_ARM")
        p0=float(src0["beamA_frac_ge_0_75"]);p1=float(src1["beamA_frac_ge_0_75"])
        if not all(math.isfinite(y) and 0<=y<=1 for y in (p0,p1)):
            raise ValueError("INVALID_POOL_FRACTION")
        rows.append({"pool":x["pool"],"Vehicle_SNG":n0,"DEX_SNG":n1,
                     "Vehicle_BeamA_fraction":p0,"DEX_BeamA_fraction":p1,
                     "within_pool_delta_percentage_points":100*(p1-p0)})
    n0=sum(x["Vehicle_SNG"] for x in rows)
    n1=sum(x["DEX_SNG"] for x in rows)
    observed0=int(source_arms["no-dex"]["total_sng"])
    observed1=int(source_arms["dex"]["total_sng"])
    if n0!=observed0 or n1!=observed1:
        raise ValueError(f"POOL_COVERAGE_INCOMPLETE source={observed0}/{observed1} pool={n0}/{n1}")
    out0=source_arms["no-dex"]["beamA_fraction_ge_0.75"]
    out1=source_arms["dex"]["beamA_fraction_ge_0.75"]
    total0=sum(x["Vehicle_SNG"]/n0*x["Vehicle_BeamA_fraction"] for x in rows)
    total1=sum(x["DEX_SNG"]/n1*x["DEX_BeamA_fraction"] for x in rows)
    # P5e only emits rounded 6-decimal fractions; don't imply exact recovered counts.
    if max(abs(total0-float(out0)),abs(total1-float(out1)))>2e-6:
        raise ValueError("POOL_AND_GLOBAL_FRACTIONS_DISAGREE")
    within=mix=0.
    for x in rows:
        w0=x["Vehicle_SNG"]/n0;w1=x["DEX_SNG"]/n1
        p0=x["Vehicle_BeamA_fraction"];p1=x["DEX_BeamA_fraction"]
        x["Vehicle_weight"]=w0;x["DEX_weight"]=w1
        x["symmetric_within_contribution_pp"]=100*(w0+w1)/2*(p1-p0)
        x["symmetric_mix_contribution_pp"]=100*(p0+p1)/2*(w1-w0)
        within+=x["symmetric_within_contribution_pp"]
        mix+=x["symmetric_mix_contribution_pp"]
    raw_delta=100*(total1-total0)
    if abs(raw_delta-within-mix)>1e-8:
        raise ValueError("DECOMPOSITION_ARITHMETIC_FAILURE")
    if abs(raw_delta-2.4064)>0.005:
        raise ValueError("PUBLISHED_P5E_GLOBAL_DELTA_CHANGED")
    paired=sel["donor_paired_differences"]["BeamA_0.75"]
    paired_sign=paired["paired_sign_test"]
    if (paired_sign["positive"],paired_sign["negative"])!=(29,30):
        raise ValueError("ORIGINAL_DONOR_SIGN_COUNTS_DRIFT")
    median=100*paired["paired_fraction_delta_dex_minus_vehicle"]["median"]
    if abs(median+0.421)>0.005:
        raise ValueError("ORIGINAL_DONOR_MEDIAN_CHANGED")
    return {
        "status":"PASS_FULL_POOL_COVERAGE_ROUNDED_FRACTIONS",
        "total_sng_by_arm":{"Vehicle":n0,"DEX":n1},
        "weighted_pool_fraction_by_arm":{"Vehicle":total0,"DEX":total1},
        "pooled_BeamA_Dex_minus_Vehicle_percentage_points":raw_delta,
        "symmetric_within_pool_fraction_contribution_pp":within,
        "symmetric_pool_arm_weight_mix_contribution_pp":mix,
        "equal_pool_mean_delta_pp":sum(x["within_pool_delta_percentage_points"] for x in rows)/9,
        "positive_pool_deltas":sum(x["within_pool_delta_percentage_points"]>0 for x in rows),
        "negative_pool_deltas":sum(x["within_pool_delta_percentage_points"]<0 for x in rows),
        "paired_donor_n":59,
        "paired_donor_BeamA_fraction_median_delta_pp":median,
        "paired_donor_sign_positive":29,"paired_donor_sign_negative":30,
        "paired_sign_test_two_sided_p":paired_sign["two_sided_p"],
        "per_pool":rows,
        "decomposition_method":("Symmetric Kitagawa: sum avg(arm weights)*within-pool change "
                                "+ sum avg(pool proportions)*arm-weight change."),
        "caution":("Pool-level aggregate selection rates and SNG capture weighting, "
                   "not paired donor-level biological treatment effects; "
                   "fractional source rates rounded to 6 decimals; "
                   "no unique causal attribution."),
    }


def inspect_models(p5b,p5d):
    if p5b.get("status")!="PASS_NINE_LPO_EDGE_R_SENSITIVITY_NUMERIC":
        raise ValueError("P5B_NOT_QC_PASS")
    if p5d.get("status")!="PASS_9_POOL_PAIRED_EDGE_R_QC":
        raise ValueError("P5D_NOT_QC_PASS")
    md=p5b["mds"]
    if md.get("n_samples")!=114 or md.get("n_paired_donors")!=57:
        raise ValueError("MDS_SAMPLE_COUNT_CHANGED")
    r2=float(md["first_axis_treatment_r2"])
    if abs(r2-0.95557721)>0.00001:
        raise ValueError("MDS_AXIS_R2_SOURCE_DRIFT")
    if p5d["pooled_libraries"]!=18 or p5d["residual_df"]!=8 or p5d["common_genes"]!=16049:
        raise ValueError("NINE_POOL_UNIT_SOURCE_CHANGED")
    if p5d["FDR_genes_shared"]!=11079 or p5d["primary_FDR_genes_common"]!=13716:
        raise ValueError("P5D_FDR_COUNTS_SOURCE_DRIFT")
    if float(p5b["loo_logFC_pearson_min"])<.9972:
        raise ValueError("P5B_NINE_POOL_SENSITIVITY_CHANGED")
    return {
        "primary_donor_paired_model_n":57,
        "pooled_source_10x_library_pairs":9,
        "pooled_18_library_model_residual_df":8,
        "donor_vs_pool_model_common_tested_genes":p5d["common_genes"],
        "donor_significant_and_pool_significant_common":p5d["FDR_genes_shared"],
        "primary_significant_in_common":p5d["primary_FDR_genes_common"],
        "share_of_primary_common_FDR_retained":p5d["primary_FDR_retain_frac"],
        "pooled_vs_donor_model_gene_logFC_Pearson":p5d["pearson_logFC_vs_donor_model"],
        "nine_omit_one_pool_Pearson_min":p5b["loo_logFC_pearson_min"],
        "114_pseudobulk_samples_MDS1_treatment_single_variable_R2":r2,
        "114_pseudobulk_samples_MDS1_pool_single_variable_R2":md["first_axis_pool_r2"],
        "mds_r2_is_not_causal_fraction_explained":True,
        "pooled_9_library_pair_model_does_not_identify_drug_vs_batch":True,
    }


def fixture_test():
    result=full_rank_audit()
    assert result["new_rank_from_treatment_after_unrestricted_library_effects"]==0
    # Synthetic two-pool exact Kitagawa identity (and a pure mix case).
    p0=np.array([.25,.75]);p1=np.array([.25,.75])
    w0=np.array([.5,.5]);w1=np.array([.75,.25])
    within=sum((w0+w1)/2*(p1-p0))
    mix=sum((p0+p1)/2*(w1-w0))
    assert abs(within)<1e-12 and abs(mix-(-.125))<1e-12
    assert abs((w1@p1-w0@p0)-within-mix)<1e-12
    # A pure within-pool case.
    w1=w0.copy();p1=p0+.1
    within=sum((w0+w1)/2*(p1-p0))
    mix=sum((p0+p1)/2*(w1-w0))
    assert abs(within-.1)<1e-12 and abs(mix)<1e-12
    print("OCULAR_CONFOUNDING_RANK_AND_KITAGAWA_SELFTEST_PASS",flush=True)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--selection",type=Path)
    p.add_argument("--pool-loo",type=Path)
    p.add_argument("--pool-model",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    if a.self_test:
        fixture_test()
        return
    if not all((a.selection,a.pool_loo,a.pool_model,a.out)):
        p.error("Requires --selection --pool-loo --pool-model --out")
    sources=[a.selection,a.pool_loo,a.pool_model]
    vals=[json.loads(x.read_text(encoding="utf8")) for x in sources]
    alias=full_rank_audit()
    pool=selection_pool_decomposition(vals[0])
    models=inspect_models(vals[1],vals[2])
    sha={name:hashlib.sha256(path.read_bytes()).hexdigest() for name,path in
         zip(("P5e_selection","P5b_pool_LOO","P5d_library_unit"),sources)}
    report={
        "status":"PASS_AGGREGATE_IDENTIFIABILITY_AND_SELECTION_DECOMPOSITION",
        "dataset":"GSE313736","run_utc":datetime.now(timezone.utc).isoformat(),
        "original_run_ids":{"P5e":41091146021,"P5b":41085982185,"P5d":41090827791},
        "source_sha256":sha,
        "model_design_exact_rank":alias,
        "pool_level_fraction_decomposition":pool,
        "original_paired_gene_effect_model_diagnostics":models,
        "conditional_on_original_study_and_its_library_design":True,
        "no_randomized_library_cross_condition_replicates":True,
        "not_a_new_independent_donor_or_source_cohort":True,
        "causal_dex_effect_vs_original_library_nuisance":"NOT_IDENTIFIABLE_FROM_SOURCE",
        "cell_selection_bias_causal_contribution":"NOT_IDENTIFIABLE_FROM_AGGREGATE",
        "author80_gene_overlap":"NONBLOCKING_HOLD_NO_OUTREACH",
        "clinical_steroid_induced_IOP":"NOT_MEASURED",
    }
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"GSE313736_technical_confounding_design_and_pool_decomposition.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    lines=[
        "# GSE313736 technical confounding source audit",
        "",
        "**IDENTIFIABILITY HOLD.** Treatment-column rank is not incremented by ",
        "adding a DEX arm to 18 unrestricted original-library fixed effects.",
        "The source has 9 paired pools × 2 original treatment-specific 10x libraries.",
        "Pairing estimates a contrast under modelling assumptions; it does not",
        "make drug and arbitrary library-preparation effects separately estimable.",
        "",
        "## Real previously generated metadata: pooled versus paired fractions",
        "",
        f"- Global selection difference (DEX−Vehicle): {pool['pooled_BeamA_Dex_minus_Vehicle_percentage_points']:+.6f} percentage points",
        f"- Symmetric within-pool contribution: {pool['symmetric_within_pool_fraction_contribution_pp']:+.6f} pp",
        f"- Symmetric change in arm-dependent pool weights: {pool['symmetric_pool_arm_weight_mix_contribution_pp']:+.6f} pp",
        f"- Equal-pool mean within fraction delta: {pool['equal_pool_mean_delta_pp']:+.6f} pp",
        f"- Median of 59 paired donor fraction delta: {pool['paired_donor_BeamA_fraction_median_delta_pp']:+.6f} pp; positive/negative=29/30; two-sided sign p={pool['paired_sign_test_two_sided_p']}",
        "",
        "| Pool | Vehicle SNG | DEX SNG | Vehicle BeamA | DEX BeamA | Within-pool Δ (pp) |",
        "|---|---:|---:|---:|---:|---:|"
    ]
    for x in pool["per_pool"]:
        lines.append(f"| {x['pool']} | {x['Vehicle_SNG']:,} | {x['DEX_SNG']:,} | "
                     f"{x['Vehicle_BeamA_fraction']:.5f} | {x['DEX_BeamA_fraction']:.5f} | "
                     f"{x['within_pool_delta_percentage_points']:+.3f} |")
    lines+=["","## Gene-level and MDS evidence does not resolve aliasing",
       f"- Primary-versus-library gene model shared FDR: {models['donor_significant_and_pool_significant_common']:,}/{models['primary_significant_in_common']:,} (common original significance denominator)",
       f"- MDS first-axis treatment-only R²: {models['114_pseudobulk_samples_MDS1_treatment_single_variable_R2']:.6f} (not a causal variance estimate)",
       f"- Omit-one-pool log2FC minimum Pearson: {models['nine_omit_one_pool_Pearson_min']:.6f} (within same study)",
       "", "No genetic or clinical participant IDs, individual counts or raw sequencing were accessed.",
       "Technical and drug biological effects remain inseparable with unrestricted source-library effects.",
       "The observed pool-composition component is a descriptive accounting identity, NOT a fraction of DEX effect mediated by batch or cell selection."]
    (a.out/"GSE313736_technical_confounding_public_report.md").write_text(
        "\n".join(lines)+"\n",encoding="utf8")
    print("OCULAR_TECH_CONFOUNDING_REAL_PASS",json.dumps({
        "rank":alias,"composition_pp":{
            "pooled":pool["pooled_BeamA_Dex_minus_Vehicle_percentage_points"],
            "within_pool":pool["symmetric_within_pool_fraction_contribution_pp"],
            "pool_weight_mix":pool["symmetric_pool_arm_weight_mix_contribution_pp"],
            "donor_median":pool["paired_donor_BeamA_fraction_median_delta_pp"],
        },"models":models
    },ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
