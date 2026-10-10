#!/usr/bin/env python3
"""Fail-closed science gate after two real Japanese alcohol locus SuSiE tests.

Never translate exploratory CS/PIP into conditional independence, alcohol
mediation, AIS stroke causality, drugs, or genomic 2,225 candidate shrinkage.
"""
import csv,json,math
from pathlib import Path
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/alcohol_conditional_susie_v1")
def records(file):
    with file.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def audit(root):
    qcs={g:json.loads((root/g/"input_qc.json").read_text()) for g in ("ALDH2","ADH1B")}
    d={g:json.loads((root/g/"ALCOHOL_SUSIE_DIAGNOSTIC_AND_GATE.json").read_text()) for g in qcs}
    for g in qcs:
        if qcs[g]["status"]!="ALCOHOL_REGIONAL_REAL_GWAS_EAS_LD_INPUT_PASS":
            raise ValueError("Source GWAS/LD provenance missing")
        if not qcs[g]["source_beta_SE_per_SNP_verified"] or not qcs[g]["LD_input_order_exact_variant_table"]:
            raise ValueError("GWAS beta SE and dosage LD order verification missing")
        if qcs[g]["reference_EAS_n"]!=504 or qcs[g]["variants"]!=d[g]["snps"]:
            raise ValueError("Input variant count not aligned to SuSiE diagnostic")
        if d[g]["reference_n"]!=504 or d[g]["final_finemapping_validated"] or d[g]["causal_alcohol_to_BP_AIS_mediation_estimated"]:
            raise ValueError("Causal status source discrepancy")
    filters=records(root/"ALCOHOL_ADH1B_ALDH2_GWAS_LD_MISMATCH_FILTER_SENSITIVITY.tsv")
    models=records(root/"ALCOHOL_ADH1B_SUSIE_MODEL_STABILITY.tsv")
    if len(filters)!=8 or len(models)!=6:raise ValueError("Incomplete filter and L sensitivity")
    for g in qcs:
        f=[x for x in filters if x["locus"]==g and x["scenario"]=="SOURCE_QC_DEFAULT"]
        if len(f)!=1 or abs(float(f[0]["LD_mismatch_s"])-float(d[g]["LD_mismatch_s"]))>1e-6:
            raise ValueError("Primary GWAS LD mismatch value not reproduced")
    if float(d["ALDH2"]["LD_mismatch_s"])<=.05 or not d["ALDH2"]["LD_gate"].startswith("BLOCKED"):
        raise ValueError("ALDH2 mismatch gate violated; no false causal promotion")
    if d["ALDH2"]["converged"] or d["ALDH2"]["n_credible_sets"]!=0:
        raise ValueError("ALDH2 should not run source-unmatched SuSiE")
    if d["ADH1B"]["LD_gate"]!="EXPLORATORY_MODEL_RUN_ALLOWED_NOT_CAUSAL" or not d["ADH1B"]["converged"]:
        raise ValueError("ADH1B no admissible exploratory fit")
    if any(x["causal_SNP_validated"]!="FALSE" or x["causal_gene_validated"]!="FALSE" for x in filters):
        raise ValueError("Filtered mismatch falsely made causal assertion")
    if any(x["biological_target_or_mediation_validated"]!="FALSE" or x["cross_model_independent_causal_signals_proven"]!="FALSE" for x in models):
        raise ValueError("SuSiE sensitivity model falsely validates biology")
    num=set(int(x["exploratory_95pct_credible_sets"]) for x in models)
    leads=set(x["top_pip_snp"] for x in models)
    if len(num)<2 or len(leads)<2:
        raise ValueError("Unexpected model stability result; audit interpretation")
    if any(x["converged"]!="TRUE" for x in models):
        raise ValueError("Unstable ADH1B convergence in modeled grid")
    if any(int(x["SNP_count"])<50 or not 0<=float(x["LD_mismatch_s"])<=1 for x in models):
        raise ValueError("Invalid ADH1B model diagnostics")
    result={
      "status":"ALDH2_BLOCKED_ADH1B_EXPLORATORY_UNSTABLE",
      "source":"Koyanagi_2024_Japanese_unstratified_daily_alcohol_full_MD5_verified",
      "reference":"1000G_EAS504_grch37_external_to_Japanese_GWAS",
      "ALDH2":{"source_input_variants":d["ALDH2"]["snps"],
        "LD_mismatch_s":d["ALDH2"]["LD_mismatch_s"],
        "source_p_numeric_underflow":qcs["ALDH2"]["source_p_zero_count"],
        "safety_result":"BLOCKED_GWAS_REFERENCE_LD_MISMATCH",
        "credible_sets_run":False,
        "independent_causal_effects_justified":0},
      "ADH1B":{"source_input_variants":d["ADH1B"]["snps"],
        "LD_mismatch_s":d["ADH1B"]["LD_mismatch_s"],
        "converged":True,
        "default_exploratory_credible_sets":d["ADH1B"]["n_credible_sets"],
        "L_by_filter_model_count":len(models),
        "exploratory_credible_set_counts":sorted(num),
        "top_lead_variant_varies_with_filtering":len(leads)>1,
        "top_pip_variants_by_filter":sorted(leads),
        "causal_credible_sets_validated":0},
      "mismatch_s_gate_is_inhouse_conservative_threshold_not_universal_community_standard":True,
      "per_variant_GWAS_N_metadata_available_but_rss_uses_locus_median":True,
      "different_Japanese_and_mixed_EAS_LD_reference_populations":True,
      "GWAS_source_p_below_precision_does_not_equal_true_p_zero":True,
      "single_region_6CS_are_not_six_independent_causal_variants":True,
      "independent_ischemic_stroke_replication":False,
      "exclusion_restriction_alcohol_specific_instruments_verified":False,
      "horizontal_pleiotropy_absent":False,
      "two_sample_cohort_independence_verified":False,
      "MR_IVW_or_Egger_calculated":False,
      "ALDH2_enzymatic_stroke_mediation_claim":False,
      "maintain_IS_2225_positional_gene_discovery_universe":True,
      "next_valid_gate":"ancestry_and_cohort_matched_large_LD_or_sufficient_snp_cohort_statistics_then_multisignal_conditional_finemapping_and_pleiotropy_assessment"}
    (root/"IS_ALCOHOL_ADH1B_ALDH2_CONDITIONAL_SCIENTIFIC_GATE.json").write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    return result
if __name__=="__main__":audit(ROOT)
