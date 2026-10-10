#!/usr/bin/env python3
"""Gate isolated susieR 0.16.6 source-based finite reference LD experiments.

This emits EXPLORATORY method sensitivity only even if model converges
and maximum PIP=1. False/absent reference reliability and larger cohort
requirements cannot be silently promoted to a causal conclusion.
"""
import csv,hashlib,json,math
from pathlib import Path
DEFAULT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/susie_rss_finite_ref_sandbox_v1")
SOURCE=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/alcohol_conditional_susie_v1")
LOCI=("ADH1B","ALDH2")
MODES=("FINITE_REF_504","FINITE_REF_504_EB_MISMATCH")
def build(root=DEFAULT,source=SOURCE):
    package=root/"packages"/"susieR_0.16.6.tar.gz"
    dep=root/"packages"/"cpp11armadillo_0.5.4.tar.gz"
    if not all(p.is_file() for p in (package,dep)):
        raise ValueError("Isolated package source pinning missing")
    installation={name:hashlib.sha256(path.read_bytes()).hexdigest()
       for name,path in (("susieR_0.16.6",package),("cpp11armadillo_0.5.4",dep))}
    rows=[]
    for locus in LOCI:
        src=json.loads((source/locus/"ALCOHOL_SUSIE_DIAGNOSTIC_AND_GATE.json").read_text())
        qc=json.loads((source/locus/"input_qc.json").read_text())
        if qc["status"]!="ALCOHOL_REGIONAL_REAL_GWAS_EAS_LD_INPUT_PASS" or qc["reference_EAS_n"]!=504:
            raise ValueError("Original GWAS signed LD QC missing")
        for mode in MODES:
            p=root/"models"/locus/f"{mode}.qc.json"
            data=json.loads(p.read_text())
            q=data["diagnostics"]
            if not data["finished_without_error"] or not data["converged"]:
                raise ValueError("Method experiment failed to converge; requires review")
            if data["susieR_version"]!="0.16.6" or data["reference_actual_B"]!=504 or not data["R_finite_enabled"]:
                raise ValueError("Unexpected susieR package version or finite-LD model")
            if data["R_mismatch"]!=("eb" if mode==MODES[1] else "none"):
                raise ValueError("Wrong LD mismatch option")
            if data["original_snp_count"]!=src["snps"] or data["source_median_n"]!=src["median_sample_n"]:
                raise ValueError("Variant set/study sample N changed relative to frozen baseline")
            if not data["source_alleles_signed_and_order_QC"] or not data["finite_LD_reference_correction_sandbox_only"]:
                raise ValueError("Missing signed dosage input QC")
            if not q.get("diagnostics_provided") or not math.isfinite(float(q.get("r_over_B",float("nan")))):
                raise ValueError("New finite-LD correction diagnostics absent")
            if q.get("B")!=504 or not isinstance(q.get("R_reliability_flag"),bool):
                raise ValueError("Missing or invalid reliability flag")
            if not isinstance(data["credible_sets_exploratory"],int) or data["credible_sets_exploratory"]<0:
                raise ValueError("Missing exploratory credible set count")
            if not isinstance(q.get("R_sensitivity_flag"),bool):
                raise ValueError("Missing or invalid sensitivity flag")
            if mode==MODES[1] and (q.get("B_corrected") is None or q.get("lambda_bias") is None):
                raise ValueError("EB-specific uncertainty correction source parameters not recorded")
            warninglist=data.get("warnings",[])
            if not warninglist:raise ValueError("Missing source QC warnings unexpectedly")
            rows.append({
              "locus":locus,"model":mode,
              "software":"susieR_0.16.6_isolated",
              "SNP_n":data["original_snp_count"],
              "original_GWAS_ld_mismatch_s":data["original_no_correction_s"],
              "reference_n":504,"GWAS_median_n":data["source_median_n"],
              "R_finite":504,
              "R_mismatch_mode":data["R_mismatch"],
              "converged":data["converged"],
              "exploratory_95pct_CS_n":data["credible_sets_exploratory"],
              "top_variant":data["top_variant"],
              "top_PIP":data["max_PIP"],
              "R_reliability_flag":q["R_reliability_flag"],
              "R_sensitivity_flag":q["R_sensitivity_flag"],
              "R_r_over_B":q["r_over_B"],
              "reference_penalty_median":q.get("penalty_median"),
              "reference_penalty_max":q.get("penalty_max"),
              "EB_lambda_bias":q.get("lambda_bias",""),
              "EB_B_corrected":q.get("B_corrected",""),
              "warnings_count":len(warninglist),
              "credible_sets_scientifically_validated":0,
              "causal_IS_gene_validated":False,
              "alcohol_BP_AIS_causal_mediation_estimated":False,
              "reported_as_new_independent_causal_signals":False,
              "publication_ready":False})
    with (root/"IS_ALCOHOL_FINITE_REF_LD_CORRECTION_MODEL_COMPARISON.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
    status={
      "status":"EXPERIMENTAL_FINITE_REFERENCE_SUSIE_RUNS_PASS_BUT_RELIABILITY_FAIL",
      "software_isolated_and_pinned":installation,
      "server_default_susieR_package_unmodified":True,
      "total_models":4,
      "ADH1B_models":{x["model"]:x["exploratory_95pct_CS_n"] for x in rows if x["locus"]=="ADH1B"},
      "ALDH2_models":{x["model"]:x["exploratory_95pct_CS_n"] for x in rows if x["locus"]=="ALDH2"},
      "reliability_warning_all_4_models":all(x["R_reliability_flag"] for x in rows),
      "reliability_warning_count":sum(bool(x["R_reliability_flag"]) for x in rows),
      "sensitivity_warning_all_4_models":all(x["R_sensitivity_flag"] for x in rows),
      "sensitivity_warning_count":sum(bool(x["R_sensitivity_flag"]) for x in rows),
      "ALDH2_reference_penalty_median_in_original_reference":max(x["reference_penalty_median"] for x in rows if x["locus"]=="ALDH2"),
      "EB_extra_mismatch_lambda_estimates":{x["locus"]:x["EB_lambda_bias"] for x in rows if x["model"]==MODES[1]},
      "top_PIP_is_not_direct_causal_evidence":True,
      "ALDH2_prior_BLOCKED_gate_lifted":False,
      "ADH1B_causal_variant_verified":False,
      "valid_independent_alcohol_MR_instrument_set":False,
      "Koyanagi_BBJ_sample_overlap_resolved":False,
      "corresponding_independent_stroke_cohort_validated":False,
      "final_causal_variants":0,
      "final_causal_genes":0,
      "clinical_alcohol_mediation_effect_reported":False,
      "maintain_original_IS_positional_genes":2225,
      "need_adequate_Japanese_cohort_ld":True,
      "candidate_reference_jMorp_61KJPN_contains_individual_ld_confirmed":False}
    (root/"IS_ALCOHOL_FINITE_REF_LD_CORRECTION_SCIENCE_GATE.json").write_text(json.dumps(status,indent=2))
    print(json.dumps(status,indent=2))
    return status
if __name__=="__main__":build()
