#!/usr/bin/env python3
"""Integrate whole-locus G0022 SuSiE RSS + sensitivity + old QTL provenance.

Fail closed: a numerically converged posterior is not validated causal mapping.
"""
import argparse,csv,json
from pathlib import Path

ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
PILOT=ROOT.parent/"susie_g0022_pilot"
def load(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter="\t"))

def audit(root,pilot):
    summary=json.loads((root/"G0022_FULL_AIS_SUSIE_SUMMARY.json").read_text())
    reference=json.loads((root/"FULL_LOCUS_INPUT_QC.json").read_text())
    proximal=json.loads((root/"G0022_FULL_LOCUS_POSITIONAL_QTL_AUDIT_SUMMARY.json").read_text())
    sensitivity=load(root/"G0022_FULL_AIS_SENSITIVITY.tsv")
    old=load(pilot/"G0022_SUSIE_PILOT_EVIDENCE_GATES.tsv")
    if len(sensitivity)!=3 or len(old)!=5:
        raise RuntimeError("Sensitivity/pilot provenance is incomplete")
    if summary["status"]!="FULL_LOCUS_EXPLORATORY_ONLY" or reference["status"]!="FULL_G0022_SIGNED_LD_INPUT_QC_PASS":
        raise ValueError("Incomplete full locus data")
    if summary["n_snps"]!=reference["genotype_snp_count"] or reference["genotype_n"]!=504:
        raise ValueError("LD/genotype dimension mismatch")
    if summary["approximate_n_eff"]!=70474 or summary["n_eff_is_per_snp"]:
        raise ValueError("Unexpected sample size interpretation")
    if summary["final_causal_signal_validated"] or summary["final_causal_gene_validated"]:
        raise ValueError("Unjustified causal elevation")
    if proximal["new_AIS_QTL_coloc_completed"]!=0:
        raise RuntimeError("Stale molecular evidence provenance")
    tests={r["config"]:r for r in sensitivity}
    expected={"LOWER_N20000_L8","APPROX_NEFF_L3","APPROX_NEFF_L8_LD_SHRINK1PCT"}
    if set(tests)!=expected:raise RuntimeError("Missing required full locus sensitivity")
    stable=(all(r["converged"]=="TRUE" for r in sensitivity)
        and len({r["purity_filtered_cs"] for r in sensitivity})==1
        and len({r["cs_sizes"] for r in sensitivity})==1
        and len({r["max_pip_variant"] for r in sensitivity})==1
        and next(iter({r["max_pip_variant"] for r in sensitivity}))==summary["max_pip_snp"])
    local_ais=[r for r in old if r["trait"]=="AIS"]
    if len(local_ais)!=3:raise ValueError("Missing original 3 AIS local clump windows")
    records=[]
    for r in sensitivity:
        records.append({"analysis":"FULL_LOCUS_SUSIE_RSS",
           "configuration":r["config"],
           "n_ref":reference["genotype_n"],"n_snps":summary["n_snps"],
           "n_assumed":r["n_assumed"],
           "converged":r["converged"],
           "purity_filtered_cs":r["purity_filtered_cs"],
           "top_variant":r["max_pip_variant"],
           "top_pip":r["max_pip"],
           "causal_signal_validated":"NO",
           "multi_ancestry_replication":"NOT_TESTED",
           "colocalization_with_AIS":"NOT_TESTED",
           "interpretation":"EXPLORATORY_SENSITIVITY_ONLY"})
    with (root/"G0022_FULL_AIS_PUBLICATION_GATE.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,delimiter="\t",fieldnames=list(records[0]))
        w.writeheader();w.writerows(records)
    decision={"status":"EXPLORATORY_FULL_LOCUS_SUSIE_QC_ONLY",
       "whole_locus_converged":summary["converged"],
       "source_SNP_count":reference["source_harmonized_variants"],
       "full_locus_SNP_count":summary["n_snps"],
       "reference_n":504,
       "effective_N_study_approx":summary["approximate_n_eff"],
       "LD_matrix_rank_upper_bound":503,
       "LD_matrix_singular_due_to_n_ref_less_than_n_snps":True,
       "local_AIS_pilot_windows":len(local_ais),
       "full_locus_purity_filtered_cs":summary["n_credible_sets"],
       "full_locus_cs_sizes":summary["credible_set_sizes"],
       "top_pip_snp":summary["max_pip_snp"],
       "top_pip":summary["max_pip"],
       "sensitivity_configurations":len(sensitivity),
       "all_sensitivity_converged_and_same_CS_top_snp":bool(stable),
       "full_AIS_molecular_QTL_coloc_done":0,
       "old_BBJ_QTL_max_PP_H4_above_0_8":proximal["legacy_BBJ_coloc_h4_above_0_8"],
       "sRss_finiteLD_reference_correction":"NOT_APPLIED",
       "causal_signal_publication_ready":False,
       "causal_gene_publication_ready":False,
       "independent_genetic_replication_done":False,
       "next_needed":["PER_VARIANT_EFFECTIVE_N","LARGER_ANCESTRY_MATCHED_LD_REFERENCE",
          "FORMAL_FINITE_REFERENCE_LD_CORRECTION","QTL_FINE_MAPPING_AND_COLOCALIZATION",
          "INDEPENDENT_COHORT_REPLICATION"]}
    (root/"G0022_FULL_AIS_PUBLICATION_GATE_SUMMARY.json").write_text(json.dumps(decision,indent=2))
    print(json.dumps(decision,indent=2))
    return decision
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=ROOT)
    p.add_argument("--pilot",type=Path,default=PILOT)
    a=p.parse_args()
    audit(a.root,a.pilot)
