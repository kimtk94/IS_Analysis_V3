#!/usr/bin/env python3
"""Fail-closed final G0022 mismatch diagnostic ledger; NOT causal inference."""
import csv,json,argparse,math
from pathlib import Path

ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/susie_g0022_pilot")
def tsv(p):
    with p.open() as h:return list(csv.DictReader(h,delimiter="\t"))
def build(root):
    neff=tsv(root/"G0022_GIGASTROKE_EAS_EFFECTIVE_N_APPROX.tsv")
    susie=tsv(root/"effective_n_sensitivity/G0022_SUSIE_EFFECTIVE_N_SENSITIVITY.tsv")
    kriging=tsv(root/"G0022_KRIGING_OUTLIER_SUMMARY.tsv")
    filt=tsv(root/"G0022_LD_MISMATCH_FILTER_SENSITIVITY.tsv")
    diagnostic=tsv(root/"matched_as_ais_ld_sensitivity/G0022_SHARED_SNP_AS_AIS_LD_DIAGNOSTIC.tsv")
    random=tsv(root/"matched_as_ais_ld_sensitivity/G0022_AS_MATCHED_SNP_VS_RANDOM_DROPOUT.tsv")
    source_audit=json.loads((root/"G0022_AS_ONLY_SNP_CANONICAL_SOURCE_SUMMARY.json").read_text())
    if source_audit.get("as_only_marker_rows")!=78 or source_audit.get("source_presence_status",{}).get("NO_AIS_CANONICAL_POSITION")!=78:
        raise RuntimeError("AS-only canonical source audit is not complete")
    gates=tsv(root/"G0022_SUSIE_PILOT_EVIDENCE_GATES.tsv")
    if [len(neff),len(susie),len(kriging),len(filt),len(diagnostic),len(random),len(gates)] != [2,5,5,15,4,2,5]:
        raise RuntimeError("Incomplete research evidence inputs")
    if {x["trait"] for x in neff}!={"AS","AIS"}:
        raise ValueError("Incorrect population provenance")
    ns={x["trait"]:int(x["working_n_rounded"]) for x in neff}
    if ns!={"AS":98294,"AIS":70474}:
        raise ValueError("Unexpected EAS case/control effective N approximations")
    g={x["signal"]:x for x in gates}
    if len(g)!=5 or any(x["claim"]!="NOT_FINE_MAPPING_VALIDATED" for x in gates):
        raise RuntimeError("Original safety gates missing")
    kr={r["signal"]:r for r in kriging}
    runs={r["signal"]:r for r in susie}
    if set(g)!=set(kr) or set(g)!=set(runs):
        raise RuntimeError("G0022 pilot study mismatched")
    for r in runs.values():
        trait=r["trait"]
        if int(r["approximate_n_eff"])!=ns[trait]:
            raise RuntimeError("Approximate effective-N mismatch")
    random_by={r["window"]:r for r in random}
    shared={}
    for r in diagnostic:
        shared[(r["window"],r["trait"])]=r
    output=[]
    for name in sorted(g):
        gate=g[name]
        krg=kr[name]
        run=runs[name]
        is_as=run["trait"]=="AS"
        paired=random_by.get(name) if is_as else None
        if is_as and (not paired or (name,"AS") not in shared or (name,"AIS") not in shared):
            raise ValueError("Missing matched AS/AIS LD comparison")
        if is_as and float(paired["matched_study_shared_s"])>=0.05:
            raise RuntimeError("Unexpected shared SNP LD discrepancy")
        if is_as and int(paired["random_subset_s_le_matched"])>0:
            raise RuntimeError("Random SNP deletion not a clean negative control")
        if int(krg["flagged_loglr2_absz2"])<0:
            raise RuntimeError("Invalid Kriging results")
        row={
          "signal":name,"trait":run["trait"],
          "prior_gate":gate["scientific_gate"],
          "approximate_study_n_eff":run["approximate_n_eff"],
          "effective_n_pilot_converged":run["converged"],
          "effective_n_purity_cs_count":run["purity_filtered_cs_count"],
          "effective_n_top_pip":run["top_pip"],
          "kriging_flagged_isolated_flip_count":krg["flagged_loglr2_absz2"],
          "as_full_s":paired["full_s"] if paired else "",
          "as_shared_s":paired["matched_study_shared_s"] if paired else "",
          "as_shared_variant_n":paired["shared_n"] if paired else "",
          "as_random_s_min":paired["random_dropout_s_min"] if paired else "",
          "as_random_s_median":paired["random_dropout_s_median"] if paired else "",
          "as_n_random_deletions":paired["random_dropout_replicates"] if paired else "",
          "overall_status":"BLOCKED_AS_SPECIFIC_MARKER_QC" if is_as
                else "EXPLORATORY_AIS_NOT_FINE_MAPPING_VALIDATED",
          "causal_gene_status":"UNRESOLVED"}
        output.append(row)
    out=root/"G0022_DIAGNOSTIC_DECISION_V2.tsv"
    with out.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(output[0]),delimiter="\t")
        w.writeheader();w.writerows(output)
    summary={
      "windows_total":5,
      "case_control_effective_n_approx":ns,
      "effective_n_converged":sum(x["effective_n_pilot_converged"]=="TRUE" for x in output),
      "suspected_isolated_flips":sum(int(x["kriging_flagged_isolated_flip_count"]) for x in output),
      "AS_only_source_variants_missing_in_AIS_canonical":source_audit["source_presence_status"]["NO_AIS_CANONICAL_POSITION"],
      "AS_only_source_variants_gws":source_audit["as_only_gws_p_le_5e8"],
      "AS_windows_blocked":sum(x["overall_status"]=="BLOCKED_AS_SPECIFIC_MARKER_QC" for x in output),
      "AIS_windows_exploratory":sum(x["overall_status"].startswith("EXPLORATORY_AIS") for x in output),
      "random_negative_control_replicates_per_AS_window":20,
      "high_AS_LD_mismatch_on_matched_variants":sum(float(r["estimated_s"])>.05 for r in diagnostic if r["trait"]=="AS"),
      "root_cause_interpretation":"AS-only GWAS/reference SNP inclusion is associated with mismatch. Specific QC or sampling cause not established.",
      "fine_mapping_validated":0,
      "independent_replication":0,
      "scientific_warning":"No per-SNP effective N; 504-person reference; 250kb pilot; related AS/AIS samples"}
    (root/"G0022_DIAGNOSTIC_DECISION_V2_SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    return summary
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=ROOT)
    a=p.parse_args()
    build(a.root)
