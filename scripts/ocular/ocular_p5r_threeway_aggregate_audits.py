#!/usr/bin/env python3
"""P5r: three parallel aggregate-only, non-causal TM target evidence audits.

Modes:
- thresholds: P5p post-hoc cutpoint sensitivity over all 96 candidates.
- timecourse: timepoint significance vs raw sign for the exact 12 candidates.
- consensus: all 96 P5o donor/pool descriptives, P5p abundance, P5m external FDR.

No raw reads, donor identifiers, clinical outcomes, or independent replication.
"""
import argparse
import csv
import json
import math
from pathlib import Path
from collections import Counter

LOW = ("ERFE", "CHAC1", "PARM1")
ADEQUATE = ("ADAMTS15", "EGR3", "LPAR1", "ABHD3", "DDIT4L",
            "TFPI", "STK17B", "RAB32", "SDC1")
TARGETS = set(LOW + ADEQUATE)

def number(v):
    try:
        f = float(v)
        return f if math.isfinite(f) else None
    except (TypeError, ValueError):
        return None

def required_number(row, field):
    value = number(row.get(field))
    if value is None:
        raise ValueError("NONFINITE_REQUIRED_FIELD " + field + " " + row.get("gene_id", ""))
    return value

def read96(path):
    with open(path, "r", newline="", encoding="utf8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "gene_id" not in reader.fieldnames or "gene_symbol" not in reader.fieldnames:
            raise ValueError("MISSING_GENE_ID_OR_SYMBOL " + str(path))
        data = list(reader)
    if len(data) != 96 or len({r["gene_id"] for r in data}) != 96:
        raise ValueError("SOURCE_MUST_HAVE_96_UNIQUE_IDS " + str(path))
    if len({r["gene_symbol"] for r in data}) != 96:
        raise ValueError("DUPLICATED_SYMBOLS_REQUIRE_MANUAL_REVIEW")
    return data

def target_check(data):
    items = [r for r in data if r["gene_symbol"] in TARGETS]
    if len(items) != 12 or {r["gene_symbol"] for r in items} != TARGETS:
        raise ValueError("P5P_12_TARGET_SET_CHANGED")
    return items

def save_tsv(path, rows):
    if not rows:
        raise ValueError("EMPTY_PUBLIC_RESULT")
    with path.open("w", newline="", encoding="utf8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

def save_json(path, content):
    path.write_text(json.dumps(content, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                    encoding="utf8")

def thresholds(p5p, out):
    dat = read96(p5p)
    target_check(dat)
    base = []
    for r in dat:
        primary = required_number(r, "primary_BeamA_500nM_48h_log2fc")
        external = required_number(r, "external_TM_100nM_48h_log2fc")
        pool = required_number(r, "conservative_poolunit_log2fc")
        fdr = required_number(r, "conservative_poolunit_FDR")
        logcpm = required_number(r, "primary_mean_logCPM")
        if not primary * external < 0:
            raise ValueError("SOURCE_NOT_CROSS_STUDY_OPPOSITE " + r["gene_id"])
        base.append((r, primary, external, pool, fdr, logcpm))
    grid = []
    baseline_names = set()
    for abundance in (0, 1, 2, 3, 4, 5):
        for min_effect in (0.25, 0.5, 0.75, 1.0):
            names = sorted(r["gene_symbol"] for r, p, e, pool, fdr, cp in base
                if cp >= abundance and fdr < 0.05 and p*pool > 0
                and min(abs(p), abs(e), abs(pool)) >= min_effect)
            if abundance == 2 and min_effect == 0.5:
                baseline_names = set(names)
            grid.append({"logCPM_cutoff": abundance, "min_abs_log2fc_all_three": min_effect,
                         "conservative_pool_FDR_cutoff": 0.05,
                         "n_genes": len(names), "genes": ",".join(names),
                         "posthoc_screen_not_hypothesis_test": True})
    if baseline_names != set(ADEQUATE):
        raise ValueError("P5P_NINE_GENE_BASELINE_DRIFT " + repr(baseline_names))
    detailed = []
    for r, p, e, pool, fdr, cp in base:
        if r["gene_symbol"] not in TARGETS:
            continue
        detailed.append({"gene_id": r["gene_id"], "gene_symbol": r["gene_symbol"],
            "tier": "LOW_EXPRESSION" if r["gene_symbol"] in LOW else "ADEQUATE_EXPRESSION",
            "primary_mean_logCPM": cp, "primary_log2fc": p,
            "external_100nM_48h_log2fc": e, "poolunit_log2fc": pool,
            "poolunit_FDR": fdr,
            "passes_logCPM3_effect05": cp>=3 and min(abs(p),abs(e),abs(pool))>=.5 and fdr<.05 and p*pool>0,
            "passes_logCPM2_effect075": cp>=2 and min(abs(p),abs(e),abs(pool))>=.75 and fdr<.05 and p*pool>0})
    save_tsv(out/"P5r_threshold_cutpoint_grid.tsv", grid)
    save_tsv(out/"P5r_12_threshold_sensitivity.tsv", detailed)
    result = {"status":"PASS_AGGREGATE_POSTHOC_THRESHOLD_SENSITIVITY",
      "source":"P5p 96 gene-level public aggregate",
      "grid_combinations":len(grid), "n_input":len(dat),
      "baseline_logCPM2_min_abs05_names":sorted(baseline_names),
      "baseline_nine_reproduced":True,
      "note":"Changing cutpoints is exploratory robustness description, NOT an independent p-value or preregistered selection.",
      "dose_causality_established":False, "clinical_IOP_prediction_allowed":False}
    save_json(out/"P5r_threshold_summary.json",result)
    return result

def timecourse(p5m, p5p, out):
    tm = read96(p5m)
    p = read96(p5p)
    if {r["gene_id"] for r in tm}!={r["gene_id"] for r in p}:
        raise ValueError("P5M_P5P_GENE_JOIN_DRIFT")
    target_check(p)
    tm_by = {r["gene_id"]:r for r in tm}
    output = []
    for x in sorted((r for r in p if r["gene_symbol"] in TARGETS),key=lambda r:r["gene_symbol"]):
        t=tm_by[x["gene_id"]]
        if t["gene_symbol"]!=x["gene_symbol"]:
            raise ValueError("SYMBOL_GENE_ID_MISMATCH")
        external48=required_number(t,"external_TM_100nM_48h_log2fc")
        external48_fdr=required_number(t,"external_TM_100nM_48h_FDR")
        primary=required_number(t,"primary_BeamA_500nM_48h_log2fc")
        if primary*external48>=0 or external48_fdr>=.05:
            raise ValueError("NOT_QUALIFIED_48H_OPPOSITE")
        row={"gene_id":x["gene_id"],"gene_symbol":x["gene_symbol"],
             "tier":"LOW_EXPRESSION" if x["gene_symbol"] in LOW else "ADEQUATE_EXPRESSION",
             "primary_500nM_48h_log2fc":primary,
             "external_100nM_48h_log2fc":external48,
             "external_100nM_48h_FDR":external48_fdr}
        early_sig=0
        raw_flip=0
        tested=0
        for time in ("1h","6h"):
            logfc=number(t.get("external_TM_100nM_"+time+"_log2fc"))
            padj=number(t.get("external_TM_100nM_"+time+"_FDR"))
            valid=logfc is not None
            significant=valid and padj is not None and padj<.05
            flip=valid and logfc*external48<0
            row[time+"_available"]=valid
            row[time+"_log2fc"]=logfc if valid else ""
            row[time+"_FDR"]=padj if padj is not None else ""
            row[time+"_FDR_significant"]=bool(significant)
            row[time+"_raw_sign_opposite_48h"]=bool(flip)
            tested+=int(valid)
            early_sig+=int(significant)
            raw_flip+=int(flip)
        row["early_significant_timepoints"]=early_sig
        row["early_raw_sign_flips"]=raw_flip
        row["early_timepoints_with_log2fc"]=tested
        row["formal_time_interaction_available"]=False
        output.append(row)
    if len(output)!=12:
        raise ValueError("MISSING_12_TIMECOURSE_ROWS")
    save_tsv(out/"P5r_12_gene_external_timecourse.tsv",output)
    summary={"status":"PASS_PUBLICATION_TIMECOURSE_DESCRIPTIVE",
      "n_targets":12,
      "one_hour_fdr_significant":sum(x["1h_FDR_significant"] for x in output),
      "six_hour_fdr_significant":sum(x["6h_FDR_significant"] for x in output),
      "with_any_early_significant":sum(x["early_significant_timepoints"]>0 for x in output),
      "with_any_raw_sign_flip":sum(x["early_raw_sign_flips"]>0 for x in output),
      "external_is_same_study_one_dose":True,
      "formal_time_interaction_estimated":False,
      "causal_100_vs_500nM_reversal_established":False,
      "note":"Early non-significant direction fluctuations are NOT validated temporal reversals. Published contrasts lack their sampling covariance."}
    save_json(out/"P5r_timecourse_summary.json",summary)
    return summary

def boolfield(v):
    if v in ("TRUE","True","true","1"):return True
    if v in ("FALSE","False","false","0"):return False
    raise ValueError("NON_BOOLEAN_VALUE "+str(v))

def consensus(p5p,p5o,out):
    p=read96(p5p)
    o=read96(p5o)
    target_check(p)
    om={r["gene_id"]:r for r in o}
    if set(om)!={r["gene_id"] for r in p}:
        raise ValueError("P5P_P5O_96_GENE_ID_MISMATCH")
    result=[]
    for r in p:
        a=om[r["gene_id"]]
        if a["gene_symbol"]!=r["gene_symbol"]:
            raise ValueError("P5O_SYMBOL_MISMATCH")
        sign_same=int(required_number(a,"paired_57_sign_same"))
        pool_sign=int(required_number(a,"pool_sign_same"))
        loo_sign=int(required_number(a,"leave_one_pool_out_mean_sign_same"))
        ci_lo=required_number(a,"pool_bootstrap_95_lo")
        ci_hi=required_number(a,"pool_bootstrap_95_hi")
        interval_excludes=boolfield(a["pool_bootstrap_95_excludes_zero"])
        if sign_same<0 or sign_same>57 or pool_sign<0 or pool_sign>9 or loo_sign<0 or loo_sign>9:
            raise ValueError("INVALID_SOURCE_POOL_DONOR_COUNT")
        if interval_excludes!=(ci_lo*ci_hi>0):
            raise ValueError("BOOTSTRAP_FLAG_BOUNDARY_INCONSISTENT")
        pfc=required_number(r,"primary_BeamA_500nM_48h_log2fc")
        efc=required_number(r,"external_TM_100nM_48h_log2fc")
        poolfc=required_number(r,"conservative_poolunit_log2fc")
        if pfc*efc>=0:raise ValueError("NOT_CROSS_STUDY_DISCORDANT")
        tier=("LOW_EXP_LARGE" if r["gene_symbol"] in LOW
              else "ADEQUATE_EXP" if r["gene_symbol"] in ADEQUATE
              else "OTHER_96")
        result.append({"gene_id":r["gene_id"],"gene_symbol":r["gene_symbol"],
          "tier":tier,"logCPM":required_number(r,"primary_mean_logCPM"),
          "primary_log2fc":pfc,"external_log2fc":efc,"poolunit_log2fc":poolfc,
          "poolunit_FDR":required_number(r,"conservative_poolunit_FDR"),
          "donor_signed_n":sign_same,"donor_total":57,
          "pool_signed_n":pool_sign,"pool_total":9,
          "omitted_pool_mean_signed_n":loo_sign,
          "pool_bootstrap_lo":ci_lo,"pool_bootstrap_hi":ci_hi,
          "pool_bootstrap_excludes_zero":interval_excludes,
          "descriptive_only_no_ql_inference":True})
    if sum(r["pool_signed_n"]==9 for r in result)!=33:
        raise ValueError("P5O_33_ALL_9_POOL_SIGN_CHECK_CHANGED")
    if sum(r["omitted_pool_mean_signed_n"]==9 for r in result)!=96:
        raise ValueError("P5O_96_ALL_9_LOO_MEANS_CHANGED")
    if sum(r["pool_bootstrap_excludes_zero"] for r in result)!=75:
        raise ValueError("P5O_75_CI_EXCLUDE_ZERO_CHANGED")
    result.sort(key=lambda x:(0 if x["tier"]=="ADEQUATE_EXP" else
                              1 if x["tier"]=="LOW_EXP_LARGE" else 2, x["gene_symbol"]))
    save_tsv(out/"P5r_96_aggregate_evidence_reconciliation.tsv",result)
    subset=[r for r in result if r["tier"]!="OTHER_96"]
    save_tsv(out/"P5r_12_candidate_aggregate_evidence.tsv",subset)
    summary={"status":"PASS_AGGREGATE_96_GENE_CROSS_SOURCE_RECONCILIATION",
      "n_genes":96,"n_target":12,
      "all_nine_source_pools_signed":33,
      "all_nine_omitted_pool_means_signed":96,
      "pool_bootstrap_intervals_exclude_zero":75,
      "candidate_12_all_9_pool_signs":sum(x["pool_signed_n"]==9 for x in subset),
      "candidate_12_all_9_loo_mean_signs":sum(x["omitted_pool_mean_signed_n"]==9 for x in subset),
      "candidate_12_pool_bootstrap_excludes_zero":sum(x["pool_bootstrap_excludes_zero"] for x in subset),
      "not_formal_9_edgeR_QL_refits":True,
      "not_independent_replication":True,
      "note":"Original P5o pool bootstrap and LOO means are descriptive; P5q formal edgeR refits must be interpreted separately."}
    save_json(out/"P5r_consensus_summary.json",summary)
    return summary

def self_test():
    assert number("0.05")==.05 and number("") is None and number("NaN") is None
    assert boolfield("TRUE") and not boolfield("False")
    assert len(TARGETS)==12 and len(set(ADEQUATE))==9
    print("P5R_SELF_TEST_PASS")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mode",choices=("thresholds","timecourse","consensus"))
    ap.add_argument("--p5p",type=Path)
    ap.add_argument("--p5m",type=Path)
    ap.add_argument("--p5o",type=Path)
    ap.add_argument("--output-dir",type=Path)
    ap.add_argument("--self-test",action="store_true")
    args=ap.parse_args()
    if args.self_test:
        self_test()
        return
    if not args.mode or not args.output_dir or not args.p5p:
        ap.error("--mode, --output-dir and --p5p are required")
    args.output_dir.mkdir(parents=True,exist_ok=True)
    if args.mode=="thresholds":
        report=thresholds(args.p5p,args.output_dir)
    elif args.mode=="timecourse":
        if not args.p5m:ap.error("--p5m required for timecourse")
        report=timecourse(args.p5m,args.p5p,args.output_dir)
    else:
        if not args.p5o:ap.error("--p5o required for consensus")
        report=consensus(args.p5p,args.p5o,args.output_dir)
    print("P5R_RESULT",args.mode,json.dumps(report,ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()
