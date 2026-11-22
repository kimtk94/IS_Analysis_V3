#!/usr/bin/env python3
"""Replay 15 BBJ exploratory ABF panels from SOURCE-DERIVED LOGSUM statistics.

GitHub Actions does NOT contain original SNP-level Drive tables or execute R coloc.
The 15 original source TSVs were actually read outside GitHub (2026-11-22);
each group retains source SHA256. The 225 scenarios below independently re-evaluate
coloc::combine.abf mathematics from sufficient log-sums (not study validity).
"""
from __future__ import annotations
import argparse
import json
import math
import re
from pathlib import Path

SCHEMA = "is-bbj-legacy-P01-P15-coloc-sufficient-20261122-v1"
N_FACTORS = [0.75, 1, 1.25]
P12_VALUES = [1e-6, 3e-6, 1e-5, 3e-5, 1e-4]

def check(p, error):
    if not p:
        raise ValueError(error)

def logsumexp(values):
    pivot = max(values)
    return pivot + math.log(math.fsum(math.exp(v - pivot) for v in values))

def posterior(s, p12, p1=1e-4, p2=1e-4):
    S1, S2, S12 = s
    check(S12 < S1 + S2, "H3 two-distinct-SNP exclusion violated")
    l3 = math.log(p1*p2) + S1 + S2 + math.log1p(-math.exp(S12-S1-S2))
    values = [0, math.log(p1)+S1, math.log(p2)+S2, l3, math.log(p12)+S12]
    denom = logsumexp(values)
    pp = [math.exp(v - denom) for v in values]
    check(abs(math.fsum(pp)-1) < 1e-9, "posterior normalization failed")
    return pp

def reconcile(source, shard=0, shards=1):
    check(source.get("schema") == SCHEMA, "source model/schema mismatch")
    check(source.get("source_run_real_local") is True
          and source.get("ci_raw_source_available") is False
          and source.get("study_N_and_alleles_approved") is False
          and source.get("statistical_status") == "EXPLORATORY_HOLD",
          "scientific provenance falsely upgraded")
    check(source.get("p12_values") == P12_VALUES
          and source.get("N_factor_values") == N_FACTORS, "prior/N sensitivity model changed")
    panels = source.get("panels")
    check(isinstance(panels,list) and len(panels) == 15, "expected 15 source-panel records")
    check(0 <= shard < shards <= 15, "invalid shard")
    rows = []
    for i, raw in enumerate(panels):
        check(isinstance(raw,list) and len(raw) == 10, "invalid sufficient statistic row")
        key,gene,tissue,n,source_bytes,sha,source_N,pal,gap,suff = raw
        check(key == f"P{i+1:02d}" and gene and tissue, "panel order/identity changed")
        check(type(n) is int and n > 200 and type(source_bytes) is int and source_bytes > 10000
              and type(source_N) is int and source_N > 0
              and type(pal) is int and 0 <= pal <= n
              and type(gap) is int and 0 <= gap <= n
              and isinstance(sha,str) and re.fullmatch("[a-f0-9]{64}",sha),
              "invalid source hash or original counts")
        check(isinstance(suff,list) and len(suff) == 3, "missing N scenarios")
        for triple in suff:
            check(isinstance(triple,list) and len(triple) == 3 and all(
                type(v) in (int,float) and math.isfinite(v) for v in triple),
                "invalid original logsum model")
        check(max(abs(x[0]-suff[0][0]) for x in suff) < 1e-8,
              "GWAS logsum changed when only QTL hypothetical N changed")
        if i % shards != shard:
            continue
        grid=[]
        for nfactor, triple in zip(N_FACTORS,suff):
            for p12 in P12_VALUES:
                pp = posterior(triple,p12)
                grid.append({
                    "hypothetical_qtl_N_factor": nfactor,
                    "p12": p12,
                    "H0": pp[0], "H1": pp[1], "H2": pp[2],
                    "H3": pp[3], "H4": pp[4],
                    "H4_given_H3_or_H4": pp[4]/(pp[3]+pp[4]) if pp[3]+pp[4] else None
                })
        baseline=next(x for x in grid if x["hypothetical_qtl_N_factor"]==1
                      and x["p12"]==1e-5)
        rows.append({
            "panel": key, "gene": gene, "tissue": tissue, "original_source_SNP_rows": n,
            "source_sha256": sha, "archived_scalar_QTL_N_NOT_APPROVED": source_N,
            "palindromic_input_rows": pal, "abs_cross_study_MAF_gap_gt_0p2_rows": gap,
            "baseline": baseline, "sensitivity_grid": grid,
            "min_H4": min(x["H4"] for x in grid),
            "max_H4": max(x["H4"] for x in grid),
            "scientific_status": "EXPLORATORY_HOLD",
            "new_independent_stroke_replication": False
        })
    return {
        "schema": "IS-legacy-15-assay-source-sufficient-shard-v1",
        "technical_status": "PASS_ARCHIVED_SUFFICIENT_MATH",
        "scientific_status": "EXPLORATORY_HOLD",
        "raw_source_files_re_read_in_ci": False,
        "R_coloc_executed": False,
        "P3A_36_group_executed": False,
        "shard": shard, "shards": shards,
        "panel_count": len(rows),
        "scenario_count": sum(len(x["sensitivity_grid"]) for x in rows),
        "panels": rows,
    }

def merge(files):
    check(len(files) == 3,"exactly three shard files required")
    all_panels=[]
    for i, path in enumerate(files):
        obj=json.loads(path.read_text(encoding="utf-8"))
        check(obj.get("schema")=="IS-legacy-15-assay-source-sufficient-shard-v1"
              and obj.get("shard")==i and obj.get("shards")==3
              and obj.get("scientific_status")=="EXPLORATORY_HOLD"
              and obj.get("raw_source_files_re_read_in_ci") is False,
              "unverified or incorrectly ordered shard")
        all_panels.extend(obj["panels"])
    all_panels.sort(key=lambda x:x["panel"])
    check([x["panel"] for x in all_panels]==[f"P{i:02d}" for i in range(1,16)],
          "panel IDs missing, duplicated or incorrectly reconstructed")
    check(sum(x["original_source_SNP_rows"] for x in all_panels)==29661,
          "original source SNP universe changed")
    check(sum(len(x["sensitivity_grid"]) for x in all_panels)==225,
          "225 scenario grid incomplete")
    check(len(set(x["gene"] for x in all_panels))==11, "11 gene identities not preserved")
    by={x["panel"]:x for x in all_panels}
    check(abs(by["P01"]["baseline"]["H4"]-0.7790043)<1e-5, "FGF5 source posterior drift")
    check(abs(by["P06"]["baseline"]["H4"]-0.7979356)<1e-5, "CALHM2 source posterior drift")
    check(all(x["baseline"]["H4"]<0.001 for x in all_panels if x["gene"]=="ALDH2"),
          "ALDH2 group-specific source result drift")
    return {
        "schema":"IS-legacy-15-assay-batch-reviewed-v1",
        "technical_status":"PASS_ALL_15_GROUPS_225_SCENARIOS",
        "scientific_status":"EXPLORATORY_HOLD",
        "actual_private_source_files_read_in_prior_local_stage":15,
        "actual_private_source_files_read_by_GitHub_runner":0,
        "original_SNP_gene_tissue_rows_QC":29661,
        "independently_assayed_genes":11,
        "gene_tissue_assays":15,
        "hypothetical_prior_N_models":225,
        "new_P3A_36_group_coloc_executed":0,
        "genuine_independent_stroke_replication":False,
        "R_coloc_package_executed":False,
        "candidate_ranking_is_not_causal_inference":True,
        "ranked_baseline_H4_review_only":sorted(
            [{"panel":x["panel"],"gene":x["gene"],"tissue":x["tissue"],
              "baseline_H4":x["baseline"]["H4"],"H3":x["baseline"]["H3"],
              "min_H4":x["min_H4"],"max_H4":x["max_H4"]}
             for x in all_panels], key=lambda x:-x["baseline_H4"]),
        "groups":all_panels
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--shards", type=int, default=1)
    ap.add_argument("--merge", nargs="*", type=Path)
    args=ap.parse_args()
    try:
        if args.merge is not None:
            result=merge(args.merge)
        else:
            check(args.source is not None,"source sufficient-statistics JSON required")
            result=reconcile(json.loads(args.source.read_text(encoding="utf-8")),
                             args.shard,args.shards)
    except (OSError,ValueError,TypeError,KeyError) as err:
        result={"schema":"IS-legacy-15-assay-batch-HOLD-v1",
                "technical_status":"HOLD_UNVERIFIED_SOURCE_OR_MODEL",
                "scientific_status":"EXPLORATORY_HOLD","error":str(err)[:250]}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k in (
        "technical_status","scientific_status","original_SNP_gene_tissue_rows_QC",
        "hypothetical_prior_N_models","error")}))
    return 0 if result["technical_status"].startswith("PASS_") else 1

if __name__ == "__main__":
    raise SystemExit(main())
