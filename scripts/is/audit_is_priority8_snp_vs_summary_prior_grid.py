#!/usr/bin/env python3
"""Reconcile 8 direct SNP-level coloc prior grids vs 646-row summary reweight.

These two computations need not be mathematically identical at changed priors.
Fails if baseline replication is not exact, or prior-adjust deviance is large.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

GENE_ID={"FGF5":"ENSG00000138675",
         "SH3PXD2A":"ENSG00000107957",
         "COL4A2":"ENSG00000134871",
         "CALHM2":"ENSG00000138172",
         "NEURL1":"ENSG00000107954",
         "C4orf22":"ENSG00000197826",
         "INA":"ENSG00000148798",
         "COL4A1":"ENSG00000187498"}


def read(file):
    with file.open(encoding="utf8",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))


def run(priors,replay,out):
    if out.exists():
        raise FileExistsError("REFUSE_OVERWRITE")
    orig=read(priors)
    newly=read(replay)
    if len(orig)!=3230 or len(newly)!=120:
        raise ValueError("UNEXPECTED_GRID_DIMENSIONS")
    keys={}
    for v in orig:
        k=(v["locus"],v["dataset_key"],v["gene_base"],float(v["conditional_p12"]))
        if k in keys:raise ValueError("DUPLICATE_ORIGINAL_POSTERIOR_GRID")
        keys[k]=v
    vals=[]
    for r in newly:
        if r["qtl_n_mode"]!="FIRST":
            continue
        k=(r["locus"],r["tissue"],GENE_ID[r["gene"]],float(r["p12"]))
        old=keys.get(k)
        if old is None:raise ValueError("SNP_REPLAY_MISSING_PRIOR_MATCH")
        diffs=[float(old[f"PP_H{i}"])-float(r[f"PP_H{i}"]) for i in range(5)]
        if abs(float(r["p12"])-1e-5)<1e-12 and max(map(abs,diffs))>1e-8:
            raise ValueError("BASELINE_SOURCE_POSTERIOR_DISAGREEMENT")
        vals.append(dict(gene=r["gene"],locus=r["locus"],tissue=r["tissue"],
            p12=float(r["p12"]),direct_SNP_H4=float(r["PP_H4"]),
            posterior_reweighted_H4=float(old["PP_H4"]),
            delta_reweighted_minus_direct_H4=diffs[4],
            max_abs_H0_to_H4_delta=max(abs(x) for x in diffs)))
    if len(vals)!=40:
        raise ValueError("EXPECTED_3_GENES_5_PRIORS")
    largest=max(x["max_abs_H0_to_H4_delta"] for x in vals)
    if largest > 0.001:
        raise ValueError("PRIOR_REWEIGHT_DISCREPANCY_UNEXPECTEDLY_LARGE")
    out.parent.mkdir(parents=True,exist_ok=True)
    result=dict(
        schema="IS_DIRECT_SNP_VS_SUMMARY_PRIOR_REWEIGHT_V1",
        status="BASELINE_EXACT_CHANGED_PRIOR_APPROXIMATION",
        pairs=8,prior_grid=5,paired_analyses=len(vals),
        baseline_all_five_H0_to_H4_reproduced=True,
        max_abs_posterior_delta=largest,
        reweight_method="coloc_5.2.3_prior.adjust_historical_ABF_H0_H4",
        direct_method="coloc_5.2.3_coloc.abf_original_per_SNP_beta_SE_MAF",
        policy="Direct per-SNP coloc.abf takes precedence for retained genes; summary reweighting remains descriptive",
        input_hashes={"summary_prior_grid_sha256":hashlib.sha256(priors.read_bytes()).hexdigest(),
                      "direct_snp_grid_sha256":hashlib.sha256(replay.read_bytes()).hexdigest()},
        comparisons=vals
    )
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print(f"SNP_P12_REPLAY_VS_REWEIGHT_PASS comparisons={len(vals)} maximum_difference={largest:.9g}")
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--priors",type=Path,required=True)
    p.add_argument("--replay",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    run(a.priors,a.replay,a.out)


if __name__=="__main__":
    main()
