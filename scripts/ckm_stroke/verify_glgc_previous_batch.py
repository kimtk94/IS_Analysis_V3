#!/usr/bin/env python3
"""Verify all three completed GLGC chromosome summaries before next batch.

Designed for CI artifacts from one prior bounded-batch Actions run.
This checks technical evidence only; it NEVER declares statistical or
biological readiness.
"""
import argparse
import json
import re
from pathlib import Path
from scripts.ckm_stroke.plan_glgc_batch import chromosomes, TRAITS, SCOPE

REQUIRED_CHECKS={
    "observed_header","shard_complete","row_accounting","variant_qc_rows",
    "no_duplicate_variant_keys","no_p_diagnostic_mismatches",
    "source_identity","technical_gate_passed","source_region_only"
}


def verify_previous(folder, prior_text, next_text):
    prior=chromosomes(prior_text)
    nxt=chromosomes(next_text)
    if prior!=[str(i) for i in range(int(prior[0]),int(prior[0])-len(prior),-1)]:
        raise ValueError("prior must be descending contiguous autosomes")
    if nxt!=[str(i) for i in range(int(nxt[0]),int(nxt[0])-len(nxt),-1)]:
        raise ValueError("next must be descending contiguous autosomes")
    if int(nxt[0])!=int(prior[-1])-1:
        raise ValueError("next batch must immediately follow previous batch")
    folder=Path(folder)
    names=sorted(folder.glob("glgc_chr*-aggregate.json"))
    expected={f"glgc_chr{c}-aggregate.json" for c in prior}
    if {p.name for p in names}!=expected:
        raise ValueError("previous run is missing aggregate(s) or includes unexpected chromosome evidence")
    commits=set()
    counts={}
    for chrom in prior:
        report=json.loads((folder/f"glgc_chr{chrom}-aggregate.json").read_text(encoding="utf-8"))
        if report.get("scope")!=SCOPE or report.get("chromosome")!=chrom:
            raise ValueError("invalid previous aggregate scope/chromosome "+chrom)
        if report.get("status")!="PASS" or report.get("scientific_status")!="HOLD":
            raise ValueError("previous chromosome not technically verified "+chrom)
        if report.get("missing_traits") or report.get("errors") or set(report.get("traits",{}))!=set(TRAITS):
            raise ValueError("previous chromosome trait set incomplete "+chrom)
        commit=report.get("source_commit","")
        if not re.fullmatch("[0-9a-f]{40}",commit):
            raise ValueError("aggregate source commit missing or malformed "+chrom)
        commits.add(commit)
        total=0
        for trait in TRAITS:
            x=report["traits"][trait]
            if x.get("source_commit")!=commit or x.get("technical_status")!="PASS":
                raise ValueError("trait provenance or technical status mismatch "+chrom+" "+trait)
            checks=x.get("checks",{})
            if not REQUIRED_CHECKS.issubset(set(checks)) or any(checks.get(k) is not True for k in REQUIRED_CHECKS):
                raise ValueError("missing/failed mandatory technical check "+chrom+" "+trait)
            if x.get("original_rows",0)<=0 or x.get("valid_snp_rows",0)<=0 or x.get("duplicate_keys")!=0:
                raise ValueError("invalid SNP evidence "+chrom+" "+trait)
            total+=x["valid_snp_rows"]
        if total!=report.get("snp_rows_total"):
            raise ValueError("aggregate SNP accounting mismatch "+chrom)
        counts[chrom]=total
    if len(commits)!=1:
        raise ValueError("previous batches mixed code revisions")
    return {
        "status":"PASS","scientific_status":"HOLD",
        "prior_chromosomes":prior,"next_chromosomes":nxt,
        "prior_commit":next(iter(commits)),
        "prior_chr_snp_totals":counts,
        "matrix":{"include":[{"chromosome":c,"trait":t} for c in nxt for t in TRAITS]},
        "scope":"TECHNICAL_QC_ONLY",
        "note":"Validated prior per-chromosome evidence; whole-genome scientific status remains HOLD."
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--prior-dir",required=True)
    parser.add_argument("--prior-chromosomes",required=True)
    parser.add_argument("--next-chromosomes",required=True)
    parser.add_argument("--output",required=True)
    args=parser.parse_args()
    report=verify_previous(args.prior_dir,args.prior_chromosomes,args.next_chromosomes)
    output=Path(args.output)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report))
    if "GITHUB_OUTPUT" in __import__("os").environ:
        with open(__import__("os").environ["GITHUB_OUTPUT"],"a",encoding="utf-8") as out:
            out.write("matrix="+json.dumps(report["matrix"],separators=(",",":"))+"\n")
            out.write("chromosomes="+json.dumps(report["next_chromosomes"],separators=(",",":"))+"\n")
    return 0


if __name__=="__main__":
    raise SystemExit(main())
