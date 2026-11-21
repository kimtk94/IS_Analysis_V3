#!/usr/bin/env python3
"""Fail-closed config parser for scoped GLGC follow-on batch CI pushes.

The config is not scientific evidence: every supplied run ID and expected
source count must be verified again against actual prior GitHub artifacts.
"""
import argparse
import json
import os
import re
from pathlib import Path
from scripts.ckm_stroke.plan_glgc_batch import chromosomes

SCHEMA="ckm-glgc-next-batch-approval-v2"


def parse_approval(data):
    if not isinstance(data,dict) or data.get("schema")!=SCHEMA:
        raise ValueError("unsupported GLGC approval schema")
    if data.get("scientific_status")!="HOLD":
        raise ValueError("scientific status must remain HOLD")
    run_id=data.get("prior_run_id")
    if type(run_id) is not int or run_id<=0:
        raise ValueError("positive integer prior_run_id required")
    prior_raw=data.get("prior_chromosomes")
    next_raw=data.get("next_chromosomes")
    if not isinstance(prior_raw,list) or not isinstance(next_raw,list):
        raise ValueError("chromosome lists required")
    if any(not isinstance(v,str) for v in prior_raw+next_raw):
        raise ValueError("chromosome values must be strings")
    prior=chromosomes(",".join(prior_raw))
    nxt=chromosomes(",".join(next_raw))
    if prior!=[str(int(prior[0])-k) for k in range(len(prior))]:
        raise ValueError("previous autosomes must be contiguous descending")
    if nxt!=[str(int(nxt[0])-k) for k in range(len(nxt))]:
        raise ValueError("next autosomes must be contiguous descending")
    if int(nxt[0])!=int(prior[-1])-1:
        raise ValueError("next batch does not follow prior coverage")
    sha=data.get("expected_prior_source_commit","")
    if not isinstance(sha,str) or re.fullmatch(r"[0-9a-f]{40}",sha) is None:
        raise ValueError("prior source SHA must be 40 hex characters")
    passes=data.get("prior_aggregate_statuses")
    counts=data.get("expected_prior_snp_rows")
    if not isinstance(passes,dict) or set(passes)!=set(prior):
        raise ValueError("one status per prior chromosome required")
    if any(passes[c]!="PASS" for c in prior):
        raise ValueError("prior aggregate status cannot be HOLD or missing")
    if not isinstance(counts,dict) or set(counts)!=set(prior):
        raise ValueError("one prior total per chromosome required")
    if any(type(counts[c]) is not int or counts[c]<=0 for c in prior):
        raise ValueError("invalid prior expected SNP count")
    if data.get("max_parallel_jobs")!=3 or data.get("expected_next_matrix_jobs")!=len(nxt)*3:
        raise ValueError("unsafe matrix or parallelism request")
    return {
        "schema":SCHEMA,
        "prior_run_id":str(run_id),
        "prior_chromosomes":",".join(prior),
        "next_chromosomes":",".join(nxt),
        "expected_prior_commit":sha,
        "expected_prior_snp_rows_json":json.dumps(counts,separators=(",",":"),sort_keys=True),
        "status":"READY_FOR_REAL_EVIDENCE_GATE_ONLY",
        "scientific_status":"HOLD",
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--input",required=True)
    args=parser.parse_args()
    report=parse_approval(json.loads(Path(args.input).read_text(encoding="utf-8")))
    print(json.dumps(report,ensure_ascii=False))
    if os.getenv("GITHUB_OUTPUT"):
        with open(os.environ["GITHUB_OUTPUT"],"a",encoding="utf-8") as out:
            for k in ("prior_run_id","prior_chromosomes","next_chromosomes",
                      "expected_prior_commit","expected_prior_snp_rows_json"):
                out.write(k+"="+str(report[k])+"\n")
    return 0


if __name__=="__main__":
    raise SystemExit(main())
