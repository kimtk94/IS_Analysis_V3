#!/usr/bin/env python3
"""Gate bounded GLGC chromosome batches on an independently completed chr22 run.

The gate is a genuine artifact from a previous GitHub Actions run. It must
show 3/3 real-trait SNP-only technical PASS. This is NOT scientific validation.
"""
import argparse
import json
import os
from pathlib import Path

TRAITS=("LDL_C","HDL_C","TG")
SCOPE="ONE_FULL_CHROMOSOME_ONLY_NOT_WHOLE_GWAS"


def chromosomes(text, limit=3):
    if not isinstance(text,str) or not text.strip():
        raise ValueError("no chromosomes provided")
    values=[x.strip() for x in text.split(",")]
    if len(values)>limit or not values or any(not x.isdigit() for x in values):
        raise ValueError("batch must have 1-"+str(limit)+" autosome numbers")
    if any(not (1<=int(x)<=21) for x in values):
        raise ValueError("chr22 is verified separately; select only chr1-21")
    normalized=[str(int(x)) for x in values]
    if len(set(normalized))!=len(normalized):
        raise ValueError("duplicate chromosomes")
    return normalized


def verify_gate(report):
    data=json.loads(Path(report).read_text(encoding="utf-8"))
    if data.get("status")!="PASS" or data.get("chromosome")!="22":
        raise ValueError("chr22 technical aggregate PASS not established")
    if data.get("scope")!=SCOPE or data.get("scientific_status")!="HOLD":
        raise ValueError("unexpected chr22 aggregate schema/status")
    if data.get("missing_traits") or data.get("errors"):
        raise ValueError("chr22 gate reports missing or incorrect evidence")
    traits=data.get("traits",{})
    if set(traits)!=set(TRAITS):
        raise ValueError("chr22 aggregate must contain all 3 GLGC traits")
    commits=set()
    for trait in TRAITS:
        r=traits[trait]
        if r.get("technical_status")!="PASS" or not all(r.get("checks",{}).values()):
            raise ValueError("chr22 trait has incomplete technical QC: "+trait)
        if r.get("original_rows",0)<=0 or r.get("valid_snp_rows",0)<=0:
            raise ValueError("chr22 trait evidence missing real source rows")
        commit=r.get("source_commit","")
        if not commit:
            raise ValueError("source commit missing")
        commits.add(commit)
    if len(commits)!=1:
        raise ValueError("chr22 data came from different code commits")
    return next(iter(commits))


def plan(text, report, limit=3):
    chroms=chromosomes(text, limit)
    gate_commit=verify_gate(report)
    return {
        "chromosomes":chroms,
        "matrix":{"include":[{"chromosome":chrom,"trait":trait}
                             for chrom in chroms for trait in TRAITS]},
        "required_chr22_gate_commit":gate_commit,
        "status":"READY_FOR_TECHNICAL_QC_ONLY",
        "scientific_status":"NOT_RUN",
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--chromosomes",required=True)
    p.add_argument("--chr22-gate",required=True)
    p.add_argument("--output",default="artifacts/ckm_stroke/batch_plan.json")
    a=p.parse_args()
    data=plan(a.chromosomes,a.chr22_gate)
    path=Path(a.output)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(data,ensure_ascii=False))
    ghout=os.getenv("GITHUB_OUTPUT")
    if ghout:
        with open(ghout,"a",encoding="utf-8") as f:
            f.write("matrix="+json.dumps(data["matrix"],separators=(",",":"))+"\n")
            f.write("chromosomes="+json.dumps(data["chromosomes"],separators=(",",":"))+"\n")
    return 0


if __name__=="__main__":
    raise SystemExit(main())
