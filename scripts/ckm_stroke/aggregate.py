#!/usr/bin/env python3
"""Summarize seven evidence artifacts without confusing CI success and science."""
import argparse
import datetime
import json
import os
from pathlib import Path

EXPECTED = tuple("p"+str(i) for i in range(7))

def collect(folder):
    records={}
    for path in folder.rglob("*.json"):
        data=json.loads(path.read_text(encoding="utf-8"))
        stage=data.get("stage")
        if stage in EXPECTED:
            if stage in records: raise ValueError("duplicate stage evidence: "+stage)
            records[stage]=data
    missing=[x for x in EXPECTED if x not in records]
    technical={x:records[x]["status"] if x in records else "MISSING" for x in EXPECTED}
    scientific={x:records[x]["scientific_status"] if x in records else "NOT_RUN" for x in EXPECTED}
    fully_validated=sum(v=="PASS" for v in scientific.values())
    return {"as_of_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "expected_stages":len(EXPECTED),"observed_stages":len(records),
            "technical":technical,"scientific":scientific,"missing":missing,
            "fully_validated_scientific_stages":fully_validated,
            "overall_scientific_status":"PASS" if fully_validated==len(EXPECTED) else "HOLD"}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--folder",default="artifacts/ckm_stroke")
    p.add_argument("--out",default="artifacts/ckm_stroke_summary.json")
    a=p.parse_args()
    summary=collect(Path(a.folder))
    dest=Path(a.out)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    lines=["## CKM→Stroke staged evidence status","","| Stage | Technical | Scientific |",
           "|---|---|---|"]
    for stage in EXPECTED:
        lines.append("| "+stage.upper()+" | "+summary["technical"][stage]+
                     " | "+summary["scientific"][stage]+" |")
    lines.extend(["",
                  "**Scientifically validated stages: "+str(summary["fully_validated_scientific_stages"])+"/7**",
                  "",
                  "GitHub Actions SUCCESS means the evidence gates executed, not that GWAS/MR/coloc research completed.",
                  ""])
    output="\n".join(lines)
    print(output)
    if os.getenv("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"],"a",encoding="utf-8") as f:
            f.write(output)
    if summary["missing"] or "FAIL" in summary["technical"].values():
        return 1
    return 0


if __name__=="__main__":
    raise SystemExit(main())
