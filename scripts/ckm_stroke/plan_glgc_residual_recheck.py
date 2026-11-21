#!/usr/bin/env python3
"""Build fail-closed real-source recheck matrix for all unverified GLGC other rows."""
import argparse
import json
import os
from pathlib import Path
from scripts.ckm_stroke.inventory_glgc_zero_p_residual import inventory

SOURCE_RUNS={
  22:43431051168,21:43848671279,20:43848671279,19:43848671279,
  18:43865163287,17:43865163287,16:43865163287,
  15:43866716532,14:43866716532,13:43866716532,
  12:43867261217,11:43867261217,10:43867261217,
  9:43867962230,8:43867962230,7:43867962230,
  6:43869303834,5:43869303834,4:43869303834,
  3:43870196218,2:43870196218,1:43870196218
}
def plan(root):
    report=inventory(root)
    matrix=[]
    for row in report["rows"]:
        if row["unverified_other_candidates"]==0:
            continue
        chrom=int(row["chromosome"])
        matrix.append({
          "chrom":str(chrom),"trait":row["trait"],
          "qc_run":SOURCE_RUNS[chrom],"source_sha":row["source_commit"],
          "expected_other":row["original_other_candidates"],
          "expected_total_pzero":row["reported_source_zero_p"],
          "previously_verified":row["previously_validated_zero_p"],
          "candidate_residual":row["unverified_other_candidates"],
          "schema":row["underflow_schema"]
        })
    assert sum(x["candidate_residual"] for x in matrix)==1096
    assert sum(x["candidate_residual"] for x in matrix if x["schema"]=="LEGACY_UNREPORTED")==37
    assert len(matrix)>0 and len(matrix)<=66
    return matrix
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--evidence",required=True)
    ap.add_argument("--matrix",required=True)
    a=ap.parse_args()
    matrix=plan(a.evidence)
    Path(a.matrix).parent.mkdir(parents=True,exist_ok=True)
    Path(a.matrix).write_text(json.dumps({"include":matrix},separators=(",",":"))+"\n")
    print(json.dumps({"selected_trait_chromosome_jobs":len(matrix),
       "candidate_total":sum(x["candidate_residual"] for x in matrix),
       "legacy_candidate_total":sum(x["candidate_residual"] for x in matrix if x["schema"]=="LEGACY_UNREPORTED"),
       "reported_candidate_total":sum(x["candidate_residual"] for x in matrix if x["schema"]=="REPORTED")}))
    with open(os.environ.get("GITHUB_OUTPUT","/dev/null"),"a") as output:
        output.write("matrix="+json.dumps({"include":matrix},separators=(",",":"))+"\n")
if __name__=="__main__":
    main()
