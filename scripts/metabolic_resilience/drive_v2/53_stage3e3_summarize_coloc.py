#!/usr/bin/env python3
from pathlib import Path
import csv
import json

ROOT=Path("/srv/is-analysis")
IN=ROOT/"results/metabolic_resilience/stage3_coloc/STAGE3E2_COLOC_ABF.tsv"
OUT=ROOT/"results/metabolic_resilience/stage3_coloc/STAGE3E3_COLOC_SUMMARY.tsv"

rows=[]
with IN.open("r",encoding="utf-8",newline="") as f:
    for r in csv.DictReader(f,delimiter="\t"):
        if r["status"]!="PASS": continue
        try:
            p12=float(r["prior_p12"])
            h4=float(r["PP.H4"])
            h3=float(r["PP.H3"])
        except Exception: continue
        if abs(p12-1e-5)>1e-12: continue
        rr=dict(r)
        if h4>=0.8:
            rr["evidence"]="STRONG_H4"
        elif h4>=0.5:
            rr["evidence"]="MODERATE_H4"
        elif h3>h4:
            rr["evidence"]="H3_GT_H4"
        else:
            rr["evidence"]="WEAK_OR_UNRESOLVED"
        rows.append(rr)

with OUT.open("w",encoding="utf-8",newline="") as f:
    fields=list(rows[0].keys()) if rows else ["status"]
    wr=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n")
    wr.writeheader();wr.writerows(rows)

print(json.dumps({
    "pairs":len(rows),
    "strong_h4":sum(r["evidence"]=="STRONG_H4" for r in rows),
    "moderate_h4":sum(r["evidence"]=="MODERATE_H4" for r in rows),
    "h3_gt_h4":sum(r["evidence"]=="H3_GT_H4" for r in rows),
},indent=2))
