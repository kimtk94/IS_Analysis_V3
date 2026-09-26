#!/usr/bin/env python3
from pathlib import Path
import csv
import math
import json

ROOT = Path("/srv/is-analysis")
OUT = ROOT / "results/metabolic_resilience/stage3_confirmatory_mr"
MR = OUT / "STAGE3D3_MR_RESULTS.tsv"

FAVORABLE = {
    "HDL": 1,
    "TG": -1,
    "SBP": -1,
    "DBP": -1,
    "BMI": -1,
    "WHR": -1,
    "FG": -1,
    "HBA1C": -1,
    "T2D": -1,
}

PRIMARY_METHOD_ORDER = ["WALD","IVW_MRE","IVW_FIXED"]

def bh(pairs):
    # pairs: [(index,p)]
    pairs = [(i,p) for i,p in pairs if p is not None and math.isfinite(p)]
    pairs.sort(key=lambda x:x[1])
    n = len(pairs)
    q = {}
    prev = 1.0
    for rank in range(n,0,-1):
        i,p = pairs[rank-1]
        val = min(prev, p*n/rank)
        q[i] = val
        prev = val
    return q

rows = []
with MR.open("r", encoding="utf-8", newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

for r in rows:
    for k in ["beta","se","p","Q","Q_p","egger_intercept_p"]:
        try:
            r[k] = float(r[k]) if r[k] not in ("","NA") else None
        except Exception:
            r[k] = None

# Primary result per gene/mode/trait.
groups = {}
for r in rows:
    key=(r["gene"],r["mode"],r["trait"],r["domain"])
    groups.setdefault(key,[]).append(r)

primary=[]
for key, arr in groups.items():
    chosen=None
    for method in PRIMARY_METHOD_ORDER:
        hit=[r for r in arr if r["method"]==method]
        if hit:
            chosen=hit[0]; break
    if chosen:
        rr=dict(chosen)
        s=FAVORABLE.get(rr["trait"])
        rr["direction_favorable"] = (
            int(rr["beta"]*s > 0) if s is not None and rr["beta"] is not None else ""
        )
        primary.append(rr)

# FDR separately by mode across candidate x trait primary tests.
for mode in sorted(set(r["mode"] for r in primary)):
    idx=[i for i,r in enumerate(primary) if r["mode"]==mode and r["p"] is not None]
    q=bh([(i,primary[i]["p"]) for i in idx])
    for i in idx:
        primary[i]["FDR_mode"] = q.get(i)
for r in primary:
    r.setdefault("FDR_mode",None)

out_primary = OUT / "STAGE3D5_PRIMARY_MR_SUMMARY.tsv"
fields = [
    "gene","mode","trait","domain","k","method","beta","se","p","FDR_mode",
    "direction_favorable","Q","Q_df","Q_p","egger_intercept_p","notes"
]
with out_primary.open("w", encoding="utf-8", newline="") as f:
    wr=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore")
    wr.writeheader(); wr.writerows(primary)

cand=[]
for gene in sorted(set(r["gene"] for r in primary)):
    for mode in sorted(set(r["mode"] for r in primary)):
        x=[r for r in primary if r["gene"]==gene and r["mode"]==mode]
        if not x: continue
        metabolic=[r for r in x if r["domain"]!="disease_validation"]
        t2d=[r for r in x if r["trait"]=="T2D"]
        cand.append({
            "gene":gene,
            "mode":mode,
            "traits_available":len(x),
            "metabolic_traits_available":len(metabolic),
            "nominal_p_lt_0p05":sum(r["p"] is not None and r["p"]<0.05 for r in metabolic),
            "fdr_lt_0p05":sum(r["FDR_mode"] is not None and r["FDR_mode"]<0.05 for r in metabolic),
            "favorable_direction":sum(r["direction_favorable"]==1 for r in metabolic),
            "t2d_available":len(t2d),
            "t2d_beta":t2d[0]["beta"] if t2d else "",
            "t2d_p":t2d[0]["p"] if t2d else "",
            "t2d_favorable":t2d[0]["direction_favorable"] if t2d else "",
        })

out_cand=OUT/"STAGE3D5_CANDIDATE_CONFIRMATORY_SUMMARY.tsv"
with out_cand.open("w",encoding="utf-8",newline="") as f:
    fields2=list(cand[0].keys()) if cand else ["gene"]
    wr=csv.DictWriter(f,fieldnames=fields2,delimiter="\t",lineterminator="\n")
    wr.writeheader();wr.writerows(cand)

print(json.dumps({
    "primary_rows":len(primary),
    "candidate_mode_rows":len(cand),
    "primary_file":str(out_primary),
    "candidate_file":str(out_cand),
},indent=2))
