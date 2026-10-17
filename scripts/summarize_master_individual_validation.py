#!/usr/bin/env python3
"""Summarize one or more individual-level validation result tables.

Input manifest TSV:
gene_symbol, predictor, result_file

The underlying participant-level data never enter this script. Only aggregate
model summaries are combined, preserving the KoGES restricted-data boundary.
"""
from __future__ import annotations
import argparse,csv
from pathlib import Path

def read_tsv(p):
    with Path(p).open("r",encoding="utf-8",newline="") as f:
        yield from csv.DictReader(f,delimiter="\t")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--manifest",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    rows=[]
    for m in read_tsv(a.manifest):
        gene=m.get("gene_symbol","").strip()
        predictor=m.get("predictor","").strip()
        rf=Path(m["result_file"])
        for r in read_tsv(rf):
            x={"gene_symbol":gene,"predictor":predictor}
            x.update(r)
            p=r.get("p","")
            try: pv=float(p)
            except Exception: pv=None
            x["nominal_p05"]=int(pv is not None and pv<0.05)
            rows.append(x)
    fields=["gene_symbol","predictor","model","endpoint","method","status","n","events",
            "beta","se","statistic","p","nominal_p05","effect_scale","effect_ratio",
            "ci95_low","ci95_high","ratio_ci95_low","ratio_ci95_high"]
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore")
        w.writeheader();w.writerows(rows)
    print(f"PASS rows={len(rows)} output={a.output}")

if __name__=="__main__":main()
