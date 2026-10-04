#!/usr/bin/env python3
"""Summarize protein -> risk-factor -> disease mechanism evidence.

This is a mechanistic prioritization layer, not formal mediation estimation.
It combines disease MR with risk-factor MR results and reports whether the
protein affects established disease risk factors in a direction compatible
with the observed disease association.
"""
from __future__ import annotations
import argparse,csv,gzip,math
from pathlib import Path
from collections import defaultdict

def fnum(x):
    try:v=float(str(x).strip())
    except Exception:return None
    return v if math.isfinite(v) else None

def read_tsv(path):
    path=Path(path);op=gzip.open if path.suffix==".gz" else open
    with op(path,"rt",encoding="utf-8",newline="") as f:yield from csv.DictReader(f,delimiter="\t")

def first(r,names):
    for n in names:
        if n in r and str(r[n]).strip()!="":return r[n]
    return None

def load_mr(path):
    out={}
    for r in read_tsv(path):
        gene=str(r.get("gene_symbol","")).strip();protein=str(r.get("protein_id","")).strip()
        ent=gene or protein;pheno=str(r.get("phenotype","")).strip()
        if not ent or not pheno:continue
        out[(ent,pheno)]={
          "gene_symbol":gene,"protein_id":protein,
          "beta":fnum(first(r,["beta","wald_beta","anchor_wald_beta","ivw_beta"])),
          "p":fnum(first(r,["p","wald_p","anchor_wald_p","ivw_p"]))
        }
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--disease-mr",type=Path,required=True)
    ap.add_argument("--risk-mr",type=Path,required=True)
    ap.add_argument("--disease-phenotype",required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--risk-sign-map",type=Path,help="TSV phenotype,disease_risk_sign (+1 or -1)")
    ap.add_argument("--p-threshold",type=float,default=0.05)
    args=ap.parse_args()

    d=load_mr(args.disease_mr); r=load_mr(args.risk_mr)
    signs={}
    if args.risk_sign_map:
        for x in read_tsv(args.risk_sign_map):
            signs[str(x["phenotype"]).strip()]=float(x["disease_risk_sign"])

    entities=sorted({e for e,p in d if p==args.disease_phenotype})
    rows=[]
    riskphenos=sorted({p for e,p in r})
    for ent in entities:
        ds=d.get((ent,args.disease_phenotype))
        db=ds["beta"] if ds else None
        for ph in riskphenos:
            rr=r.get((ent,ph))
            if not rr:continue
            rb=rr["beta"];rp=rr["p"];risk_sign=signs.get(ph)
            compatible=""
            if db not in (None,0) and rb not in (None,0) and risk_sign in (-1,1):
                # risk_sign=+1 means higher risk factor raises disease risk.
                compatible=int((rb*risk_sign)*db>0)
            rows.append({
              "entity":ent,
              "gene_symbol":ds["gene_symbol"] if ds else rr["gene_symbol"],
              "protein_id":ds["protein_id"] if ds else rr["protein_id"],
              "disease_phenotype":args.disease_phenotype,
              "disease_beta":db if db is not None else "",
              "disease_p":ds["p"] if ds else "",
              "risk_factor":ph,
              "risk_factor_beta":rb if rb is not None else "",
              "risk_factor_p":rp if rp is not None else "",
              "risk_factor_disease_sign":risk_sign if risk_sign is not None else "",
              "risk_factor_p05":int(rp is not None and rp<args.p_threshold),
              "direction_compatible":compatible,
              "mechanism_support":int(rp is not None and rp<args.p_threshold and compatible==1),
            })

    args.output.parent.mkdir(parents=True,exist_ok=True)
    fields=["entity","gene_symbol","protein_id","disease_phenotype","disease_beta","disease_p","risk_factor","risk_factor_beta","risk_factor_p","risk_factor_disease_sign","risk_factor_p05","direction_compatible","mechanism_support"]
    with args.output.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n")
        w.writeheader();w.writerows(rows)
    print(f"PASS rows={len(rows)} output={args.output}")

if __name__=="__main__":main()
