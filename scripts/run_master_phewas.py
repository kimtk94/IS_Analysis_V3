#!/usr/bin/env python3
"""Phenome-wide association evidence standardizer.

This stage ingests one or more candidate-level PheWAS result tables and applies
multiple-testing correction within the supplied result universe.

Manifest TSV:
source_file    source_name    phenotype_group

Recognized source columns:
gene_symbol/protein_id, phenotype/trait, beta/effect, se, p/p_value,
optional effect_scale, direction_is_adverse, category.

The module does not infer clinical safety from association alone. Instead it
separates statistically significant phenome associations from user-supplied
adverse-direction annotations.
"""
from __future__ import annotations
import argparse,csv,gzip,math
from pathlib import Path
from collections import defaultdict

def fnum(x):
    if x is None:return None
    try:v=float(str(x).strip())
    except Exception:return None
    return v if math.isfinite(v) else None

def read_tsv(path):
    path=Path(path);op=gzip.open if path.suffix==".gz" else open
    with op(path,"rt",encoding="utf-8",newline="") as f:
        yield from csv.DictReader(f,delimiter="\t")

def write_tsv(path,rows,fields):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore")
        w.writeheader();w.writerows(rows)

def first(r,names):
    for n in names:
        if n in r and str(r[n]).strip()!="":return r[n]
    return None

def boolish(x):
    s=str(x or "").strip().lower()
    if s in {"1","true","yes","y","adverse","harmful"}:return 1
    if s in {"0","false","no","n","beneficial","protective"}:return 0
    return ""

def bh(ps):
    valid=sorted((p,i) for i,p in enumerate(ps) if p is not None and math.isfinite(p))
    out=[None]*len(ps);m=len(valid);prev=1.0
    for rank in range(m,0,-1):
        p,i=valid[rank-1]
        q=min(prev,p*m/rank,1.0);out[i]=q;prev=q
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--manifest",type=Path,required=True)
    ap.add_argument("--long-output",type=Path,required=True)
    ap.add_argument("--summary-output",type=Path,required=True)
    ap.add_argument("--fdr-threshold",type=float,default=0.05)
    args=ap.parse_args()

    rows=[]
    for m in read_tsv(args.manifest):
        src=Path(m["source_file"]);sname=str(m.get("source_name","")).strip()
        group=str(m.get("phenotype_group","")).strip()
        for r in read_tsv(src):
            gene=str(first(r,["gene_symbol","gene","Gene"]) or "").strip().upper()
            protein=str(first(r,["protein_id","protein","Protein"]) or "").strip()
            entity=gene or protein
            trait=str(first(r,["phenotype","trait","outcome"]) or "").strip()
            if not entity or not trait:continue
            beta=fnum(first(r,["beta","effect","estimate"]))
            se=fnum(first(r,["se","stderr","standard_error"]))
            p=fnum(first(r,["p","p_value","P","pval"]))
            adverse=boolish(first(r,["direction_is_adverse","adverse","safety_flag"]))
            category=str(first(r,["category","phenotype_category"]) or group).strip()
            rows.append({
              "entity":entity,"gene_symbol":gene,"protein_id":protein,
              "source_name":sname,"phenotype_group":group,"category":category,
              "phenotype":trait,"beta":beta if beta is not None else "",
              "se":se if se is not None else "","p":p if p is not None else "",
              "effect_scale":str(r.get("effect_scale","")).strip(),
              "direction_is_adverse":adverse,
            })

    qs=bh([fnum(r["p"]) for r in rows])
    for r,q in zip(rows,qs):
        r["fdr_bh"]=q if q is not None else ""
        sig=int(q is not None and q<args.fdr_threshold)
        r["fdr_significant"]=sig
        adverse=r["direction_is_adverse"]
        r["significant_adverse_flag"]=int(sig and adverse==1)
        r["significant_nonadverse_flag"]=int(sig and adverse==0)

    long_fields=["entity","gene_symbol","protein_id","source_name","phenotype_group","category","phenotype",
                 "beta","se","p","fdr_bh","fdr_significant","effect_scale","direction_is_adverse",
                 "significant_adverse_flag","significant_nonadverse_flag"]
    write_tsv(args.long_output,rows,long_fields)

    by=defaultdict(list)
    for r in rows:by[r["entity"]].append(r)
    summary=[]
    for ent,xs in sorted(by.items()):
        sig=[x for x in xs if x["fdr_significant"]==1]
        adv=[x for x in sig if x["significant_adverse_flag"]==1]
        nonadv=[x for x in sig if x["significant_nonadverse_flag"]==1]
        top=min(xs,key=lambda x:fnum(x["p"]) if fnum(x["p"]) is not None else 1.0)
        summary.append({
          "entity":ent,
          "gene_symbol":next((x["gene_symbol"] for x in xs if x["gene_symbol"]),""),
          "protein_id":next((x["protein_id"] for x in xs if x["protein_id"]),""),
          "phewas_traits_tested":len(xs),
          "phewas_fdr_significant_n":len(sig),
          "phewas_significant_adverse_n":len(adv),
          "phewas_significant_nonadverse_n":len(nonadv),
          "top_phenotype":top["phenotype"],
          "top_p":top["p"],
          "top_fdr":top["fdr_bh"],
          "safety_attention":int(len(adv)>0),
        })
    sfields=["entity","gene_symbol","protein_id","phewas_traits_tested","phewas_fdr_significant_n",
             "phewas_significant_adverse_n","phewas_significant_nonadverse_n",
             "top_phenotype","top_p","top_fdr","safety_attention"]
    write_tsv(args.summary_output,summary,sfields)
    print(f"PASS entities={len(summary)} associations={len(rows)} output={args.summary_output}")

if __name__=="__main__":main()
