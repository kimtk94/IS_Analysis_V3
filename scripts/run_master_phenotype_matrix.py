#!/usr/bin/env python3
"""Build a standardized phenotype/subtype evidence matrix from MASTER MR outputs.

Inputs are one or more MR summary tables. Each row must contain gene_symbol or
protein_id, phenotype, and an effect estimate compatible with MASTER aliases.
Optional coloc summaries can be supplied through a manifest.

Manifest schema (TSV):
mr_file    label    group    coloc_file
where group is optional (e.g. kidney_function, stroke_subtype).
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
    path=Path(path); op=gzip.open if path.suffix==".gz" else open
    with op(path,"rt",encoding="utf-8",newline="") as f:
        yield from csv.DictReader(f,delimiter="\t")

def write_tsv(path,rows,fields):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore")
        w.writeheader();w.writerows(rows)

def first(r,names):
    for n in names:
        if n in r and str(r[n]).strip()!="":return r[n]
    return None

def load_coloc(path):
    if not path:return {}
    out={}
    for r in read_tsv(path):
        locus=str(first(r,["locus","gene_symbol","protein_id"]) or "").strip()
        if not locus:continue
        h4=fnum(first(r,["PP.H4","PP.H4.abf","max_PP.H4"]))
        p12=fnum(r.get("p12"))
        score=0 if p12 is None else abs(math.log10(p12)-math.log10(1e-5))
        if locus not in out or score<out[locus][1]:
            out[locus]=(h4,score)
    return {k:v[0] for k,v in out.items()}

def load_mr(path,label,group,coloc):
    rows=[]
    for r in read_tsv(path):
        gene=str(r.get("gene_symbol","")).strip()
        protein=str(r.get("protein_id","")).strip()
        key=gene or protein
        if not key:continue
        pheno=str(r.get("phenotype","")).strip() or label
        beta=fnum(first(r,["beta","wald_beta","anchor_wald_beta","ivw_beta"]))
        se=fnum(first(r,["se","wald_se","anchor_wald_se","ivw_se"]))
        p=fnum(first(r,["p","wald_p","anchor_wald_p","ivw_p"]))
        fdr=fnum(first(r,["fdr_bh","wald_fdr_bh"]))
        rows.append({
          "entity":key,"gene_symbol":gene,"protein_id":protein,
          "phenotype":pheno,"label":label,"group":group,
          "beta":beta,"se":se,"p":p,"fdr":fdr,
          "PP.H4":coloc.get(gene or protein)
        })
    return rows

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--manifest",type=Path,required=True)
    ap.add_argument("--long-output",type=Path,required=True)
    ap.add_argument("--wide-output",type=Path,required=True)
    ap.add_argument("--p-threshold",type=float,default=0.05)
    ap.add_argument("--coloc-h4-threshold",type=float,default=0.80)
    args=ap.parse_args()

    allrows=[]
    for m in read_tsv(args.manifest):
        mr=Path(m["mr_file"])
        coloc=load_coloc(Path(m["coloc_file"])) if str(m.get("coloc_file","")).strip() else {}
        allrows.extend(load_mr(mr,str(m.get("label","")).strip(),str(m.get("group","")).strip(),coloc))

    for r in allrows:
        r["direction"]="positive" if r["beta"] is not None and r["beta"]>0 else "negative" if r["beta"] is not None and r["beta"]<0 else ""
        r["p05"]=int(r["p"] is not None and r["p"]<args.p_threshold)
        r["coloc_strong"]=int(r["PP.H4"] is not None and r["PP.H4"]>=args.coloc_h4_threshold)
        if r["p05"] and r["coloc_strong"]:r["evidence"]="MR+COLOC"
        elif r["p05"]:r["evidence"]="MR"
        elif r["coloc_strong"]:r["evidence"]="COLOC_ONLY"
        else:r["evidence"]="NONE"

    fields=["entity","gene_symbol","protein_id","phenotype","label","group","beta","se","p","fdr","PP.H4","direction","p05","coloc_strong","evidence"]
    write_tsv(args.long_output,allrows,fields)

    phenos=sorted({r["label"] or r["phenotype"] for r in allrows})
    by=defaultdict(dict); meta={}
    for r in allrows:
        ent=r["entity"]; lab=r["label"] or r["phenotype"]
        by[ent][lab]=r
        meta[ent]=(r["gene_symbol"],r["protein_id"])
    wide=[]
    for ent in sorted(by):
        row={"entity":ent,"gene_symbol":meta[ent][0],"protein_id":meta[ent][1]}
        signs=[]
        for lab in phenos:
            r=by[ent].get(lab)
            row[f"{lab}__beta"]=r["beta"] if r else ""
            row[f"{lab}__p"]=r["p"] if r else ""
            row[f"{lab}__PP.H4"]=r["PP.H4"] if r and r["PP.H4"] is not None else ""
            row[f"{lab}__evidence"]=r["evidence"] if r else ""
            if r and r["beta"] not in (None,0):signs.append(1 if r["beta"]>0 else -1)
        row["direction_consistent"]=int(bool(signs) and len(set(signs))==1)
        row["n_phenotypes"]=len(by[ent])
        row["n_mr_p05"]=sum(x["p05"] for x in by[ent].values())
        row["n_coloc_strong"]=sum(x["coloc_strong"] for x in by[ent].values())
        wide.append(row)
    wide_fields=["entity","gene_symbol","protein_id","direction_consistent","n_phenotypes","n_mr_p05","n_coloc_strong"]
    for lab in phenos:
        wide_fields += [f"{lab}__beta",f"{lab}__p",f"{lab}__PP.H4",f"{lab}__evidence"]
    write_tsv(args.wide_output,wide,wide_fields)
    print(f"PASS entities={len(wide)} phenotypes={len(phenos)} long={args.long_output} wide={args.wide_output}")

if __name__=="__main__":main()
