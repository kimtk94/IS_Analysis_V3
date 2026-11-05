#!/usr/bin/env python3
"""Integrate tissue, single-cell, and spatial localization evidence.

This module standardizes heterogeneous annotation tables into one candidate-level
evidence matrix. It does not infer causality from expression/localization.

Input manifests are TSV files with:
  source_file    source_name    tissue_or_region

Supported source-file columns are recognized through aliases.

Tissue-level source rows should contain:
  gene_symbol + expression/value
Optional:
  tissue, compartment, specificity/tau, evidence_type

Single-cell source rows should contain:
  gene_symbol + cell_type + expression/value
Optional:
  compartment, specificity/tau

Spatial source rows should contain:
  gene_symbol + region/compartment + expression/value
Optional:
  spot_fraction, enrichment, evidence_type
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
    path=Path(path)
    op=gzip.open if path.suffix==".gz" else open
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

def compartment_from_cell(cell):
    x=(cell or "").lower()
    groups={
      "immune":["macrophage","monocyte","neutrophil","t-cell","t cell","b-cell","b cell","cdc","dendritic","mast","nk cell","plasma cell","microglia"],
      "vascular":["endothelial","vascular smooth muscle","pericyte"],
      "renal_epithelial":["proximal tubule","distal tubule","loop of henle","collecting duct","connecting tubule","papillary tip","podocyte","urothelial","epithelial"],
      "neural":["neuron","astrocyte","oligodendro","opc","ependymal"],
      "stromal":["fibroblast","mesangial","stromal"],
    }
    for g,terms in groups.items():
        if any(t in x for t in terms):return g
    return "other"

def load_manifest(path):
    return list(read_tsv(path)) if path else []

def load_tissue(manifest):
    out=[]
    for m in load_manifest(manifest):
        source=Path(m["source_file"]); sname=str(m.get("source_name","")).strip()
        default_region=str(m.get("tissue_or_region","")).strip()
        for r in read_tsv(source):
            gene=str(first(r,["gene_symbol","gene","Gene"]) or "").strip().upper()
            if not gene:continue
            value=fnum(first(r,["expression","value","nTPM","TPM","nCPM","mean_expression"]))
            tissue=str(first(r,["tissue","region","organ"]) or default_region).strip()
            tau=fnum(first(r,["tau_specificity","tau","specificity"]))
            out.append({
              "gene_symbol":gene,"source_name":sname,"tissue":tissue,
              "expression":value,"tau_specificity":tau,
              "evidence_type":str(r.get("evidence_type","")).strip()
            })
    return out

def load_cell(manifest):
    out=[]
    for m in load_manifest(manifest):
        source=Path(m["source_file"]);sname=str(m.get("source_name","")).strip()
        default_region=str(m.get("tissue_or_region","")).strip()
        for r in read_tsv(source):
            gene=str(first(r,["gene_symbol","gene","Gene"]) or "").strip().upper()
            cell=str(first(r,["cell_type","celltype","cluster","cell"]) or "").strip()
            if not gene or not cell:continue
            value=fnum(first(r,["expression","value","nCPM","mean_expression","avg_expr","avg_log2FC"]))
            tau=fnum(first(r,["tau_specificity","tau","specificity"]))
            comp=str(first(r,["compartment"]) or compartment_from_cell(cell)).strip()
            tissue=str(first(r,["tissue","region","organ"]) or default_region).strip()
            out.append({
              "gene_symbol":gene,"source_name":sname,"tissue":tissue,
              "cell_type":cell,"compartment":comp,"expression":value,
              "tau_specificity":tau
            })
    return out

def load_spatial(manifest):
    out=[]
    for m in load_manifest(manifest):
        source=Path(m["source_file"]);sname=str(m.get("source_name","")).strip()
        default_region=str(m.get("tissue_or_region","")).strip()
        for r in read_tsv(source):
            gene=str(first(r,["gene_symbol","gene","Gene"]) or "").strip().upper()
            region=str(first(r,["spatial_region","region","compartment","anatomical_region"]) or default_region).strip()
            if not gene or not region:continue
            value=fnum(first(r,["expression","value","mean_expression","nCPM","TPM"]))
            enrich=fnum(first(r,["enrichment","fold_enrichment","log2fc","avg_log2FC"]))
            frac=fnum(first(r,["spot_fraction","fraction","pct_expressing","pct.1"]))
            out.append({
              "gene_symbol":gene,"source_name":sname,"spatial_region":region,
              "expression":value,"enrichment":enrich,"spot_fraction":frac,
              "evidence_type":str(r.get("evidence_type","")).strip()
            })
    return out

def summarize_tissue(rows):
    by=defaultdict(list)
    for r in rows:by[r["gene_symbol"]].append(r)
    out={}
    for gene,xs in by.items():
        ranked=sorted(xs,key=lambda r:(r["expression"] is not None, r["expression"] if r["expression"] is not None else -math.inf),reverse=True)
        top=ranked[0]
        values=[r["expression"] for r in xs if r["expression"] is not None and r["expression"]>=0]
        total=sum(values)
        share=(top["expression"]/total) if total>0 and top["expression"] is not None else None
        max_tau=max([r["tau_specificity"] for r in xs if r["tau_specificity"] is not None],default=None)
        out[gene]={
          "tissue_source_n":len({r["source_name"] for r in xs}),
          "top_tissue":top["tissue"],
          "top_tissue_expression":top["expression"],
          "top_tissue_share":share,
          "max_tissue_specificity":max_tau,
        }
    return out

def summarize_cell(rows):
    by=defaultdict(list)
    for r in rows:by[r["gene_symbol"]].append(r)
    out={}
    for gene,xs in by.items():
        ranked=sorted(xs,key=lambda r:(r["expression"] is not None,r["expression"] if r["expression"] is not None else -math.inf),reverse=True)
        top=ranked[0]
        vals=[r["expression"] for r in ranked if r["expression"] is not None]
        second=vals[1] if len(vals)>1 else None
        ratio=None
        if top["expression"] is not None:
            if second is not None and second>0:ratio=top["expression"]/second
            elif top["expression"]>0:ratio=float("inf")
        max_tau=max([r["tau_specificity"] for r in xs if r["tau_specificity"] is not None],default=None)
        out[gene]={
          "cell_source_n":len({r["source_name"] for r in xs}),
          "top_cell_type":top["cell_type"],
          "top_cell_compartment":top["compartment"],
          "top_cell_expression":top["expression"],
          "top_second_ratio":ratio,
          "max_cell_specificity":max_tau,
          "positive_celltype_n":sum((r["expression"] or 0)>0 for r in xs),
        }
    return out

def summarize_spatial(rows):
    by=defaultdict(list)
    for r in rows:by[r["gene_symbol"]].append(r)
    out={}
    for gene,xs in by.items():
        def score(r):
            if r["enrichment"] is not None:return (2,r["enrichment"])
            if r["expression"] is not None:return (1,r["expression"])
            if r["spot_fraction"] is not None:return (0,r["spot_fraction"])
            return (-1,-math.inf)
        top=max(xs,key=score)
        out[gene]={
          "spatial_source_n":len({r["source_name"] for r in xs}),
          "top_spatial_region":top["spatial_region"],
          "top_spatial_expression":top["expression"],
          "top_spatial_enrichment":top["enrichment"],
          "top_spatial_spot_fraction":top["spot_fraction"],
        }
    return out

def evidence_class(t,c,s):
    n=sum(x is not None for x in (t,c,s))
    if n==3:return "TISSUE_CELL_SPATIAL"
    if t is not None and c is not None:return "TISSUE_CELL"
    if c is not None and s is not None:return "CELL_SPATIAL"
    if t is not None and s is not None:return "TISSUE_SPATIAL"
    if n==1:return "SINGLE_LAYER"
    return "NO_LOCALIZATION"

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--candidates",type=Path,required=True,help="TSV containing gene_symbol")
    ap.add_argument("--tissue-manifest",type=Path)
    ap.add_argument("--cell-manifest",type=Path)
    ap.add_argument("--spatial-manifest",type=Path)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    candidates=[]
    for r in read_tsv(args.candidates):
        g=str(first(r,["gene_symbol","gene","Gene"]) or "").strip().upper()
        if g and g not in candidates:candidates.append(g)

    trows=load_tissue(args.tissue_manifest)
    crows=load_cell(args.cell_manifest)
    srows=load_spatial(args.spatial_manifest)
    ts=summarize_tissue(trows);cs=summarize_cell(crows);ss=summarize_spatial(srows)

    rows=[]
    for gene in candidates:
        t=ts.get(gene);c=cs.get(gene);s=ss.get(gene)
        row={"gene_symbol":gene}
        row.update(t or {
          "tissue_source_n":0,"top_tissue":"","top_tissue_expression":"",
          "top_tissue_share":"","max_tissue_specificity":""
        })
        row.update(c or {
          "cell_source_n":0,"top_cell_type":"","top_cell_compartment":"",
          "top_cell_expression":"","top_second_ratio":"","max_cell_specificity":"",
          "positive_celltype_n":0
        })
        row.update(s or {
          "spatial_source_n":0,"top_spatial_region":"","top_spatial_expression":"",
          "top_spatial_enrichment":"","top_spatial_spot_fraction":""
        })
        row["localization_evidence_class"]=evidence_class(t,c,s)
        rows.append(row)

    fields=[
      "gene_symbol",
      "tissue_source_n","top_tissue","top_tissue_expression","top_tissue_share","max_tissue_specificity",
      "cell_source_n","top_cell_type","top_cell_compartment","top_cell_expression","top_second_ratio","max_cell_specificity","positive_celltype_n",
      "spatial_source_n","top_spatial_region","top_spatial_expression","top_spatial_enrichment","top_spatial_spot_fraction",
      "localization_evidence_class"
    ]
    write_tsv(args.output,rows,fields)
    print(f"PASS candidates={len(rows)} tissue_rows={len(trows)} cell_rows={len(crows)} spatial_rows={len(srows)} output={args.output}")

if __name__=="__main__":main()
