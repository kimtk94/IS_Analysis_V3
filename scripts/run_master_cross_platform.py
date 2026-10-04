#!/usr/bin/env python3
"""Cross-platform proteomic replication for MASTER pipeline.

Primary use:
Olink discovery MR -> SomaScan replication MR (or vice versa).

Replication is evaluated at the mapped protein/gene level, not by requiring
the same sentinel SNP. This is intentional because assay platforms can have
different cis-pQTL instruments.

Tier definitions:
A: direction concordant + replication P<0.05 + replication PP.H4>=0.80
B: direction concordant + replication P<0.05
C: direction concordant, replication P>=0.05
D: direction discordant
U: unavailable / non-estimable

Effect-size meta-analysis is deliberately omitted because platform-specific
protein scales (e.g. NPX vs aptamer RFU/standardized traits) are not assumed
to be directly commensurate.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import math
from pathlib import Path


def fnum(x):
    if x is None:
        return None
    try:
        v=float(str(x).strip())
    except Exception:
        return None
    return v if math.isfinite(v) else None


def read_tsv(path: Path):
    op=gzip.open if path.suffix==".gz" else open
    with op(path,"rt",encoding="utf-8",newline="") as fh:
        yield from csv.DictReader(fh,delimiter="\t")


def write_tsv(path: Path, rows: list[dict], fields: list[str]):
    path.parent.mkdir(parents=True,exist_ok=True)
    op=gzip.open if path.suffix==".gz" else open
    with op(path,"wt",encoding="utf-8",newline="") as fh:
        w=csv.DictWriter(
            fh,fieldnames=fields,delimiter="\t",lineterminator="\n",
            extrasaction="ignore"
        )
        w.writeheader(); w.writerows(rows)


def first_present(row, names):
    for n in names:
        if n in row and str(row[n]).strip()!="":
            return row[n]
    return None


def normalize_mr(path: Path, id_column: str | None):
    rows={}
    for r in read_tsv(path):
        gene=str(r.get("gene_symbol","")).strip()
        protein=str(r.get("protein_id","")).strip()
        phenotype=str(r.get("phenotype","")).strip()
        assay=str(r.get("assay_id","")).strip()

        if id_column:
            raw=str(r.get(id_column,"")).strip()
        else:
            raw=gene or protein or assay
        if not raw or not phenotype:
            continue

        beta=fnum(first_present(r,["beta","wald_beta","anchor_wald_beta","ivw_beta"]))
        se=fnum(first_present(r,["se","wald_se","anchor_wald_se","ivw_se"]))
        p=fnum(first_present(r,["p","wald_p","anchor_wald_p","ivw_p"]))
        fdr=fnum(first_present(r,["fdr_bh","wald_fdr_bh"]))
        nsnp=fnum(first_present(r,["n_instruments","n_harmonized_instruments","nsnp"]))
        key=(raw,phenotype)
        if key in rows:
            raise ValueError(f"duplicate MR key in {path}: {key}")
        rows[key]={
            "raw_id":raw,
            "gene_symbol":gene,
            "protein_id":protein,
            "assay_id":assay,
            "phenotype":phenotype,
            "beta":beta,"se":se,"p":p,"fdr":fdr,"nsnp":nsnp,
        }
    return rows


def load_mapping(path: Path | None):
    """Return discovery_id -> replication_id and optional canonical gene.

    Supported columns:
    discovery_id, replication_id
    optional: gene_symbol
    """
    if path is None:
        return {}
    out={}
    for r in read_tsv(path):
        d=str(r.get("discovery_id","")).strip()
        q=str(r.get("replication_id","")).strip()
        if not d or not q:
            continue
        if d in out:
            raise ValueError(f"duplicate discovery_id in mapping: {d}")
        out[d]={
            "replication_id":q,
            "gene_symbol":str(r.get("gene_symbol","")).strip(),
        }
    return out


def normalize_coloc(path: Path | None):
    if path is None:
        return {}
    out={}
    for r in read_tsv(path):
        locus=str(first_present(r,["locus","gene_symbol","protein_id","assay_id"]) or "").strip()
        if not locus:
            continue
        h4=fnum(first_present(r,["PP.H4","PP.H4.abf","max_PP.H4"]))
        h3=fnum(first_present(r,["PP.H3","PP.H3.abf","PP.H3_at_best_pair"]))
        p12=fnum(r.get("p12"))
        score=0 if p12 is None else abs(math.log10(p12)-math.log10(1e-5))
        prev=out.get(locus)
        if prev is None or score<prev["_prior_distance"]:
            out[locus]={"PP.H4":h4,"PP.H3":h3,"_prior_distance":score}
    return out


def tier(d, r, repl_h4, p_threshold, h4_threshold):
    if d is None or r is None:
        return "U"
    if d["beta"] is None or r["beta"] is None or d["beta"]==0 or r["beta"]==0:
        return "U"
    if d["beta"]*r["beta"]<0:
        return "D"
    if r["p"] is not None and r["p"]<p_threshold:
        if repl_h4 is not None and repl_h4>=h4_threshold:
            return "A"
        return "B"
    return "C"


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--discovery-mr",type=Path,required=True)
    ap.add_argument("--replication-mr",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--discovery-platform",required=True)
    ap.add_argument("--replication-platform",required=True)
    ap.add_argument("--mapping",type=Path)
    ap.add_argument("--discovery-id-column")
    ap.add_argument("--replication-id-column")
    ap.add_argument("--discovery-coloc",type=Path)
    ap.add_argument("--replication-coloc",type=Path)
    ap.add_argument("--replication-p-threshold",type=float,default=0.05)
    ap.add_argument("--coloc-h4-threshold",type=float,default=0.80)
    args=ap.parse_args()

    disc=normalize_mr(args.discovery_mr,args.discovery_id_column)
    repl=normalize_mr(args.replication_mr,args.replication_id_column)
    mapping=load_mapping(args.mapping)
    dc=normalize_coloc(args.discovery_coloc)
    rc=normalize_coloc(args.replication_coloc)

    if mapping:
        pairs=[]
        for (did,pheno),d in disc.items():
            m=mapping.get(did)
            if not m:
                pairs.append((did,None,pheno,d,None,m))
                continue
            rid=m["replication_id"]
            pairs.append((did,rid,pheno,d,repl.get((rid,pheno)),m))
    else:
        # Without explicit mapping, matching identifier must be identical.
        allkeys=sorted(set(disc)|set(repl))
        pairs=[(k[0],k[0],k[1],disc.get(k),repl.get(k),None) for k in allkeys]

    rows=[]
    for did,rid,pheno,d,r,m in pairs:
        gene=""
        if m and m.get("gene_symbol"):
            gene=m["gene_symbol"]
        elif d and d.get("gene_symbol"):
            gene=d["gene_symbol"]
        elif r and r.get("gene_symbol"):
            gene=r["gene_symbol"]

        dlocus=(d.get("gene_symbol") if d else "") or did or ""
        rlocus=(r.get("gene_symbol") if r else "") or rid or ""
        dh4=dc.get(dlocus,{}).get("PP.H4")
        rh4=rc.get(rlocus,{}).get("PP.H4")

        status=tier(d,r,rh4,args.replication_p_threshold,args.coloc_h4_threshold)
        direction=""
        if d and r and d["beta"] is not None and r["beta"] is not None:
            if d["beta"]==0 or r["beta"]==0:
                direction="zero_effect"
            elif d["beta"]*r["beta"]>0:
                direction="concordant"
            else:
                direction="discordant"

        rows.append({
            "gene_symbol":gene,
            "phenotype":pheno,
            "discovery_platform":args.discovery_platform,
            "replication_platform":args.replication_platform,
            "discovery_id":did,
            "replication_id":rid or "",
            "discovery_beta":d["beta"] if d else "",
            "discovery_se":d["se"] if d else "",
            "discovery_p":d["p"] if d else "",
            "discovery_fdr":d["fdr"] if d else "",
            "discovery_n_instruments":d["nsnp"] if d else "",
            "discovery_PP.H4":dh4 if dh4 is not None else "",
            "replication_beta":r["beta"] if r else "",
            "replication_se":r["se"] if r else "",
            "replication_p":r["p"] if r else "",
            "replication_fdr":r["fdr"] if r else "",
            "replication_n_instruments":r["nsnp"] if r else "",
            "replication_PP.H4":rh4 if rh4 is not None else "",
            "direction":direction,
            "platform_replication_tier":status,
        })

    rank={"A":0,"B":1,"C":2,"D":3,"U":4}
    rows.sort(key=lambda x:(
        rank[x["platform_replication_tier"]],
        x["replication_p"] if isinstance(x["replication_p"],float) else 1.0,
        x["gene_symbol"],x["phenotype"]
    ))

    fields=[
        "gene_symbol","phenotype",
        "discovery_platform","replication_platform",
        "discovery_id","replication_id",
        "discovery_beta","discovery_se","discovery_p","discovery_fdr",
        "discovery_n_instruments","discovery_PP.H4",
        "replication_beta","replication_se","replication_p","replication_fdr",
        "replication_n_instruments","replication_PP.H4",
        "direction","platform_replication_tier"
    ]
    write_tsv(args.output,rows,fields)

    counts={k:sum(r["platform_replication_tier"]==k for r in rows) for k in ["A","B","C","D","U"]}
    print("PASS "+" ".join(f"tier_{k}={v}" for k,v in counts.items())+f" output={args.output}")


if __name__=="__main__":
    main()
