#!/usr/bin/env python3
"""Cross-ancestry replication summary for MASTER pipeline.

Compares ancestry-specific MR estimates without requiring the same lead SNP.
Optional coloc summaries can strengthen replication classification.

Default replication classes
A: direction-concordant + replication P<0.05 + replication PP.H4>=0.80
B: direction-concordant + replication P<0.05
C: direction-concordant but replication P>=0.05
D: opposite causal-effect direction
U: unavailable / non-estimable

The script also reports a two-study fixed-effect meta estimate and a
between-ancestry effect-difference test. These are secondary summaries;
replication is not defined by meta-analysis significance.
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


def pnorm(z):
    return math.erfc(abs(z)/math.sqrt(2.0))


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


def normalize_mr(path: Path):
    out={}
    for r in read_tsv(path):
        protein=str(r.get("protein_id","")).strip()
        gene=str(r.get("gene_symbol","")).strip()
        phenotype=str(r.get("phenotype","")).strip()
        key_id=protein or gene
        if not key_id or not phenotype:
            continue
        beta=fnum(first_present(r,["beta","wald_beta","anchor_wald_beta","ivw_beta"]))
        se=fnum(first_present(r,["se","wald_se","anchor_wald_se","ivw_se"]))
        p=fnum(first_present(r,["p","wald_p","anchor_wald_p","ivw_p"]))
        fdr=fnum(first_present(r,["fdr_bh","wald_fdr_bh"]))
        nsnp=fnum(first_present(r,["n_instruments","n_harmonized_instruments","nsnp"]))
        key=(key_id,phenotype)
        if key in out:
            raise ValueError(f"duplicate MR key in {path}: {key}")
        out[key]={
            "protein_id":protein,
            "gene_symbol":gene,
            "phenotype":phenotype,
            "beta":beta,"se":se,"p":p,"fdr":fdr,"nsnp":nsnp,
        }
    return out


def normalize_coloc(path: Path | None):
    if path is None:
        return {}
    out={}
    for r in read_tsv(path):
        locus=str(first_present(r,["locus","gene_symbol","protein_id"]) or "").strip()
        if not locus:
            continue
        h4=fnum(first_present(r,["PP.H4","PP.H4.abf","max_PP.H4"]))
        h3=fnum(first_present(r,["PP.H3","PP.H3.abf","PP.H3_at_best_pair"]))
        p12=fnum(r.get("p12"))
        # Prefer default prior p12=1e-5 if multiple rows exist.
        score=(0 if p12 is None else abs(math.log10(p12)-math.log10(1e-5)))
        prev=out.get(locus)
        if prev is None or score < prev["_prior_distance"]:
            out[locus]={"PP.H4":h4,"PP.H3":h3,"_prior_distance":score}
    return out


def fixed_meta(b1,s1,b2,s2):
    if None in (b1,s1,b2,s2) or s1<=0 or s2<=0:
        return None
    w1=1/(s1*s1); w2=1/(s2*s2)
    beta=(w1*b1+w2*b2)/(w1+w2)
    se=math.sqrt(1/(w1+w2))
    q=w1*(b1-beta)**2+w2*(b2-beta)**2
    return {
        "beta":beta,"se":se,"p":pnorm(beta/se),
        "q":q,"q_p":pnorm(math.sqrt(max(q,0.0))),
    }


def diff_test(b1,s1,b2,s2):
    if None in (b1,s1,b2,s2) or s1<=0 or s2<=0:
        return (None,None)
    sed=math.sqrt(s1*s1+s2*s2)
    z=(b1-b2)/sed
    return z,pnorm(z)


def classify(d, r, rep_h4, p_threshold, h4_threshold):
    if d is None or r is None:
        return "U"
    b1,b2=d["beta"],r["beta"]
    if b1 is None or b2 is None or b1==0 or b2==0:
        return "U"
    concordant=(b1*b2)>0
    if not concordant:
        return "D"
    rp=r["p"]
    if rp is not None and rp < p_threshold:
        if rep_h4 is not None and rep_h4 >= h4_threshold:
            return "A"
        return "B"
    return "C"


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--discovery-mr",type=Path,required=True)
    ap.add_argument("--replication-mr",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--discovery-ancestry",required=True)
    ap.add_argument("--replication-ancestry",required=True)
    ap.add_argument("--discovery-coloc",type=Path)
    ap.add_argument("--replication-coloc",type=Path)
    ap.add_argument("--replication-p-threshold",type=float,default=0.05)
    ap.add_argument("--coloc-h4-threshold",type=float,default=0.80)
    args=ap.parse_args()

    disc=normalize_mr(args.discovery_mr)
    repl=normalize_mr(args.replication_mr)
    dc=normalize_coloc(args.discovery_coloc)
    rc=normalize_coloc(args.replication_coloc)

    rows=[]
    for key in sorted(set(disc)|set(repl)):
        d=disc.get(key); r=repl.get(key)
        base=d or r
        gene=(d or {}).get("gene_symbol") or (r or {}).get("gene_symbol") or ""
        protein=(d or {}).get("protein_id") or (r or {}).get("protein_id") or ""
        locus=gene or protein or key[0]
        dh4=dc.get(locus,{}).get("PP.H4")
        rh4=rc.get(locus,{}).get("PP.H4")

        tier=classify(d,r,rh4,args.replication_p_threshold,args.coloc_h4_threshold)
        direction=""
        if d and r and d["beta"] is not None and r["beta"] is not None:
            if d["beta"]==0 or r["beta"]==0:
                direction="zero_effect"
            elif d["beta"]*r["beta"]>0:
                direction="concordant"
            else:
                direction="discordant"

        meta=fixed_meta(
            d["beta"] if d else None,d["se"] if d else None,
            r["beta"] if r else None,r["se"] if r else None
        )
        dz,dp=diff_test(
            d["beta"] if d else None,d["se"] if d else None,
            r["beta"] if r else None,r["se"] if r else None
        )

        rows.append({
            "protein_id":protein,
            "gene_symbol":gene,
            "phenotype":base["phenotype"],
            "discovery_ancestry":args.discovery_ancestry,
            "replication_ancestry":args.replication_ancestry,
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
            "effect_difference_z":dz if dz is not None else "",
            "effect_difference_p":dp if dp is not None else "",
            "fixed_meta_beta":meta["beta"] if meta else "",
            "fixed_meta_se":meta["se"] if meta else "",
            "fixed_meta_p":meta["p"] if meta else "",
            "fixed_meta_q":meta["q"] if meta else "",
            "fixed_meta_q_p":meta["q_p"] if meta else "",
            "replication_tier":tier,
        })

    tier_rank={"A":0,"B":1,"C":2,"D":3,"U":4}
    rows.sort(key=lambda x:(tier_rank[x["replication_tier"]],
                            x["replication_p"] if isinstance(x["replication_p"],float) else 1.0,
                            x["gene_symbol"],x["protein_id"],x["phenotype"]))

    fields=[
        "protein_id","gene_symbol","phenotype",
        "discovery_ancestry","replication_ancestry",
        "discovery_beta","discovery_se","discovery_p","discovery_fdr",
        "discovery_n_instruments","discovery_PP.H4",
        "replication_beta","replication_se","replication_p","replication_fdr",
        "replication_n_instruments","replication_PP.H4",
        "direction","effect_difference_z","effect_difference_p",
        "fixed_meta_beta","fixed_meta_se","fixed_meta_p","fixed_meta_q","fixed_meta_q_p",
        "replication_tier"
    ]
    write_tsv(args.output,rows,fields)

    counts={k:sum(r["replication_tier"]==k for r in rows) for k in ["A","B","C","D","U"]}
    print("PASS "+" ".join(f"tier_{k}={v}" for k,v in counts.items())+f" output={args.output}")


if __name__=="__main__":
    main()
