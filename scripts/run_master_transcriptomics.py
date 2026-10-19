#!/usr/bin/env python3
"""Integrate transcriptomic evidence (SMR/HEIDI/eQTL coloc) with protein MR.

This stage does not re-implement the external SMR executable. It standardizes
its result table and combines it with protein-level causal evidence.

Typical inputs
--------------
protein MR:
  gene_symbol, phenotype, beta/se/p (or compatible MASTER aliases)

SMR output:
  gene/probe identifier, b_SMR, se_SMR, p_SMR, p_HEIDI
  Common aliases are recognized.

Optional eQTL coloc:
  locus/gene_symbol plus PP.H4 / PP.H3.

Evidence classes
----------------
T1: SMR significant + HEIDI non-significant + eQTL coloc H4>=threshold
    + protein/SMR causal-effect direction concordant.
T2: SMR significant + HEIDI non-significant + at least one of
    (coloc support, direction concordance).
T3: SMR significant but orthogonal evidence incomplete.
CONFLICT: significant HEIDI heterogeneity or protein/SMR direction discordance.
U: insufficient transcriptomic evidence.

HEIDI P>=threshold is treated as "no evidence of heterogeneity", not proof of
a shared causal variant. Strong shared-variant evidence should come from coloc.
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
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore")
        w.writeheader(); w.writerows(rows)


def first_present(row, names):
    for n in names:
        if n in row and str(row[n]).strip()!="":
            return row[n]
    return None


def load_probe_map(path: Path | None):
    if path is None:
        return {}
    out={}
    for r in read_tsv(path):
        probe=str(first_present(r,["probe_id","ProbeID","probe"]) or "").strip()
        gene=str(first_present(r,["gene_symbol","Gene","gene"]) or "").strip()
        if probe and gene:
            out[probe]=gene
    return out


def normalize_protein_mr(path: Path):
    out={}
    for r in read_tsv(path):
        gene=str(r.get("gene_symbol","")).strip()
        pheno=str(r.get("phenotype","")).strip()
        if not gene or not pheno:
            continue
        beta=fnum(first_present(r,["beta","wald_beta","anchor_wald_beta","ivw_beta"]))
        se=fnum(first_present(r,["se","wald_se","anchor_wald_se","ivw_se"]))
        p=fnum(first_present(r,["p","wald_p","anchor_wald_p","ivw_p"]))
        out[(gene,pheno)]={"beta":beta,"se":se,"p":p}
    return out


def normalize_smr(path: Path, probe_map: dict, phenotype: str):
    out={}
    for r in read_tsv(path):
        probe=str(first_present(r,["ProbeID","probe_id","probe","probeID"]) or "").strip()
        gene=str(first_present(r,["gene_symbol","Gene","gene","GENE"]) or "").strip()
        if not gene and probe:
            gene=probe_map.get(probe,"")
        if not gene:
            continue
        beta=fnum(first_present(r,["b_SMR","beta_smr","SMR_beta","beta"]))
        se=fnum(first_present(r,["se_SMR","se_smr","SMR_se","se"]))
        p=fnum(first_present(r,["p_SMR","p_smr","SMR_p","p"]))
        heidi=fnum(first_present(r,["p_HEIDI","p_heidi","HEIDI_p","heidi_p"]))
        nsnp=fnum(first_present(r,["nsnp_HEIDI","nsnp_heidi","HEIDI_nsnp"]))
        key=(gene,phenotype)
        if key in out:
            # Keep strongest SMR result if multiple probes map to one gene.
            old=out[key]
            if old["p_smr"] is not None and (p is None or old["p_smr"]<=p):
                continue
        out[key]={
            "probe_id":probe,"beta_smr":beta,"se_smr":se,
            "p_smr":p,"p_heidi":heidi,"nsnp_heidi":nsnp,
        }
    return out


def normalize_coloc(path: Path | None):
    if path is None:
        return {}
    out={}
    for r in read_tsv(path):
        locus=str(first_present(r,["locus","gene_symbol","Gene","gene"]) or "").strip()
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


def classify(protein, smr, h4, smr_p_threshold, heidi_p_threshold, h4_threshold):
    if smr is None or smr["p_smr"] is None:
        return "U"

    smr_sig=smr["p_smr"]<smr_p_threshold
    if not smr_sig:
        return "U"

    heidi=smr["p_heidi"]
    heidi_conflict=(heidi is not None and heidi<heidi_p_threshold)
    heidi_ok=(heidi is not None and heidi>=heidi_p_threshold)
    coloc_ok=(h4 is not None and h4>=h4_threshold)

    direction=None
    if protein and protein["beta"] is not None and smr["beta_smr"] is not None:
        if protein["beta"]==0 or smr["beta_smr"]==0:
            direction=None
        else:
            direction=(protein["beta"]*smr["beta_smr"])>0

    if heidi_conflict or direction is False:
        return "CONFLICT"
    if heidi_ok and coloc_ok and direction is True:
        return "T1"
    if heidi_ok and (coloc_ok or direction is True):
        return "T2"
    return "T3"


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--protein-mr",type=Path,required=True)
    ap.add_argument("--smr",type=Path,required=True)
    ap.add_argument("--phenotype",required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--eqtl-coloc",type=Path)
    ap.add_argument("--probe-map",type=Path)
    ap.add_argument("--smr-p-threshold",type=float,default=0.05)
    ap.add_argument("--heidi-p-threshold",type=float,default=0.05)
    ap.add_argument("--coloc-h4-threshold",type=float,default=0.80)
    args=ap.parse_args()

    probe_map=load_probe_map(args.probe_map)
    pmr=normalize_protein_mr(args.protein_mr)
    smr=normalize_smr(args.smr,probe_map,args.phenotype)
    coloc=normalize_coloc(args.eqtl_coloc)

    keys=sorted(set(k for k in pmr if k[1]==args.phenotype) | set(smr))
    rows=[]
    for key in keys:
        gene,pheno=key
        p=pmr.get(key)
        s=smr.get(key)
        h4=coloc.get(gene,{}).get("PP.H4")
        h3=coloc.get(gene,{}).get("PP.H3")

        direction=""
        if p and s and p["beta"] is not None and s["beta_smr"] is not None:
            if p["beta"]==0 or s["beta_smr"]==0:
                direction="zero_effect"
            elif p["beta"]*s["beta_smr"]>0:
                direction="concordant"
            else:
                direction="discordant"

        heidi_status=""
        if s and s["p_heidi"] is not None:
            heidi_status="no_evidence_of_heterogeneity" if s["p_heidi"]>=args.heidi_p_threshold else "heterogeneity_signal"

        evidence=classify(
            p,s,h4,args.smr_p_threshold,args.heidi_p_threshold,args.coloc_h4_threshold
        )

        rows.append({
            "gene_symbol":gene,
            "phenotype":pheno,
            "protein_mr_beta":p["beta"] if p else "",
            "protein_mr_se":p["se"] if p else "",
            "protein_mr_p":p["p"] if p else "",
            "smr_probe_id":s["probe_id"] if s else "",
            "smr_beta":s["beta_smr"] if s else "",
            "smr_se":s["se_smr"] if s else "",
            "smr_p":s["p_smr"] if s else "",
            "heidi_p":s["p_heidi"] if s else "",
            "heidi_nsnp":s["nsnp_heidi"] if s else "",
            "heidi_status":heidi_status,
            "eqtl_coloc_PP.H3":h3 if h3 is not None else "",
            "eqtl_coloc_PP.H4":h4 if h4 is not None else "",
            "protein_vs_expression_direction":direction,
            "transcript_evidence_tier":evidence,
        })

    rank={"T1":0,"T2":1,"T3":2,"CONFLICT":3,"U":4}
    rows.sort(key=lambda r:(rank[r["transcript_evidence_tier"]],
                            r["smr_p"] if isinstance(r["smr_p"],float) else 1.0,
                            r["gene_symbol"]))

    fields=[
        "gene_symbol","phenotype",
        "protein_mr_beta","protein_mr_se","protein_mr_p",
        "smr_probe_id","smr_beta","smr_se","smr_p",
        "heidi_p","heidi_nsnp","heidi_status",
        "eqtl_coloc_PP.H3","eqtl_coloc_PP.H4",
        "protein_vs_expression_direction","transcript_evidence_tier"
    ]
    write_tsv(args.output,rows,fields)
    counts={k:sum(r["transcript_evidence_tier"]==k for r in rows) for k in rank}
    print("PASS "+" ".join(f"{k}={v}" for k,v in counts.items())+f" output={args.output}")


if __name__=="__main__":
    main()
