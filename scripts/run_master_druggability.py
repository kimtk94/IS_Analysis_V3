#!/usr/bin/env python3
"""Standardize druggability / target-translation evidence.

Input table may originate from Open Targets, DGIdb, ChEMBL, DrugBank exports,
clinical-trial curation, or manual review. No external API is called here.

Required: gene_symbol
Optional aliases:
drug_name, mechanism/action, direction, indication, phase/status,
approved, target_class, source, evidence_level.

Direction logic is explicit: if a causal direction file is supplied, the
script evaluates whether the therapeutic action is directionally compatible
with the MR effect. It does not equate genetic lifelong exposure with drug
intervention.
"""
from __future__ import annotations
import argparse,csv,gzip,math
from pathlib import Path
from collections import defaultdict

def read_tsv(path):
    path=Path(path);op=gzip.open if path.suffix==".gz" else open
    with op(path,"rt",encoding="utf-8",newline="") as f:
        yield from csv.DictReader(f,delimiter="\t")

def write(path,rows,fields):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore")
        w.writeheader();w.writerows(rows)

def first(r,names):
    for n in names:
        if n in r and str(r[n]).strip()!="":return r[n]
    return None

def fnum(x):
    try:v=float(str(x).strip())
    except Exception:return None
    return v if math.isfinite(v) else None

def boolish(x):
    s=str(x or "").strip().lower()
    return int(s in {"1","true","yes","y","approved","marketed"})

def action_sign(action):
    s=str(action or "").lower()
    if any(k in s for k in ["inhibit","inhibitor","antagon","block","degrad","silenc","knockdown"]):return -1
    if any(k in s for k in ["agon","activat","stimulat","enhanc","replacement"]):return 1
    return None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--drug-table",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--causal-direction",type=Path,help="gene_symbol,beta[,phenotype]")
    args=ap.parse_args()

    causal={}
    if args.causal_direction:
        for r in read_tsv(args.causal_direction):
            g=str(first(r,["gene_symbol","gene","Gene"]) or "").strip().upper()
            b=fnum(first(r,["beta","wald_beta","anchor_wald_beta","ivw_beta"]))
            if g and b is not None:causal[g]=b

    rows=[]
    for r in read_tsv(args.drug_table):
        g=str(first(r,["gene_symbol","gene","Gene","target_gene"]) or "").strip().upper()
        if not g:continue
        drug=str(first(r,["drug_name","drug","compound"]) or "").strip()
        action=str(first(r,["mechanism","action","mode_of_action","direction"]) or "").strip()
        indication=str(first(r,["indication","disease"]) or "").strip()
        phase=str(first(r,["phase","clinical_phase","status"]) or "").strip()
        source=str(first(r,["source","database"]) or "").strip()
        approved=boolish(first(r,["approved","is_approved","status"]))
        beta=causal.get(g)
        asign=action_sign(action)
        compatible=""
        # beta>0 means genetically higher target is associated with higher disease risk;
        # inhibition is directionally compatible. beta<0 implies activation could be.
        if beta not in (None,0) and asign in (-1,1):
            compatible=int((beta>0 and asign==-1) or (beta<0 and asign==1))
        rows.append({
          "gene_symbol":g,"drug_name":drug,"mechanism":action,"action_sign":asign if asign is not None else "",
          "indication":indication,"clinical_phase_or_status":phase,"approved":approved,
          "source":source,"causal_beta":beta if beta is not None else "",
          "directionally_compatible":compatible,
          "translation_note":"genetic direction is supportive only; pharmacologic equivalence is not assumed"
        })

    by=defaultdict(list)
    for r in rows:by[r["gene_symbol"]].append(r)
    summary=[]
    for g,xs in sorted(by.items()):
        compatible=[x for x in xs if x["directionally_compatible"]==1]
        approved=[x for x in xs if x["approved"]==1]
        summary.append({
          "gene_symbol":g,
          "drug_record_n":len(xs),
          "approved_drug_n":len(approved),
          "directionally_compatible_drug_n":len(compatible),
          "has_approved_drug":int(bool(approved)),
          "has_directionally_compatible_drug":int(bool(compatible)),
          "top_drug":next((x["drug_name"] for x in compatible if x["drug_name"]),next((x["drug_name"] for x in xs if x["drug_name"]),"")),
          "translation_tier":"A" if compatible and approved else "B" if compatible else "C" if xs else "U",
        })

    fields=["gene_symbol","drug_record_n","approved_drug_n","directionally_compatible_drug_n",
            "has_approved_drug","has_directionally_compatible_drug","top_drug","translation_tier"]
    write(args.output,summary,fields)
    detail=args.output.with_name(args.output.stem+"_DETAIL.tsv")
    dfields=["gene_symbol","drug_name","mechanism","action_sign","indication","clinical_phase_or_status","approved",
             "source","causal_beta","directionally_compatible","translation_note"]
    write(detail,rows,dfields)
    print(f"PASS genes={len(summary)} records={len(rows)} output={args.output}")

if __name__=="__main__":main()
