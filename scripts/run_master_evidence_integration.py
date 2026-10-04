#!/usr/bin/env python3
"""Final MASTER evidence integration and tier assignment.

The script consumes a manifest of stage-level summary files:
stage    file    key_column

Known stage labels:
mr, robustness, coloc, susie, ancestry, platform, transcriptomics,
phenotype, risk_factor, individual, localization, phewas, druggability.

Scoring is intentionally transparent and conservative. Hard conflicts can
downgrade otherwise strong candidates.
"""
from __future__ import annotations
import argparse,csv,gzip,math
from pathlib import Path
from collections import defaultdict

def read_tsv(path):
    path=Path(path);op=gzip.open if path.suffix==".gz" else open
    with op(path,"rt",encoding="utf-8",newline="") as f:yield from csv.DictReader(f,delimiter="\t")

def fnum(x):
    try:v=float(str(x).strip())
    except Exception:return None
    return v if math.isfinite(v) else None

def first(r,names):
    for n in names:
        if n in r and str(r[n]).strip()!="":return r[n]
    return None

def truth(x):
    return str(x or "").strip().lower() in {"1","true","yes","y"}

def stage_signal(stage,r):
    if stage=="mr":
        p=fnum(first(r,["p","wald_p","anchor_wald_p","ivw_p"]))
        return (2 if p is not None and p<0.05 else 0, "")
    if stage=="robustness":
        qp=fnum(r.get("cochran_q_p")); ep=fnum(r.get("egger_intercept_p"))
        score=1 if (qp is None or qp>=0.05) and (ep is None or ep>=0.05) else 0
        return score, "" if score else "MR_ROBUSTNESS"
    if stage=="coloc":
        h4=fnum(first(r,["PP.H4","max_PP.H4"]))
        return (2 if h4 is not None and h4>=0.8 else 0, "")
    if stage=="susie":
        h4=fnum(first(r,["max_PP.H4","PP.H4"]))
        status=str(r.get("comparison_status",""))
        if "nonconverged" in status:return (0,"SUSIE_UNRESOLVED")
        return (2 if h4 is not None and h4>=0.8 else 0, "")
    if stage=="ancestry":
        t=str(r.get("replication_tier",""))
        return (2 if t=="A" else 1 if t in {"B","C"} else 0, "ANCESTRY_DIRECTION" if t=="D" else "")
    if stage=="platform":
        t=str(r.get("platform_replication_tier",""))
        return (2 if t=="A" else 1 if t in {"B","C"} else 0, "PLATFORM_DIRECTION" if t=="D" else "")
    if stage=="transcriptomics":
        t=str(r.get("transcript_evidence_tier",""))
        return (2 if t=="T1" else 1 if t in {"T2","T3"} else 0, "TRANSCRIPT_CONFLICT" if t=="CONFLICT" else "")
    if stage=="phenotype":
        n=fnum(first(r,["n_mr_p05","n_phenotypes"]))
        return (1 if n is not None and n>0 else 0, "")
    if stage=="risk_factor":
        n=sum(1 for x in [r] if truth(x.get("mechanism_support")))
        return (1 if n else 0, "")
    if stage=="individual":
        p=fnum(r.get("p"));status=str(r.get("status",""))
        return (2 if status=="ok" and p is not None and p<0.05 else 1 if status=="ok" else 0, "")
    if stage=="localization":
        c=str(r.get("localization_evidence_class",""))
        return (2 if c=="TISSUE_CELL_SPATIAL" else 1 if c not in {"","NO_LOCALIZATION"} else 0, "")
    if stage=="phewas":
        adverse=fnum(r.get("phewas_significant_adverse_n"))
        sig=fnum(r.get("phewas_fdr_significant_n"))
        return (1 if sig is not None and sig>0 else 0, "PHEWAS_SAFETY" if adverse is not None and adverse>0 else "")
    if stage=="druggability":
        t=str(r.get("translation_tier",""))
        return (2 if t=="A" else 1 if t in {"B","C"} else 0, "")
    return 0,""

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--candidates",type=Path,required=True)
    ap.add_argument("--manifest",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    genes=[]
    for r in read_tsv(args.candidates):
        g=str(first(r,["gene_symbol","gene","Gene"]) or "").strip().upper()
        if g and g not in genes:genes.append(g)

    evidence=defaultdict(lambda:defaultdict(list))
    for m in read_tsv(args.manifest):
        stage=str(m["stage"]).strip()
        path=Path(m["file"])
        keycol=str(m.get("key_column","gene_symbol")).strip() or "gene_symbol"
        for r in read_tsv(path):
            g=str(r.get(keycol,"")).strip().upper()
            if g:evidence[g][stage].append(r)

    rows=[]
    for g in genes:
        total=0;conflicts=[];stage_scores={}
        for stage,rs in evidence[g].items():
            best=0
            for r in rs:
                sc,cf=stage_signal(stage,r)
                best=max(best,sc)
                if cf:conflicts.append(cf)
            stage_scores[stage]=best
            total+=best

        core=(stage_scores.get("mr",0)>=2 and stage_scores.get("coloc",0)>=2)
        replication=max(stage_scores.get("ancestry",0),stage_scores.get("platform",0))
        biological=max(stage_scores.get("transcriptomics",0),stage_scores.get("localization",0))
        clinical=stage_scores.get("individual",0)
        hard_conflict=any(x in conflicts for x in ["ANCESTRY_DIRECTION","PLATFORM_DIRECTION","TRANSCRIPT_CONFLICT"])

        if core and replication>=1 and biological>=1 and clinical>=1 and not hard_conflict:
            tier="Tier 1"
        elif core and (replication>=1 or biological>=1) and not hard_conflict:
            tier="Tier 2"
        elif stage_scores.get("mr",0)>=2:
            tier="Tier 3"
        else:
            tier="Exploratory"
        if hard_conflict and tier=="Tier 1":tier="Tier 2"
        if hard_conflict and tier=="Tier 2":tier="Tier 3"

        row={"gene_symbol":g,"final_tier":tier,"evidence_score":total,
             "conflict_flags":";".join(sorted(set(conflicts)))}
        for s in ["mr","robustness","coloc","susie","ancestry","platform","transcriptomics","phenotype",
                  "risk_factor","individual","localization","phewas","druggability"]:
            row[f"score_{s}"]=stage_scores.get(s,0)
        rows.append(row)

    rank={"Tier 1":0,"Tier 2":1,"Tier 3":2,"Exploratory":3}
    rows.sort(key=lambda r:(rank[r["final_tier"]],-r["evidence_score"],r["gene_symbol"]))
    fields=["gene_symbol","final_tier","evidence_score","conflict_flags"]+[f"score_{s}" for s in
      ["mr","robustness","coloc","susie","ancestry","platform","transcriptomics","phenotype","risk_factor",
       "individual","localization","phewas","druggability"]]
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n")
        w.writeheader();w.writerows(rows)
    print(f"PASS candidates={len(rows)} output={args.output}")

if __name__=="__main__":main()
