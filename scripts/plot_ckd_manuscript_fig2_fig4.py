#!/usr/bin/env python3
"""Reproduce CKD Fig 2 EAS anchor MR and Fig 4 kidney tissue/cell from frozen 9-gene TSV.
Descriptive plots only. No causal conclusions, no external download, SHA guard mandatory.
"""
import argparse,csv,hashlib,math
from pathlib import Path

GENES=("ACP1","CPVL","F12","GSTA1","GSTA3","HLA-E","INHBC","SDCCAG8","UMOD")
def num(v,field,allow_na=True):
    if str(v or "").strip().upper() in ("","NA","NAN","NONE","NULL","MISSING"):
        if allow_na:return None
        raise ValueError("Missing "+field)
    try:r=float(v)
    except (TypeError,ValueError):raise ValueError("Invalid "+field)
    if not math.isfinite(r):raise ValueError("Nonfinite "+field)
    return r
def load(p):
    with p.open(encoding="utf-8-sig",newline="") as f:
        d=csv.DictReader(f,delimiter="\t")
        need={"gene_symbol","eas_egfr_wald_beta","eas_egfr_wald_se","eas_egfr_wald_p",
              "top_cell_type","kidney_eqtl_meta686_hits","kidney_eqtl_tubule356_hits",
              "kidney_eqtl_glomerulus303_hits","kidney_pqtl_maintext"}
        if not d.fieldnames or not need.issubset(d.fieldnames):raise ValueError("Missing required fields")
        rows=list(d)
    if len(rows)!=9 or set(r["gene_symbol"] for r in rows)!=set(GENES):
        raise ValueError("Not exact nine-gene matrix")
    return {r["gene_symbol"]:r for r in rows}
def build(source,out,expected_sha):
    actual=hashlib.sha256(source.read_bytes()).hexdigest()
    if actual.lower()!=expected_sha.lower():raise ValueError("SHA mismatch")
    if out.exists() and any(out.iterdir()):raise ValueError("Refusing overwrite")
    data=load(source); records=[]
    for g in GENES:
        r=data[g]
        beta=num(r["eas_egfr_wald_beta"],g+" beta")
        se=num(r["eas_egfr_wald_se"],g+" SE")
        p=num(r["eas_egfr_wald_p"],g+" P")
        if se is not None and se<=0:raise ValueError("Nonpositive SE")
        if p is not None and not 0<=p<=1:raise ValueError("Invalid P")
        if (beta is None)!=(se is None):raise ValueError("Incomplete effect/SE pair")
        hits=[]
        for key in ("kidney_eqtl_meta686_hits","kidney_eqtl_tubule356_hits","kidney_eqtl_glomerulus303_hits"):
            z=num(r[key],g+" "+key)
            if z is not None and (z<0 or z!=int(z)):raise ValueError("Invalid hit count")
            hits.append(int(z) if z is not None else None)
        records.append(dict(gene_symbol=g,beta=beta,se=se,p=p,cell=r["top_cell_type"] or "NOT_ASSESSED",
                            pqtl=r["kidney_pqtl_maintext"] or "NOT_ASSESSED",
                            meta686=hits[0],tubule356=hits[1],glomerulus303=hits[2],
                            source_sha256=actual))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out.mkdir(parents=True,exist_ok=True)
    fig,ax=plt.subplots(figsize=(9.5,6))
    for i,r in enumerate(records):
        if r["beta"] is not None:
            ax.errorbar(r["beta"],i,xerr=1.96*r["se"],fmt="o",capsize=3)
        else:
            ax.text(.02,i,"NA",transform=ax.get_yaxis_transform(),va="center")
    ax.axvline(0,ls="--",lw=1);ax.set_yticks(range(9),GENES);ax.invert_yaxis()
    ax.set_title("CKD Fig 2 | EAS eGFR anchor Wald estimates",loc="left")
    ax.set_xlabel("Wald effect estimate (source outcome and protein units)")
    ax.grid(axis="x",alpha=.2)
    fig.text(.12,.025,"Descriptive single-anchor MR; not independently validated causality. SHA "+actual[:14],fontsize=8)
    fig.subplots_adjust(left=.16,bottom=.15,top=.88)
    fig.savefig(out/"CKD_FIG2_EAS_ANCHOR_WALD.svg");plt.close(fig)
    fig,ax=plt.subplots(figsize=(11,6.6))
    y=list(range(9))
    for offset,col,label in [(-.22,"meta686","Kidney meta-eQTL"),(0,"tubule356","Tubule eQTL"),(.22,"glomerulus303","Glomerulus eQTL")]:
        xs=[r[col] if r[col] is not None else 0 for r in records]
        ax.barh([v+offset for v in y],xs,height=.19,label=label,alpha=.82)
    ax.set_yticks(y,GENES);ax.invert_yaxis();ax.set_xlabel("Significant-only source hit count (not comparable across resources)")
    ax.set_title("CKD Fig 4 | Kidney tissue QTL coverage",loc="left");ax.legend()
    ax.grid(axis="x",alpha=.2)
    for i,r in enumerate(records):
        ax.text(.99,i,("Cell: "+r["cell"])[:54],transform=ax.get_yaxis_transform(),ha="right",va="center",fontsize=7,
                bbox={"facecolor":"white","edgecolor":"none","alpha":.8})
    fig.text(.12,.022,"Zero hits may mean no significant hits or missing assay coverage; expression is not causal QTL.",fontsize=8)
    fig.subplots_adjust(left=.16,right=.98,bottom=.15,top=.88)
    fig.savefig(out/"CKD_FIG4_KIDNEY_TISSUE_QTL.svg");plt.close(fig)
    fields=["gene_symbol","beta","se","p","cell","pqtl","meta686","tubule356","glomerulus303","source_sha256"]
    with (out/"CKD_FIG2_FIG4_SOURCE.tsv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=fields,delimiter="\t");writer.writeheader()
        for row in records:writer.writerow({k:("NA" if row[k] is None else row[k]) for k in fields})
    return {"genes":9,"mr_estimates":sum(r["beta"] is not None for r in records),"sha256":actual}
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--input",required=True,type=Path);p.add_argument("--out",required=True,type=Path)
    p.add_argument("--sha256",required=True)
    a=p.parse_args();print(build(a.input,a.out,a.sha256))
