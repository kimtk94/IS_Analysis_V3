#!/usr/bin/env python3
"""Descriptive CKD Figure 1: provenance-locked analysis workflow, no invented attrition."""
import argparse,csv,hashlib,json
from pathlib import Path
def render(gene_matrix,out,sha):
    if hashlib.sha256(gene_matrix.read_bytes()).hexdigest()!=sha:raise ValueError("CKD source checksum mismatch")
    with gene_matrix.open(newline="",encoding="utf-8-sig") as f:
        rows=list(csv.DictReader(f,delimiter="\\t"))
    candidates={"ACP1","CPVL","F12","GSTA1","GSTA3","HLA-E","INHBC","SDCCAG8","UMOD"}
    if len(rows)!=9 or {r.get("gene_symbol") for r in rows}!=candidates:raise ValueError("CKD nine-gene panel mismatch")
    if out.exists() and any(out.iterdir()):raise ValueError("Output already populated")
    # Graph labels refer to analysis stages, NOT fabricated transition counts.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch
    fig,ax=plt.subplots(figsize=(12,7))
    ax.set_xlim(0,12);ax.set_ylim(0,7);ax.axis("off")
    stages=[
      (0.3,5.45,"Plasma cis-pQTL","UKB-PPP ST16",2.65),
      (4.65,5.45,"Kidney phenotype GWAS","EUR eGFR / CKD / BUN",2.75),
      (8.95,5.45,"Ancestry support","EAS eGFR / BUN",2.7),
      (0.3,3.3,"Genetic associations","Allele QC + F statistic + MR",3.25),
      (4.65,3.3,"Regional evidence","ABF coloc + SuSiE diagnostics",3.25),
      (8.95,3.3,"Tissue context","Kidney pQTL / eQTL + atlas",2.7),
      (4.25,.95,"Frozen candidates (n = 9)","Independent scientific QA pending",3.65)
    ]
    for x,y,title,sub,w in stages:
      patch=FancyBboxPatch((x,y),w,1.05,boxstyle="round,pad=0.09",fill=False,lw=1.4)
      ax.add_patch(patch);ax.text(x+w/2,y+.68,title,ha="center",va="center",fontweight="bold",fontsize=11)
      ax.text(x+w/2,y+.29,sub,ha="center",va="center",fontsize=8.8)
    def arrow(x1,y1,x2,y2):
      ax.annotate("",(x2,y2),(x1,y1),arrowprops={"arrowstyle":"->","lw":1.5})
    arrow(1.6,5.35,1.6,4.48);arrow(6,5.35,6,4.48);arrow(10.2,5.35,10.2,4.48)
    arrow(2.4,3.2,5.0,2.15);arrow(6.1,3.2,6.1,2.15);arrow(10.3,3.2,7.5,2.15)
    ax.text(6,6.85,"CKD | Multi-omics evidence prioritization workflow",ha="center",fontsize=14,fontweight="bold")
    ax.text(6,.25,"9 candidates = reported integrated analysis panel; NOT a counted flow from 1,936 proteins.",ha="center",fontsize=9)
    out.mkdir(parents=True,exist_ok=True)
    fig.savefig(out/"CKD_FIG1_ANALYSIS_WORKFLOW.svg",bbox_inches="tight");plt.close(fig)
    (out/"CKD_FIG1_PROVENANCE.json").write_text(json.dumps({"ckd_sha256":sha,"reported_panel_n":9,"attrition_counts_asserted":False,"scientific_validation":False},indent=2))
    return {"figure":"CKD_FIG1_ANALYSIS_WORKFLOW.svg","reported_panel_n":9}
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--ckd",type=Path,required=True);p.add_argument("--sha256",required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    print(render(a.ckd,a.out,a.sha256))
