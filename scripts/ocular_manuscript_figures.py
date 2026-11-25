#!/usr/bin/env python3
"""Draw descriptive manuscript figures from archived, verified aggregate stdout."""
import argparse, json, math
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import PercentFormatter

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"docs/publication/evidence/ocular_manuscript_aggregate_20261125.json"
PIN="5ad39d6fb293130492e66c713ff7d155d03d752f31cac84c98f78c42671bd606"
def require(ok,msg):
    if not ok:raise ValueError(msg)
def read(p):
    x=json.loads(p.read_text())
    require(x["schema"]==1 and x["frozen_panel_sha256"]==PIN,"Wrong aggregate source")
    src=x["sources"]
    require({k:v["run_id"] for k,v in src.items()}=={"seven_day":45383400882,"paired":45386503377,"crossdrug":45530972572},"Wrong run provenance")
    a,b,c=[src[k]["payload"] for k in ("seven_day","paired","crossdrug")]
    require(a["frozen_panel_sha256"]==c["panel_sha256"]==PIN,"Panel identity mismatch")
    require(c["biological_strains"]==2 and c["arrays"]==23,"Crossdrug design")
    require(x["paired_declared_design"]["pairs"]==8 and x["paired_declared_design"]["complete_pairs"]==5,"Paired design")
    for k,n in (("all_241",241),("ECM_16",16)):
        d=a["groups"][k];z=d["strata"][0]
        require(d["frozen_total"]==d["mapped"]==z["evaluable"]==n,"Seven-day denominator")
        require(0<=z["all_three_strains_same_sign_abs_ge_025"]<=z["all_three_strains_same_sign"]<=z["aggregate_same_sign"]<=n,"Seven-day direction counts")
        q=b[k]
        require(q["pairs"]==5 and q["positive"]+q["negative"]+q["ties"]==5,"Pair partition")
        require(q["positive"]==5 and math.isclose(q["two_sided_p"],.0625) and math.isclose(q["holm_p"],.125),"Paired original test changed")
    for k,n in (("all241",241),("ECM16",16)):
        d=c["panels"][k];m=d["mapped"]
        require(d["frozen"]==n and m<=n and m>=math.ceil(.75*n),"Crossdrug mapping")
        require(set(d["scores"])=={f"{s}_{a}" for s in ("86","93") for a in ("Dex","FA","TA")},"Crossdrug design")
        require(all(0<=v["positive_genes"]<=m and math.isfinite(v["mean_signed_percentile_shift"]) for v in d["scores"].values()),"Crossdrug effects")
    require(len(a["ecm_gene_results"])==16 and len({r["symbol"] for r in a["ecm_gene_results"]})==16,"ECM identities")
    return x,a,b,c
def save(fig,out,stem):
    for ext in ("svg","png","pdf"):
        fig.savefig(out/f"{stem}.{ext}",dpi=350,bbox_inches="tight")
    plt.close(fig)
def figures(src,out):
    x,a,b,c=read(src)
    out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,"axes.spines.right":False,"svg.fonttype":"none","pdf.fonttype":42})
    colors=["#245c78","#b46c36"]
    fig,axs=plt.subplots(1,3,figsize=(12,4.4),layout="constrained",gridspec_kw={"width_ratios":[1.25,1.1,1.0]})
    names=["Mean direction\nagrees","All 3 strains\nagree","All 3 agree;\n|effect| ≥ 0.25"]
    keys=["aggregate_same_sign","all_three_strains_same_sign","all_three_strains_same_sign_abs_ge_025"]
    for k,(g,n) in enumerate((("all_241",241),("ECM_16",16))):
        z=a["groups"][g]["strata"][0];v=[z[t] for t in keys];pos=np.arange(3)+(k-.5)*.34
        axs[0].barh(pos,np.array(v)/n,height=.31,color=colors[k],label=f"{n}-gene panel")
        for p,j in zip(pos,v):axs[0].text(j/n+.015,p,f"{j}/{n}",va="center",fontsize=8)
    axs[0].set_yticks(np.arange(3),names);axs[0].invert_yaxis()
    axs[0].set_xlim(0,1.24);axs[0].xaxis.set_major_formatter(PercentFormatter(1));axs[0].set_xticks([0,.5,1])
    axs[0].set_title("A  Seven-day DEX\nGSE301517: 3 strains",loc="left",fontweight="bold")
    axs[0].legend(frameon=False,fontsize=8,loc="lower left",bbox_to_anchor=(0,-.25))
    order=[f"{s}_{g}" for s in ("86","93") for g in ("Dex","FA","TA")]
    for k,g in enumerate(("all241","ECM16")):
        v=[c["panels"][g]["scores"][t]["mean_signed_percentile_shift"] for t in order]
        axs[1].bar(np.arange(6)+(k-.5)*.34,v,width=.31,color=colors[k])
    axs[1].set_xticks(range(6),["DEX","FA","TA"]*2,rotation=35)
    axs[1].axvline(2.5,color="0.7",lw=.8);axs[1].text(1,.107,"TM86",ha="center");axs[1].text(4,.107,"TM93",ha="center")
    axs[1].set_ylim(0,.12);axs[1].set_ylabel("Mean signed percentile shift")
    axs[1].set_title("B  Cross-drug response\nGSE16643: 2 strains",loc="left",fontweight="bold")
    vals=[b[g]["median_signed_log2_difference"] for g in ("all_241","ECM_16")]
    axs[2].bar([0,1],vals,color=colors,width=.55)
    for i,v in enumerate(vals):axs[2].text(i,v+.04,f"{v:.3f}",ha="center")
    axs[2].set_xticks([0,1],["241 genes","ECM16"]);axs[2].set_ylim(0,1.25)
    axs[2].set_ylabel("Median signed log2(CPM+1) difference")
    axs[2].set_title("C  Paired seven-day DEX\nPRJNA729873: 5/8 pairs",loc="left",fontweight="bold")
    axs[2].text(.5,1.12,"5/5 positive each\nHolm P = 0.125",ha="center",fontsize=9)
    fig.suptitle("Fixed glucocorticoid-response signature: descriptive transportability",fontsize=14,fontweight="bold")
    save(fig,out,"fig2_signature_transportability")
    fig,ax=plt.subplots(figsize=(7.4,7.1),layout="constrained")
    e=sorted(a["ecm_gene_results"],key=lambda v:v["acute_log2fc"])
    y=np.arange(16);ax.axvline(0,color=".6",lw=.9)
    for i,v in enumerate(e):
        ax.plot([v["acute_log2fc"],v["mean_7d_log2_cpm_ratio"]],[i,i],color=".75",lw=1,zorder=0)
    ax.scatter([v["acute_log2fc"] for v in e],y,color=colors[0],s=34,label="Discovery pseudobulk log2FC")
    ax.scatter([v["mean_7d_log2_cpm_ratio"] for v in e],y,color=colors[1],marker="D",s=31,label="Seven-day mean log2 CPM ratio")
    ax.set_yticks(y,[v["symbol"] for v in e]);ax.set_xlabel("Reported expression effect (study-specific definition)")
    ax.set_title("ECM16: fixed gene-level direction across studies",loc="left",fontweight="bold")
    ax.legend(frameon=False,loc="lower center",bbox_to_anchor=(.5,-.13),ncol=1)
    save(fig,out,"figS1_ecm16_direction")
    legends="""# Draft figure legends — descriptive evidence only

Figure 2. Fixed 241-gene and 16-gene extracellular-matrix signatures across public trabecular meshwork expression studies. A: Directional concordance with discovery effects in three GSE301517 strains after seven days of dexamethasone (DEX), including a descriptive absolute per-strain log2(CPM ratio) threshold of 0.25. This previously examined study contributed to candidate selection at 48 h; seven-day analysis is not a blind holdout. B: Mean signed percentile rank shifts after DEX, fluocinolone acetonide (FA), or triamcinolone acetonide (TA) in GSE16643 (two biological strains; 23 arrays). 229/241 and 16/16 candidates mapped. Rank universe consists of 30,548 complete, unique-Entrez annotated probes, not the original all-noncontrol universe. C: Paired scores in PRJNA729873: five of eight declared donor-derived pairs have both arms in the pinned DEE2 release (four responders and one nonresponder; three nonresponder pairs missing one arm). Median signed log2(CPM+1) differences shown; each panel has 5/5 positive and two-sided exact sign P=.0625, Holm P=.125. The bar heights are not a shared effect scale and no error bars are invented from aggregate-only logs. No gene-independence inference, no clinical prediction, no drug equivalence or confirmatory independent cohort claim.

Figure S1. Discovery donor pseudobulk log2FC and GSE301517 equal-strain mean seven-day log2 CPM ratio for all 16 frozen ECM genes, ordered by discovery effect. These different estimands are shown side by side to describe direction, not as a formal effect-size comparison or independent replication. No confidence intervals available from the plotted summary; no claim of causality or generalized ECM activation.

Provenance: archived aggregate stdout readbacks in docs/publication/evidence/ocular_manuscript_aggregate_20261125.json; original Actions runs 45383400882, 45386503377, 45530972572. Drafts pending independent donor validation, review and figure QC. No raw donor IDs, raw expression or new analysis.
"""
    (out/"FIGURE_LEGENDS.md").write_text(legends)
def main():
    p=argparse.ArgumentParser();p.add_argument("--source",type=Path,default=SOURCE);p.add_argument("--out",type=Path,default=Path("outputs/ocular_figures"));a=p.parse_args();figures(a.source,a.out);print("Generated 2 descriptive figures in SVG/PNG/PDF and draft legends; no new inferential result.")
if __name__=="__main__":main()
