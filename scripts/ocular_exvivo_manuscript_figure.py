#!/usr/bin/env python3
"""Descriptive tissue figure from verified aggregate stdout only."""
import argparse,json,math,statistics
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
EVID=ROOT/"docs/publication/evidence"
PIN="5ad39d6fb293130492e66c713ff7d155d03d752f31cac84c98f78c42671bd606"
def require(ok,msg):
    if not ok:raise ValueError(msg)
def read(name,run,status):
    x=json.loads((EVID/name).read_text())
    require(x["provenance"]["run_id"]==run,"Wrong source run")
    d=x["result"]
    require(d["status"]==status and d["panel_sha256"]==PIN and d["study"]=="GSE37474" and d["pairs"]==5,"Source identity")
    return d
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    primary=read("ocular_exvivo_paired_20261125.json",45535697612,"PAIRED_TISSUE_DESCRIPTIVE")
    sens=read("ocular_exvivo_gene_rank_sensitivity_20261125.json",45536989022,"POST_READOUT_GENE_UNIVERSE_SENSITIVITY")
    require(primary["pairing_verified_by_public_titles_treatment_age_sex"] and sens["primary_run"]==45535697612,"Pairing provenance")
    require(primary["matrix_sha256"]==sens["matrix_sha256"]=="d348cb1add4a01ec1b17d82da75f4b2081de4d7c3ad2433a42261fdd9543c849" and primary["annotation_sha256"]==sens["annotation_sha256"]=="d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394","Pinned public inputs")
    require(primary["rank_universe"]==sens["eligible_probes"]==27920 and sens["unique_entrez_genes"]==17686,"Rank universe drift")
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,"axes.spines.right":False,"svg.fonttype":"none","pdf.fonttype":42})
    fig,axs=plt.subplots(1,2,figsize=(9.2,4.4),layout="constrained",sharey=True)
    colors=["#245c78","#b46c36"]
    for j,(k,n,expected) in enumerate((("all241",241,199),("ECM16",16,14))):
        d=primary["panels"][k];s=sens["panels"][k]
        require(d["frozen"]==n and d["mapped"]==s["mapped"]==expected,"Fixed panel denominator")
        require(d["positive_pairs"]==s["positive_pairs"]==5 and math.isclose(d["two_sided_exact_sign_p"],.0625) and math.isclose(d["holm_p"],.125),"Primary statistic changed")
        ax=axs[j]
        ax.axhline(0,color=".65",lw=.9,ls="--")
        for pos,v,col in [(0,d["pair_score_distribution_sorted"],colors[0]),(1,s["pair_scores_sorted"],colors[1])]:
            require(len(v)==5 and all(math.isfinite(y) and y>0 for y in v),"Pair distribution")
            obj=d if pos==0 else s
            require(math.isclose(statistics.median(v),obj["median_signed_percentile_shift"]),"Median mismatch")
            ax.scatter(pos+np.array([-.12,.08,-.04,.12,0]),np.array(v)*100,color=col,s=42,alpha=.9,zorder=3)
            med=statistics.median(v)*100
            ax.plot([pos-.23,pos+.23],[med,med],color=col,lw=2.5,zorder=2)
            ax.text(pos,14.1,f"5/5 positive\nMedian {med:.2f}",ha="center",va="top",fontsize=9)
        ax.set_xticks([0,1],["Primary\nprobe ranks","Sensitivity\nunique-gene ranks"])
        ax.set_xlim(-.55,1.55);ax.set_ylim(-.7,15)
        ax.set_title(f"{'A' if j==0 else 'B'}  {n}-gene panel\n{expected}/{n} mapped",loc="left",fontweight="bold")
    axs[0].set_ylabel("Signed percentile shift (percentage points)")
    fig.suptitle("GSE37474: paired ex-vivo tissue response",fontsize=14,fontweight="bold")
    a.out.mkdir(parents=True,exist_ok=True)
    stem=a.out/"fig3_exvivo_tissue_sensitivity"
    for ext in ("svg","png","pdf"):fig.savefig(stem.with_suffix("."+ext),dpi=350,bbox_inches="tight")
    plt.close(fig)
    legend="""# Draft Figure 3 legend

Fixed glucocorticoid-responsive transcriptional signatures in ten arrays from five paired donor eyes (GSE37474; mixed trabecular meshwork and corneoscleral tissue; 100 nM dexamethasone, 10 days). A, 241-gene signature (199 mapped); B, 16-gene ECM organization subset (14 mapped). Each dot is one paired donor score: mean of the discovery-direction-signed changes in within-array percentile ranks (percentage points). Horizontal lines denote medians, not confidence intervals. Primary: rank 27,920 eligible, complete nonshared probe sets, then median mapped probe ranks by candidate gene. Post-readout sensitivity: median expression across eligible probe sets per Entrez ID, then rank 17,686 unique Entrez genes. Both panels had five positive donor scores in each analysis. Primary two-sided exact sign P = 0.0625 for each panel, Holm P = 0.125 for each. The sensitivity reuses exactly the same donor pairs; no additional confirmatory test. Dots are independently sorted/jittered within groups: no donor linkage can be inferred from positions. Previously analyzed public dataset; original donor independence and selection lineage not fully certified. No clinical phenotype or novel ECM mechanism inferred. Source: verified aggregate stdout runs 45535697612 and 45536989022; reported artifact ZIP digests archived alongside JSON. No raw expression in this figure workflow. Draft; visual review and journal layout pending.
"""
    (a.out/"FIGURE3_LEGEND.md").write_text(legend)
    im=mpimg.imread(stem.with_suffix(".png"))
    require(im.shape[0]>=800 and im.shape[1]>=1200 and float(im.std())>.01,"Blank or undersized PNG")
    for ext in ("svg","png","pdf"):require(stem.with_suffix("."+ext).stat().st_size>1000,"Missing/empty figure")
    print(json.dumps({"status":"DESCRIPTIVE_FIGURE_GENERATED","png_shape":list(im.shape),"files":[{"name":v.name,"bytes":v.stat().st_size} for v in sorted(a.out.iterdir())],"visual_inspection":False},indent=2))
if __name__=="__main__":main()
