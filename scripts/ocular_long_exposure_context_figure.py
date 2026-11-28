#!/usr/bin/env python3
"""Study-specific descriptive score panels, without pooling or a time trajectory."""
import argparse,json,math,statistics
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
PIN="5ad39d6fb293130492e66c713ff7d155d03d752f31cac84c98f78c42671bd606"
def req(ok,msg):
    if not ok:raise ValueError(msg)
def read(root,path,run,artifact,status):
    v=json.loads((root/path).read_text())
    req(v["provenance"]["run_id"]==run and v["provenance"]["artifact_id"]==artifact,"Evidence identity")
    d=v["result"];req(d["status"]==status,"Evidence status");return d
def main():
    p=argparse.ArgumentParser();p.add_argument("--evidence",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    d=read(a.evidence,"ocular_exvivo_paired_20261125.json",45535697612,14987116725,"PAIRED_TISSUE_DESCRIPTIVE")
    v=read(a.evidence,"ocular_14d_conditional_score_20261125.json",45563105992,14998946580,"CONDITIONAL_DEPOSITOR_LABEL_DESCRIPTIVE")
    g=read(a.evidence,"ocular_7d_frozen_directional_readback_20261125.json",45383400882,14916031829,"FROZEN_PANEL_DIRECTIONAL_AUDIT_COMPLETE")
    req(d["panel_sha256"]==v["panel_sha256"]==g["frozen_panel_sha256"]==PIN and v["primary_source_gate"]=="HOLD_SOURCE_NAME_CONFLICT" and v["source_conflicts"]==4 and not v["donor_independence_certified"] and d["pairs"]==5 and v["commercial_lots"]==3,"Panel/claim gate")
    for name,k1,k2,n1,n2 in [("all241","all_241","all241",199,220),("ECM16","ECM_16","ECM16",14,16)]:
        req(d["panels"][name]["mapped"]==n1 and v["panels"][k2]["mapped"]==n2,"Mapping")
        req(d["panels"][name]["positive_pairs"]==5 and v["panels"][k2]["positive_lots"]==3 and v["panels"][k2]["inferential_p"] is None,"No new inference")
        req(g["groups"][k1]["mapped"]==g["groups"][k1]["frozen_total"],"Seven-day source mapping")
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"svg.fonttype":"none","pdf.fonttype":42,"axes.spines.top":False,"axes.spines.right":False})
    a.out.mkdir(parents=True,exist_ok=True)
    fig,axes=plt.subplots(1,2,figsize=(10.6,4.6),layout="constrained",sharey=True)
    for ax,name,title in zip(axes,("all241","ECM16"),("241-gene panel","ECM16 subset")):
        x=d["panels"][name];y=v["panels"][name]
        q=[100*z for z in x["pair_score_distribution_sorted"]]
        req(len(q)==5 and math.isclose(statistics.median(q),100*x["median_signed_percentile_shift"],abs_tol=1e-9),"10d scores")
        rr=y["all_three_lots"];req(set(rr)=={"2584","3423","4973"} and rr["4973"]["source_name_conflict"] is False and all(rr[k]["source_name_conflict"] is True for k in ("2584","3423")),"14d lot labels")
        qq=[100*rr[k]["signed_percentile_shift"] for k in ("2584","3423","4973")]
        req(math.isclose(statistics.median(qq),100*y["median_conditional_shift"],abs_tol=1e-9),"14d scores")
        ax.scatter([-.16,-.08,0,.08,.16],q,c="#315e78",s=42,label="Paired tissue (10d)",zorder=3)
        ax.scatter([.91,1.0],qq[:2],facecolors="none",edgecolors="#a56a3c",linewidths=1.8,s=62,label="Source-name conflict (14d)",zorder=3)
        ax.scatter([1.09],qq[2:],c="#a56a3c",s=55,marker="D",label="Conflict-free lot (14d)",zorder=3)
        for j,vals in enumerate((q,qq)):ax.hlines(statistics.median(vals),j-.24,j+.24,lw=2,color="#222222",zorder=2)
        ax.axhline(0,color="#888888",lw=.8,ls="--")
        ax.set_xticks([0,1],[f"GSE37474\n10d | n=5 pairs\n{x['mapped']}/{x['frozen']} mapped",f"GSE65240\n14d | 3 lots, conditional\n{y['mapped']}/{y['frozen']} mapped"],fontsize=9)
        ax.set_xlim(-.45,1.45);ax.set_ylim(-.8,15)
        ax.set_title(title,loc="left",fontweight="bold")
        ax.set_ylabel("Direction-aligned percentile shift (pp)")
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="outside upper center",ncol=3,frameon=False,fontsize=9)
    fig.supxlabel("Separate studies and mapping universes; no longitudinal effect estimate or pooled inference.",fontsize=9)
    stem=a.out/"ocular_10d_14d_contexts"
    for ext in ("png","svg","pdf"):fig.savefig(stem.with_suffix("."+ext),dpi=350,bbox_inches="tight")
    plt.close(fig)
    im=mpimg.imread(stem.with_suffix(".png"));req(im.shape[0]>800 and im.shape[1]>1500 and float(im.std())>.01,"Image QC")
    (a.out/"LEGEND.md").write_text("Each point represents one reported paired tissue sample (10d) or commercial lot (14d), not a gene. Bars are separate study medians, not intervals. Open symbols mark the four source-name-conflicting samples represented by two lot pairs; diamond marks lot 4973 alone. Sample identity conflicts and donor independence remain unresolved. No cross-study pooling, trajectory, statistical comparison, or clinical interpretation. Visual QA pending.\n")
    print(json.dumps({"status":"CONTEXT_FIGURE_GENERATED","ten_day_pairs":5,"fourteen_day_lots":3,"primary_source_gate":"HOLD_SOURCE_NAME_CONFLICT","png_shape":list(im.shape),"files":[{"name":p.name,"bytes":p.stat().st_size} for p in sorted(a.out.iterdir())],"visual_inspection":False},indent=2))
if __name__=="__main__":main()
