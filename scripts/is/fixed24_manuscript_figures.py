#!/usr/bin/env python3
"""IS manuscript figures from the locked 192-contrast artifact, no source reanalysis."""
import argparse,csv,json,math
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import matplotlib.image as mpimg
import numpy as np
PAIR=[("EAS_"+t,"EUR_"+t) for t in ("AIS","CES","LAS","SVS")]+[(a+"_SVS",a+"_"+t) for a in ("EUR","EAS") for t in ("CES","LAS")]
LABEL=["AIS\nEAS-EUR","CES\nEAS-EUR","LAS\nEAS-EUR","SVS\nEAS-EUR","EUR\nSVS-CES","EUR\nSVS-LAS","EAS\nSVS-CES","EAS\nSVS-LAS"]
def require(ok,msg):
    if not ok:raise ValueError(msg)
def load(artifact,evidence):
    d=json.loads(artifact.read_text())
    e=json.loads(evidence.read_text())
    require(e["provenance"]["run_id"]==45534891070 and e["provenance"]["artifact_id"]==14987743346,"Wrong evidence provenance")
    a=e["result"]
    for k in ("status","family","anchor_verified","matched_cells","matched_contrasts","significant_zero_covariance","significant_max_variance","study_counts"):
        require(d[k]==a[k],"Full artifact aggregate differs: "+k)
    require(d["status"]=="EXPLORATORY_FIXED_24_SUBTYPE_COMPARISONS" and d["family"]==192 and d["anchor_verified"]==24 and d["matched_cells"]==170 and d["matched_contrasts"]==158,"Unexpected frozen analysis counts")
    require(len(d["significant_zero_covariance"])==11 and len(d["significant_max_variance"])==4,"Unexpected corrected count")
    m=d["manifest"]["rows"];ids=sorted(m)
    require(len(ids)==24 and all(x["status"]=="ANCHOR_VERIFIED" for x in m.values()),"Frozen anchor manifest")
    rows={}
    for row in d["contrasts"]:
        rid=row["region_id"];k=(rid,row["first"],row["second"])
        require(rid in m and row["rsid"]==m[rid]["rsid"] and (k[1],k[2]) in PAIR and k not in rows,"Contrast identity")
        if row["status"]=="CONDITIONAL_REPORTED_BETA_CONTRAST":
            delta=row["delta"];lo,hi=row["ci95_max_variance"];p0,pm=row["p_zero_covariance"],row["p_max_variance"]
            require(all(math.isfinite(v) for v in (delta,lo,hi,p0,pm)) and lo<=delta<=hi and 0<=p0<=pm<=1,"Invalid contrast fields")
            require(math.isclose(row["bonferroni_192_zero"],min(1,192*p0),abs_tol=1e-14) and math.isclose(row["bonferroni_192_max"],min(1,192*pm),abs_tol=1e-14),"Fixed family correction")
        else:require(row["status"]=="HOLD_UNMATCHED","Unexpected missing status")
        rows[k]=row
    require(len(rows)==192 and set(rows)=={(rid,a,b) for rid in ids for a,b in PAIR},"Incomplete frozen contrast grid")
    valid=[v for v in rows.values() if v["status"]=="CONDITIONAL_REPORTED_BETA_CONTRAST"]
    require(len(valid)==158 and sum(v["bonferroni_192_zero"]<.05 for v in valid)==11 and sum(v["bonferroni_192_max"]<.05 for v in valid)==4,"Readback counts")
    return d,ids,rows
def save(fig,stem):
    for ext in ("png","pdf","svg"):fig.savefig(stem.with_suffix("."+ext),dpi=350,bbox_inches="tight")
    plt.close(fig)
    p=mpimg.imread(stem.with_suffix(".png"))
    require(p.shape[0]>600 and p.shape[1]>900 and p.std()>.01,"Blank/undersized figure")
    return list(p.shape)
def main():
    p=argparse.ArgumentParser();p.add_argument("--contrasts",type=Path,required=True);p.add_argument("--evidence",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    d,ids,rows=load(a.contrasts,a.evidence);a.out.mkdir(parents=True,exist_ok=True)
    m=d["manifest"]["rows"]
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"svg.fonttype":"none","pdf.fonttype":42,"axes.spines.top":False,"axes.spines.right":False})
    arr=np.full((24,8),np.nan)
    for i,rid in enumerate(ids):
        for j,(first,second) in enumerate(PAIR):
            r=rows[(rid,first,second)]
            if r["status"]=="CONDITIONAL_REPORTED_BETA_CONTRAST":arr[i,j]=r["delta"]
    require(np.isfinite(arr).sum()==158,"Matrix coverage drift")
    vmax=float(np.nanmax(np.abs(arr)));require(vmax>0,"No numerical variation")
    fig,ax=plt.subplots(figsize=(10.7,12.2),layout="constrained")
    cmap=plt.get_cmap("RdBu_r").copy();cmap.set_bad("#dedede")
    im=ax.imshow(np.ma.masked_invalid(arr),norm=TwoSlopeNorm(vmin=-vmax,vcenter=0,vmax=vmax),cmap=cmap,aspect="auto",interpolation="none")
    ax.set_xticks(range(8),LABEL,fontsize=9);ax.xaxis.tick_top();ax.tick_params(axis="x",length=0,pad=7)
    ax.set_yticks(range(24),[rid.replace("IS_XDATA_","")+"  "+m[rid]["rsid"] for rid in ids],fontsize=9)
    ax.tick_params(axis="y",length=0)
    for i,rid in enumerate(ids):
        for j,(first,second) in enumerate(PAIR):
            x=rows[(rid,first,second)]
            if x["status"]!="CONDITIONAL_REPORTED_BETA_CONTRAST":ax.text(j,i,"—",ha="center",va="center",fontsize=8,color="#777777");continue
            if x["bonferroni_192_max"]<.05:ax.text(j,i,"*",ha="center",va="center",fontsize=16,fontweight="bold",color="black" if abs(x["delta"])<.7*vmax else "white")
    ax.axvline(3.5,color="#333333",lw=1.4);ax.axvline(5.5,color="#333333",lw=1)
    ax.set_title("A  Fixed 24-marker, 192-contrast family\n158 evaluable; 34 unresolved",loc="left",fontsize=14,fontweight="bold",pad=52)
    cb=fig.colorbar(im,ax=ax,shrink=.64,pad=.025);cb.set_label("Reported beta difference (first minus second)")
    fig.supxlabel("*  Bonferroni-192 P < 0.05 under maximum variance; gray = unresolved",fontsize=9)
    shape_a=save(fig,a.out/"is_fixed24_all_contrasts")
    sig=sorted(d["significant_max_variance"],key=lambda x:(x["bonferroni_192_max"],x["rsid"]))
    fig,(ax,info)=plt.subplots(1,2,figsize=(11.5,4.5),gridspec_kw={"width_ratios":[3.3,1.4]},layout="constrained")
    for j,x in enumerate(sig):
        y=len(sig)-1-j;b=x["delta"];lo,hi=x["ci95_max_variance"]
        ax.errorbar(b,y,xerr=[[b-lo],[hi-b]],fmt="o",color="#315e78",capsize=3,markersize=6)
        info.text(.04,y,f'{x["bonferroni_192_max"]:.3g}',va="center",fontsize=10)
    ax.set_yticks(list(range(len(sig)-1,-1,-1)),[f'{x["rsid"]} | {x["first"]} - {x["second"]}' for x in sig],fontsize=10)
    ax.set_ylim(-.6,len(sig)-.4);info.set_ylim(ax.get_ylim());info.set_xlim(0,1);info.axis("off")
    ax.axvline(0,color="#888888",ls="--",lw=.8)
    ax.set_xlabel("Reported beta difference (individual 95% CI)")
    info.set_title("Adjusted P\n(maximum variance)",loc="left",fontsize=10)
    fig.suptitle("B  Four contrasts passing the conservative fixed-family screen",fontsize=13,fontweight="bold")
    shape_b=save(fig,a.out/"is_fixed24_conservative_contrasts")
    cols=["region_id","rsid","first","second","status","delta","ci95_max_variance","p_zero_covariance","p_max_variance","bonferroni_192_zero","bonferroni_192_max"]
    with (a.out/"is_fixed24_complete_contrasts.tsv").open("w",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=cols,delimiter="\t",extrasaction="ignore");wr.writeheader()
        for rid in ids:
            for a1,a2 in PAIR:
                row=dict(rows[(rid,a1,a2)])
                if "ci95_max_variance" in row:row["ci95_max_variance"]=json.dumps(row["ci95_max_variance"])
                wr.writerow(row)
    (a.out/"LEGENDS.md").write_text("""# Draft IS figure legends

A. All fixed 192 tests in the 24 historically EAS-selected regional marker by 8 prespecified contrast grid. 158 evaluable, 34 unresolved (gray, never zero). Reported signed beta differences use the first phenotype/ancestry minus the second. Asterisks denote the four contrasts with Bonferroni-192 P < 0.05 under the maximum-variance covariance bound. Color range is symmetric over all available effects; no effect-based exclusion. Markers were not LD-pruned or independently selected. B. The four qualifying contrasts (three markers), with individual 95% intervals under the maximum-variance assumption and multiplicity-adjusted P values. Intervals are not simultaneous or multiplicity-corrected. Conditional comparisons assume sufficiently comparable reported effect scales, not certified case-case, interaction or causal estimates. rs2129981 is a known PITX2 positive-control region; rs2767060 and rs7203740 are post-readout exploratory candidates, not novel confirmed loci. Same GIGASTROKE release; no independent replication. Source run 45534891070, artifact 14987743346, verified against committed aggregate stdout. Technical checks only; visual QA pending.
""")
    print(json.dumps({"status":"FIXED_FAMILY_FIGURES_GENERATED","matched":158,"unresolved":34,"sig_zero":11,"sig_max":4,"matrix_shape":[24,8],"png_shapes":[shape_a,shape_b],"files":[{"name":p.name,"bytes":p.stat().st_size} for p in sorted(a.out.iterdir())],"visual_inspection":False},indent=2,allow_nan=False))
if __name__=="__main__":main()
