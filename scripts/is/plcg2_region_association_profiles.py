#!/usr/bin/env python3
"""Post-readout descriptive regional association profiles only; no coloc or MR."""
import argparse,json,math,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

STUDIES=("EAS_CES","EAS_SVS","EUR_CES","EUR_SVS")
TARGET="chr16_81775211_C_T"
LEAD="chr16_81878856_T_G"
EXPECTED={"EAS_CES":2230,"EAS_SVS":2398,"EUR_CES":2963,"EUR_SVS":2612}
def require(ok,msg):
    if not ok:raise ValueError(msg)
def score(p):
    require(math.isfinite(p) and 0<p<=1,"Zero/invalid P: do not fabricate log P")
    return -math.log10(p)
def ranks(x):
    """Average tied ranks; descriptive correlation, no inferential P."""
    x=np.asarray(x,float);order=np.argsort(x,kind="stable");r=np.empty(len(x),float)
    i=0
    while i<len(x):
        j=i+1
        while j<len(x) and x[order[j]]==x[order[i]]:j+=1
        r[order[i:j]]=(i+j-1)/2+1;i=j
    return r
def corr(x,y):
    x,y=ranks(x),ranks(y)
    if np.std(x)==0 or np.std(y)==0:return None
    return float(np.corrcoef(x,y)[0,1])
def top_ties(data,key,k):
    threshold=sorted(d[key] for d in data)[min(k,len(data))-1]
    return {d["variant"] for d in data if d[key]<=threshold}
def summarize(data):
    require(len(data)>0,"Empty matched set")
    q=[d["qtl_p"] for d in data];g=[d["gwas_p"] for d in data]
    qz=[d["qtl_beta_ALT"]/d["qtl_se"] for d in data];gz=[d["gwas_beta_ALT"]/d["gwas_se"] for d in data]
    result={"n":len(data),"spearman_minus_log10_p":corr([score(p) for p in q],[score(p) for p in g]),"spearman_signed_z":corr(qz,gz),"same_sign_count":sum(a*b>0 for a,b in zip(qz,gz)),"zero_effect_count":sum(a*b==0 for a,b in zip(qz,gz)),"top_qtl":min(data,key=lambda d:(d["qtl_p"],d["variant"])),"top_gwas":min(data,key=lambda d:(d["gwas_p"],d["variant"])),"gwas_p_lt_5e_8":sum(p<5e-8 for p in g),"top_overlap":{},"anchors":{}}
    for k in (10,50,100):
        a=top_ties(data,"qtl_p",k);b=top_ties(data,"gwas_p",k)
        result["top_overlap"][str(k)]={"qtl_including_boundary_ties":len(a),"gwas_including_boundary_ties":len(b),"intersection":len(a&b),"jaccard":len(a&b)/len(a|b)}
    for label,var in (("target",TARGET),("lead",LEAD)):
        rows=[d for d in data if d["variant"]==var]
        if not rows:result["anchors"][label]={"status":"NOT_IN_MATCHED_SET"};continue
        require(len(rows)==1,"Duplicate anchor")
        d=rows[0]
        result["anchors"][label]={"status":"MATCHED","variant":var,"rsid":d["rsid"],"qtl_p":d["qtl_p"],"gwas_p":d["gwas_p"],"qtl_rank_min":1+sum(p<d["qtl_p"] for p in q),"qtl_rank_max":sum(p<=d["qtl_p"] for p in q),"gwas_rank_min":1+sum(p<d["gwas_p"] for p in g),"gwas_rank_max":sum(p<=d["gwas_p"] for p in g),"qtl_alt_eaf":d["qtl_alt_eaf"],"gwas_alt_eaf":d["gwas_alt_eaf"]}
    return result
def load(path,evidence):
    full=json.loads(path.read_text())
    e=json.loads(evidence.read_text())
    require(e["provenance"]["run_id"]==45560562264 and e["provenance"]["artifact_id"]==14997647839,"Pinned artifact")
    d=e["result"]
    require(full["status"]==d["status"]=="POST_READOUT_REGION_COVERAGE_GATE" and full["formal_coloc_status"]=="HOLD" and full["comparisons_completed"]==4,"Descriptive source only")
    for k in ("GRCh38_region","gene","target","official_eGene_lead","all_four_common_variants","all_four_common_qtl_coverage","target_in_all_four","lead_in_all_four"):
        require(full[k]==d[k],"Aggregate parity "+k)
    require(full["target"]==TARGET and full["official_eGene_lead"]==LEAD and full["qtl"]["strict_nonpal_snp_count"]==3555 and len(full["qtl"]["records"])==3555,"Fixed nominal universe")
    require({k:v for k,v in full["qtl"].items() if k!="records"}==d["qtl"],"QTL aggregate parity")
    universe=full["qtl"]["records"]
    by={}
    for s in STUDIES:
        v=full["comparisons"][s]
        require({k:x for k,x in v.items() if k!="matched_rows"}==d["comparisons"][s],"GWAS aggregate parity "+s)
        require(v["status"]=="DESCRIPTIVE_ALLELE_MATCHED_REGION" and v["matched"]==EXPECTED[s] and len(v["matched_rows"])==EXPECTED[s],"Matched row count")
        rows={}
        for z in v["matched_rows"]:
            k=z["variant"];require(k not in rows and k in universe,"Matched key uniqueness")
            a=universe[k]
            require(z["rsid"] in a["rsids"] and all(z[field]==a[field] for field in ("pos","ref","alt")),"Exact variant identity")
            require(all(z[f]==a[t] for f,t in (("qtl_beta_ALT","beta_ALT"),("qtl_se","se"),("qtl_p","p"),("qtl_alt_eaf","alt_eaf"))),"Nominal QTL value drift")
            require(81275211<=z["pos"]<=82275211 and z["qtl_se"]>0 and z["gwas_se"]>0 and all(math.isfinite(z[f]) for f in ("qtl_beta_ALT","qtl_se","gwas_beta_ALT","gwas_se")) and 0<z["qtl_alt_eaf"]<1 and 0<z["gwas_alt_eaf"]<1,"Effect/AF checks")
            score(z["qtl_p"]);score(z["gwas_p"])
            rows[k]=z
        by[s]=rows
    common=set.intersection(*(set(by[s]) for s in STUDIES))
    require(len(common)==1966 and TARGET not in common and LEAD in common,"Fixed intersection")
    return full,by,common
def main():
    p=argparse.ArgumentParser();p.add_argument("--input",type=Path,required=True);p.add_argument("--evidence",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    d,by,common=load(a.input,a.evidence)
    a.out.mkdir(parents=True,exist_ok=True)
    result={"status":"DESCRIPTIVE_ASSOCIATION_PROFILES_NOT_COLOC","source_run":45560562264,"common_count":len(common),"target_in_common":TARGET in common,"lead_in_common":LEAD in common,"comparisons":{},"limits":["Post-readout descriptive, fixed +/-500 kb region. Neither exact matched set nor four-study intersection certifies the full cis universe.","Rank correlations/top-k overlaps are descriptive only; LD induces dependence. No correlation significance tests, Bayes factors, PP.H4, MR or causal inference. The common intersection excludes rs7203740 and preferentially retains cross-population variants.","The nominal and permutation lead estimates remain nonequivalent; use nominal only. ROSMAP ancestry, model N/sdY, GWAS model equivalence and matched LD remain HOLD. No independent replication."]}
    for s in STUDIES:
        rows=list(by[s].values());shared=[by[s][k] for k in sorted(common)]
        result["comparisons"][s]={"pairwise_matched":summarize(rows),"four_study_common":summarize(shared)}
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"svg.fonttype":"none","pdf.fonttype":42,"axes.spines.top":False,"axes.spines.right":False})
    for mode in ("pairwise_matched","four_study_common"):
        fig,axes=plt.subplots(2,2,figsize=(11,8.5),layout="constrained")
        for ax,s in zip(axes.flat,STUDIES):
            rows=list(by[s].values()) if mode=="pairwise_matched" else [by[s][k] for k in sorted(common)]
            ax.scatter([score(d["qtl_p"]) for d in rows],[score(d["gwas_p"]) for d in rows],s=8,alpha=.28,c="#58788a",edgecolors="none",rasterized=True)
            for key,color,marker in [(TARGET,"#b46c36","D"),(LEAD,"#237c6b","s")]:
                if key in by[s] and (mode=="pairwise_matched" or key in common):
                    v=by[s][key]
                    ax.scatter(score(v["qtl_p"]),score(v["gwas_p"]),s=55,c=color,marker=marker,edgecolors="black",linewidths=.5,zorder=4,label=v["rsid"])
            summary=result["comparisons"][s][mode]
            ax.set_title(f"{s.replace('_',' ')} | n={len(rows):,}",loc="left",fontweight="bold")
            rho=summary["spearman_minus_log10_p"]
            ax.text(.03,.96,"Spearman rho = "+(f"{rho:.3f}" if rho is not None else "NA"),transform=ax.transAxes,va="top",fontsize=9)
            ax.set_xlabel("ROSMAP PLCG2 nominal -log10(P)")
            ax.set_ylabel("GWAS -log10(P)")
            ax.set_xlim(left=-.25);ax.set_ylim(bottom=-.15)
            ax.legend(loc="upper right",fontsize=8,frameon=False)
        fig.suptitle("PLCG2 regional association profiles — "+("pairwise allele matches" if mode=="pairwise_matched" else "four-study intersection (index SNP absent)"),fontsize=13,fontweight="bold")
        stem=a.out/("plcg2_"+mode)
        for ext in ("png","pdf","svg"):fig.savefig(stem.with_suffix("."+ext),dpi=300,bbox_inches="tight")
        plt.close(fig)
        im=mpimg.imread(stem.with_suffix(".png"));require(im.shape[0]>800 and im.shape[1]>1200 and float(im.std())>.01,"Blank image")
    (a.out/"LEGENDS.md").write_text("# Draft supplementary figure legends\n\nEach dot is an exact allele- and rsID-matched nominal QTL/GWAS variant in the fixed 1 Mb interval. Diamond: rs7203740, square: rs11640696. Panels show pairwise matched sets or the common four-GWAS intersection; the latter excludes the index variant. Spearman coefficients describe patterns, not independent observations, statistical significance or colocalization. Nominal QTL estimates only; no claim of ancestry-specific causal mechanism. Visual QA pending.\n")
    result["files"]=[{"name":p.name,"bytes":p.stat().st_size} for p in sorted(a.out.iterdir())]
    result["visual_inspection"]=False
    (a.out/"association_profile_summary.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps(result,indent=2,allow_nan=False),flush=True)
if __name__=="__main__":main()
