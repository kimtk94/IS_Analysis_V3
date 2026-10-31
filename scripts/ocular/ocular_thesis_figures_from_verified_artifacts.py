#!/usr/bin/env python3
"""Reproducible ocular thesis figures from *verified aggregate* GitHub artifacts.

NEVER read or upload raw donor pseudobulk counts / private genetic records.
Study 3's CPM+1 log2 ratio is NOT an edgeR log2FC with matching estimand.
No causal drug, time, protein or clinical steroid-delta-IOP inference.
"""
import argparse,csv,hashlib,json,math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

SOURCE_FILES={
 "three_core":"GSE301517_core4_3study_direction.tsv",
 "three_emt":"GSE301517_EMT_194_3study_signed_direction.tsv",
 "third_gate":"GSE301517_3strain_vs_48h_evidence_gate.json",
 "cutoff_panel":"GATE_B_16_candidate_score050_vs075_edgeR.tsv",
 "cutoff_gate":"GATE_B_score050_075_actual_refit_concordance.json",
 "dual_panel":"GATE_B_55_matched_BeamA_vs_scPred_16_genes.tsv",
 "dual_gate":"GATE_B_55_matched_annotation_concordance.json",
}
PIN={
 "three_core":"926b7899639d9b759d03688b911b89becf6b776103927a226e4f301415a86180",
 "three_emt":"edac2a23c47e40fde419214c8c45fed6f1bd33861a72f8c37c789db2f7d88e1d",
 "cutoff_panel":"ffaf1f6376f8c1fb7e041b9890a3723e3396da3a1f365438556f07eea6794e1e",
 "dual_panel":"6b42e89ef03d8284e6875310a4a33fa193abe37fbfd50055eefb77f9d332777e",
}
CORE=("MYOC","TGFBI","MMP2","ACTA2")
NAVY="#273b57";GREEN="#207a72";BLUE="#3a6899";AMBER="#b46a32";GREY="#687586";LT="#eef1f4"
def table(path):
    with path.open("r",encoding="utf8",newline="") as f:
        out=list(csv.DictReader(f,delimiter="\t"))
    if len(out)==0:raise ValueError("EMPTY_TSV_SOURCE "+str(path))
    return out
def num(v):
    x=float(v)
    if not math.isfinite(x):raise ValueError("NONFINITE_VALUE")
    return x
def save(fig,dest,filename):
    fig.savefig(dest/(filename+".png"),dpi=240,facecolor="white",bbox_inches="tight")
    fig.savefig(dest/(filename+".svg"),facecolor="white",bbox_inches="tight")
    plt.close(fig)
def style():
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,
        "axes.titlesize":14,"axes.titleweight":"bold","axes.spines.top":False,
        "axes.spines.right":False,"savefig.pad_inches":0.18,"svg.fonttype":"none"})
def sign_color(x):return BLUE if x<0 else AMBER if x>0 else GREY
def study_heatmap(source,dest):
    c={r["gene_symbol"]:r for r in source}
    if set(c)!=set(CORE):raise ValueError("UNEXPECTED_3STUDY_CORE")
    values=[]
    for g in CORE:
        row=c[g]
        values.append([num(row[x]) for x in ("acute_500nM_48h_log2FC",
          "published_100nM_48h_log2FC","chronic_500nM_7d_CPM_log2ratio")])
    fig,ax=plt.subplots(figsize=(11.6,6.0))
    fig.subplots_adjust(left=.13,right=.98,top=.66,bottom=.29)
    for j in range(3):
        for i,g in enumerate(CORE):
            v=values[i][j]
            ax.add_patch(Rectangle((j-0.47,i-0.43),0.94,0.86,
                       facecolor=sign_color(v),alpha=0.14,edgecolor="#ccd5df"))
            ax.text(j,i-0.08,("+" if v>0 else "−"),ha="center",va="center",
                    fontsize=20,fontweight="bold",color=sign_color(v))
            ax.text(j,i+0.25,f"{v:+.3f}",ha="center",va="center",
                    fontsize=11,color=NAVY)
    ax.set_xlim(-.52,2.52);ax.set_ylim(3.56,-.56)
    ax.set_yticks(range(4),CORE,fontweight="bold")
    ax.set_xticks(range(3),["GSE313736\n500 nM · 48 h\n57 paired donors",
       "Mehrotra et al.\n100 nM · 48 h\npublished bulk TM",
       "GSE301517\n500 nM · 7 d\n3 TM strains"])
    ax.xaxis.tick_top();ax.tick_params(axis="both",length=0,pad=8)
    for s in ax.spines.values():s.set_visible(False)
    fig.suptitle("Figure 1 | Cross-study direction of DEX-associated TM transcripts",
                 x=.08,y=.97,ha="left",color=NAVY,fontsize=15,fontweight="bold")
    fig.text(.5,.075,
      "Cells display the original study estimate and sign. The 7-day normalized CPM+1 ratio\n"
      "is NOT directly effect-size comparable with the two 48-hour differential-expression log2FC values.\n"
      "ACTA2 changes sign in the 7-day study; no randomized dose × duration contrast or protein assay.",
      ha="center",va="top",fontsize=9,color=GREY)
    save(fig,dest,"Fig1_three_study_core_transcript_signs")

def label_figure(cut,dual,dest):
    a={r["gene_symbol"]:r for r in cut}
    b={r["gene_symbol"]:r for r in dual}
    if len(a)!=16 or len(b)!=16 or set(a)!=set(b):raise ValueError("PANEL_16_MISMATCH")
    fig,axs=plt.subplots(1,2,figsize=(13,6.1),sharey=True)
    fig.subplots_adjust(left=.12,right=.98,top=.63,bottom=.27,wspace=.18)
    left=("score075_log2fc","score050_log2fc")
    right=("BeamA_55_raw_edgeR_log2FC","TM_Fibro_55_raw_edgeR_log2FC")
    for ax,keys,source,toplab,title in (
        (axs[0],left,a,["BeamA ≥0.75","BeamA ≥0.50"],
          "A. Classification confidence — same 57 donors"),
        (axs[1],right,b,["BeamA ≥0.75","scPred TM_Fibro"],
          "B. Alternative annotation — same 55 donors")):
        for j,key in enumerate(keys):
            for i,g in enumerate(CORE):
                v=num(source[g][key])
                ax.add_patch(Rectangle((j-.46,i-.43),.92,.86,
                                     facecolor=sign_color(v),alpha=.14,
                                     edgecolor="#d0d7dc"))
                ax.text(j,i-.02,("+" if v>=0 else "−"),fontsize=19,
                    fontweight="bold",ha="center",color=sign_color(v))
                ax.text(j,i+.28,f"{v:+.3f}",fontsize=10,ha="center",color=NAVY)
        ax.set_xlim(-.51,1.51);ax.set_ylim(3.52,-.54)
        ax.set_yticks(range(4),CORE,fontweight="bold")
        ax.set_xticks((0,1),toplab);ax.xaxis.tick_top()
        ax.tick_params(length=0,pad=8)
        ax.set_title(title,fontsize=12,color=NAVY,pad=56)
        for sp in ax.spines.values():sp.set_visible(False)
    fig.suptitle("Figure 2 | Effect directions are robust to tested cell-annotation choices",
       fontsize=15,fontweight="bold",x=.06,ha="left",color=NAVY)
    fig.text(.50,.07,
      "Both panels use donor-paired raw-UMI edgeR log2FC. They are NOT independent cohorts:\n"
      "the 57-donor cutoff comparison and 55-donor annotation comparison reuse original 10x libraries;\n"
      "the latter labels can include overlapping cells. All four genes have FDR < 0.05 in each displayed model.",
       fontsize=9,ha="center",color=GREY)
    save(fig,dest,"Fig2_BeamA_cutoff_and_matched55_annotation")
def emt_bar(rows,dest):
    if len(rows)!=194 or len({r["gene_id"] for r in rows})!=194:
        raise ValueError("EMT_194_ID_INTEGRITY")
    pairs=[]
    for r in rows:
        a=num(r["acute_500nM_48h_log2FC"])
        b=num(r["external_100nM_48h_log2FC"])
        c=num(r["third_500nM_7d_CPM_log2ratio"])
        pairs.append((a*b>0,a*c>0,b*c>0,(a*b>0 and b*c>0)))
    counts=[sum(x[k] for x in pairs) for k in range(4)]
    if counts != [141,156,145,124]:
        raise ValueError("EMT_FROZEN_EXPECTED_COUNTS_DRIFT "+str(counts))
    labels=["500 nM 48 h vs 100 nM 48 h",
            "500 nM 48 h vs 500 nM 7 d",
            "100 nM 48 h vs 500 nM 7 d",
            "All three experimental settings"]
    fig,ax=plt.subplots(figsize=(11.5,5.3))
    fig.subplots_adjust(left=.36,right=.96,top=.80,bottom=.28)
    y=list(range(4))
    colors=[BLUE,BLUE,BLUE,GREEN]
    ax.barh(y,counts,color=colors,height=.59)
    ax.set_yticks(y,labels)
    ax.invert_yaxis()
    ax.set_xlim(0,220)
    ax.set_xlabel("EMT genes with matching signed RNA effects (of 194 tested)")
    for i,v in enumerate(counts):ax.text(v+3,i,f"{v}/194   {v/194:.1%}",va="center",color=NAVY,fontweight="bold")
    ax.set_title("Figure 3 | Signed cross-experiment concordance of Hallmark EMT genes",loc="left",color=NAVY,pad=15)
    fig.text(.51,.055,
      "Descriptive fractions only: 194 post-hoc selected EMT genes are correlated, not independent Bernoulli trials.\n"
      "Shared doses, culture conditions and biological replication are NOT equivalent across studies.\n"
      "Sign concordance alone does not validate EMT-like tissue fibrosis or aqueous-outflow changes.",
      ha="center",color=GREY,fontsize=9)
    save(fig,dest,"Fig3_EMT_194_signed_cross_experiment_concordance")
    return counts
def supplements(cut,dual,core,dest):
    a={r["gene_symbol"]:r for r in cut};b={r["gene_symbol"]:r for r in dual}
    t={r["gene_symbol"]:r for r in core}
    rows=[]
    for gene in ("MYOC","TGFBI","MMP2","ACTA2","SDC1","EGR3","LPAR1","TFPI",
       "LUM","IGFBP3","VCAM1","NTM","FKBP5","CHAC1","ERFE","PARM1"):
        r={"gene_symbol":gene,"BeamA57_score075_edgeR_log2FC":a[gene]["score075_log2fc"],
            "BeamA57_score050_edgeR_log2FC":a[gene]["score050_log2fc"],
            "BeamA55_edgeR_log2FC":b[gene]["BeamA_55_raw_edgeR_log2FC"],
            "scPred55_edgeR_log2FC":b[gene]["TM_Fibro_55_raw_edgeR_log2FC"],
            "scPred55_fdr":b[gene]["TM_Fibro_55_FDR"],
            "scPred55_filter_status":b[gene]["filter_status"],
            "external_100nM_48h_log2FC":t[gene]["published_100nM_48h_log2FC"] if gene in t else "NOT_IN_PINNED_CORE_SOURCE",
            "other_lab_500nM_7d_CPM_plus1_log2ratio":t[gene]["chronic_500nM_7d_CPM_log2ratio"] if gene in t else "NOT_IN_PINNED_CORE_SOURCE",
            "same_cohort_models_not_independent":True,
            "cross_study_7d_ratio_not_same_estimand_as_edgeR_log2FC":True}
        rows.append(r)
    with (dest/"Supplementary_Table_1_evidence_and_filters.tsv").open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--third",required=True,type=Path)
    p.add_argument("--cutoff",required=True,type=Path)
    p.add_argument("--dual",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    paths={}
    for key,name in SOURCE_FILES.items():
        d=a.third if key.startswith("three_") or key=="third_gate" else a.cutoff if key.startswith("cutoff_") else a.dual
        path=d/name
        if not path.is_file():raise FileNotFoundError(str(path))
        sha=hashlib.sha256(path.read_bytes()).hexdigest()
        if key in PIN and sha != PIN[key]:raise ValueError("IMMutable SOURCE_HASH_DRIFT "+key)
        paths[key]={"sha256":sha,"basename":name}
    third=json.loads((a.third/SOURCE_FILES["third_gate"]).read_text())
    cut=json.loads((a.cutoff/SOURCE_FILES["cutoff_gate"]).read_text())
    dual=json.loads((a.dual/SOURCE_FILES["dual_gate"]).read_text())
    if third["status"]!="PASS_SOURCE_INTEGRITY_AND_DESCRIPTIVE_CROSS_STUDY_CPM" or third["emt_all_three_studies_same_signed_n"]!=124:
        raise ValueError("THIRD_COHORT_GATE_NOT_VERIFIED")
    if cut["status"]!="PASS_GATE_B_SCORE050_VS075_COUNT_LEVEL_EDGER_SENSITIVITY" or cut["core_three_same_sign_n"]!=3:
        raise ValueError("CUTOFF_GATE_NOT_VERIFIED")
    if dual["status"]!="PASS_MATCHED_55_DONOR_DUAL_ANNOTATION_EDGER_SENSITIVITY" or dual["core_3_FDR_both_lt005"]!=3 or dual["marker_16_filtered_from_one_annotation"]!=["ERFE"]:
        raise ValueError("MATCHED55_GATE_NOT_VERIFIED")
    core=table(a.third/SOURCE_FILES["three_core"])
    emt=table(a.third/SOURCE_FILES["three_emt"])
    cutoff=table(a.cutoff/SOURCE_FILES["cutoff_panel"])
    compare=table(a.dual/SOURCE_FILES["dual_panel"])
    if (len(core),len(emt),len(cutoff),len(compare))!=(4,194,16,16):
        raise ValueError("FIGURE_INPUT_ROW_COUNT_CHANGED")
    style();study_heatmap(core,a.out);label_figure(cutoff,compare,a.out)
    counts=emt_bar(emt,a.out)
    supplements(cutoff,compare,core,a.out)
    legend=(
      "# Thesis Figures — GSE313736 / Mehrotra 2026 / GSE301517\\n\\n"
      "**Figure 1.** Direction of the DEX-associated transcript effect in three different TM experimental settings. "
      "Data from 500nM 48h donor-paired BeamA edgeR; official published 100nM 48h bulk TM log2FC; "
      "and 500nM 7d CPM log2((DEX+1)/(Control+1)) averaged equally over TM62, TM70 and TM72 strains. "
      "Positive orange cells/negative blue cells show the *sign* only. Numeric values are original study estimates, "
      "not jointly estimated treatment effects. ACTA2 is positive in the 48h cohorts but negative in the 7d cohort. "
      "No independent donor-matched time×dose trial or protein endpoint.\\n\\n"
      "**Figure 2.** Within-source cell-annotation sensitivity of original GSE313736. Panel A: same 57 paired donors, "
      "BeamA confidence 0.75 versus 0.50, raw-count paired edgeR. Panel B: same **55** donor pairs under "
      "BeamA 0.75 versus scPred TM_Fibro, separately fitted raw-count edgeR. Samples, libraries and some annotated cells overlap; "
      "these are not independent cohort replications. All four displayed core/context genes have gene-level BH-FDR<0.05 in both "
      "models per panel. MYOC/TGFBI/MMP2 are three-cohort signed transcript candidates; ACTA2 is not.\\n\\n"
      "**Figure 3.** Descriptive, post-hoc Hallmark EMT signed concordance across 194 Ensembl-matched genes: "
      "141/194 acute studies, 156/194 original vs 7d, 145/194 published 100nM vs 7d, "
      "124/194 all three. These gene-sign counts do not constitute independent testing, pathway-replication p-values "
      "or functional proof of matrix deposition/stiffness/aqueous outflow.\\n\\n"
      "**Supplementary Table S1.** 16 prespecified candidate genes with 57-donor BeamA cutoff sensitivity, 55-donor "
      "BeamA/scPred same-donor effects and explicit ERFE filtered-out status. Cross-study values are present only "
      "for the original fixed four-core comparator; no numerical values are fabricated for the other 12 candidates.\\n\\n"
      "Sources: original public GSE313736 (57-donor BeamA and 55-donor annotation models), "
      "Mehrotra et al. IOVS 2026 100nM 48h official supplementary genes, "
      "GSE301517 500nM 7d published CPM matrix. "
      "Computational source archive: GitHub Actions run 41186736618/artifact 13037335888, "
      "run 41199532259/artifact 13041944372, "
      "run 41207583296/artifact 13045350639.\\n\\n"
      "**Scientific conclusion:** priorities for prospective mechanistic experiments, NOT new causal gene discoveries "
      "or individual steroid-induced ΔIOP predictors.\\n")
    (a.out/"FIGURE_LEGENDS_AND_INTERPRETATION.md").write_text(legend.replace("\\n","\n"),encoding="utf8")
    meta={"status":"PASS_FROZEN_ARTIFACT_FIGURES",
      "source_artifacts":{"three_studies":13037335888,"BeamA_cutoffs":13041944372,"dual_55":13045350639},
      "source_file_hashes":paths,
      "all_three_emt_signed":counts[3],"acute_emt_signed":counts[0],
      "core_three_3studies_signed":["MYOC","TGFBI","MMP2"],
      "ACTA2_7day_direction_discordant":True,
      "same_55_donor_sensitivity":True,
      "ERFE_scPred_not_tested":True,
      "feature_provenance":"one analysis per source; 7d CPM+1 ratio unequal estimand",
      "not_independent_RNA_protein_or_IOP_mechanism":True}
    (a.out/"FIGURE_SOURCE_MANIFEST.json").write_text(json.dumps(meta,indent=2)+"\n",encoding="utf8")
    print("THESIS_FIGURES",json.dumps({"status":meta["status"],
         "files":sorted(z.name for z in a.out.iterdir()),"emt":counts},ensure_ascii=False),flush=True)
if __name__=="__main__":main()
