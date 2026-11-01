#!/usr/bin/env python3
"""Gate B actual 0.50 vs frozen 0.75 edgeR QL refit comparisons.

Genuine full raw-integer pseudobulk count refits are required as input, not
sample-level CPM or metadata-only fractions. Outputs gene-level statistics
and aggregate QC, never individual donor vectors, counts or names.
"""
import argparse,csv,gzip,json,math
from pathlib import Path
import numpy as np

PANEL=("MYOC","TGFBI","MMP2","ACTA2","SDC1","EGR3","LPAR1","TFPI",
       "LUM","IGFBP3","VCAM1","NTM","FKBP5","CHAC1","ERFE","PARM1")
CORE=("MYOC","TGFBI","MMP2")
def read(path):
    with gzip.open(path,"rt",encoding="utf8",newline="") as f:
        rd=csv.DictReader(f,delimiter="\t")
        if not {"gene_id","gene_symbol","logFC","FDR","logCPM"}.issubset(rd.fieldnames or []):
            raise ValueError("MISSING_EDGER_FIELDS")
        rows=list(rd)
    data={}
    for r in rows:
        gene=r["gene_id"].split(".")[0]
        if gene in data:raise ValueError("DUPLICATE_ENSEMBL_GENE")
        z={"gene_id":gene,"gene_symbol":r["gene_symbol"]}
        for k in ("logFC","FDR","logCPM"):
            z[k]=float(r[k])
            if not math.isfinite(z[k]):raise ValueError("NONFINITE_"+k)
        if not 0<=z["FDR"]<=1:raise ValueError("INVALID_FDR")
        data[gene]=z
    return data
def csvwrite(path,rows):
    if not rows:raise ValueError("EMPTY_OUTPUT")
    with path.open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        w.writeheader();w.writerows(rows)
def highset(d,cutoff=3):
    return {gene:r for gene,r in d.items() if r["FDR"]<.05 and abs(r["logFC"])>=cutoff}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--baseline-075",required=True,type=Path)
    p.add_argument("--new-050",required=True,type=Path)
    p.add_argument("--source-050",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    base=read(args.baseline_075);new=read(args.new_050)
    qc=json.loads(args.source_050.read_text(encoding="utf8"))
    if len(base)!=16050 or qc.get("score_cutoff")!=0.5:
        raise ValueError("SOURCE_SCORE_OR_BASELINE_DRIFT")
    if qc.get("metadata_paired_donors_after_actual_barcode_join")!=57:
        raise ValueError("EXPECTED_ACTUAL_57_DONORS_FOR_DIRECT_CUTOFF_COMPARISON")
    if qc.get("status")!="GO_57DONOR_FULL_RAW_PSEUDOBULK_QC":
        raise ValueError("SCORE050_FULL_RAW_COUNTS_NOT_VERIFIED")
    common=sorted(set(base)&set(new))
    if len(common)<14000:raise ValueError("INSUFFICIENT_GENE_UNIVERSE_OVERLAP")
    diffs=[]
    for gid in common:
        a,b=base[gid],new[gid]
        if a["gene_symbol"]!=b["gene_symbol"]:
            raise ValueError("GENE_SYMBOL_CHANGED_"+gid)
        diffs.append({"gene_id":gid,"gene_symbol":a["gene_symbol"],
          "logFC_score075":a["logFC"],"logFC_score050":b["logFC"],
          "FDR_score075":a["FDR"],"FDR_score050":b["FDR"],
          "mean_logCPM_score075":a["logCPM"],"mean_logCPM_score050":b["logCPM"],
          "same_sign":a["logFC"]*b["logFC"]>0,
          "both_FDR_lt_0_05":a["FDR"]<.05 and b["FDR"]<.05,
          "abs_logFC_delta":abs(a["logFC"]-b["logFC"])})
    x=np.asarray([base[k]["logFC"] for k in common],dtype=float)
    y=np.asarray([new[k]["logFC"] for k in common],dtype=float)
    corr=float(np.corrcoef(x,y)[0,1])
    if not math.isfinite(corr):raise ValueError("CORRELATION_NONFINITE")
    bysymbol={x["gene_symbol"]:x for x in diffs}
    if not set(PANEL).issubset(bysymbol):
        raise ValueError("SCORE_CUTOFF_DROPPED_P9_PANEL_GENE")
    markers=[{"gene_symbol":g,"gene_id":bysymbol[g]["gene_id"],
      "score075_log2fc":round(bysymbol[g]["logFC_score075"],6),
      "score050_log2fc":round(bysymbol[g]["logFC_score050"],6),
      "score075_BH_FDR":bysymbol[g]["FDR_score075"],
      "score050_BH_FDR":bysymbol[g]["FDR_score050"],
      "same_effect_direction":bysymbol[g]["same_sign"],
      "significant_both_thresholds":bysymbol[g]["both_FDR_lt_0_05"],
      "absolute_difference_log2fc":round(bysymbol[g]["abs_logFC_delta"],6),
      "is_three_cohort_stable_prior_candidate":g in CORE,
      "cell_state_selection_confounding_eliminated":False}
      for g in PANEL]
    csvwrite(args.out/"GATE_B_16_candidate_score050_vs075_edgeR.tsv",markers)
    hs75=highset(base);hs50=highset(new)
    hrows=[]
    for gene in sorted(set(hs75)|set(hs50)):
        a=base.get(gene);b=new.get(gene)
        hrows.append({"gene_id":gene,"gene_symbol":(a or b)["gene_symbol"],
          "high_effect_FDR_075":gene in hs75,
          "high_effect_FDR_050":gene in hs50,
          "logFC_075":a["logFC"] if a else "",
          "logFC_050":b["logFC"] if b else "",
          "FDR_075":a["FDR"] if a else "",
          "FDR_050":b["FDR"] if b else "",
          "not_original_paper_80_gene_crosswalk":True})
    csvwrite(args.out/"GATE_B_internal_high_effect_075_050_union.tsv",hrows)
    cutoff_strata=[]
    for t in (0.25,0.5,1,2,3):
        bset={g for g,r in base.items() if r["FDR"]<.05 and abs(r["logFC"])>=t}
        nset={g for g,r in new.items() if r["FDR"]<.05 and abs(r["logFC"])>=t}
        cutoff_strata.append({"abs_effect_cutoff":t,
         "score075_FDR_lt_005_and_effect":len(bset),
         "score050_FDR_lt_005_and_effect":len(nset),
         "shared_loci":len(bset&nset),
         "score075_only":len(bset-nset),
         "score050_only":len(nset-bset)})
    csvwrite(args.out/"GATE_B_FDR_effect_cutoff_sensitivity.tsv",cutoff_strata)
    summary={
      "status":"PASS_GATE_B_SCORE050_VS075_COUNT_LEVEL_EDGER_SENSITIVITY",
      "analysis":"GSE313736 high confidence BeamA paired raw UMI edgeR TMM robust QL",
      "confounding_remains":"Treatment-specific 10x libraries and post-treatment classification",
      "baseline_score":0.75,"reanalysis_score":0.5,
      "baseline_tested_genes":len(base),"reanalysis_tested_genes":len(new),
      "shared_genes":len(common),
      "baseline_bh_fdr_lt_005":sum(r["FDR"]<.05 for r in base.values()),
      "reanalysis_bh_fdr_lt_005":sum(r["FDR"]<.05 for r in new.values()),
      "pearson_gene_level_log2FC_shared":round(corr,8),
      "median_abs_effect_difference":round(float(np.median(np.abs(x-y))),7),
      "shared_genes_matching_direction_n":sum(r["same_sign"] for r in diffs),
      "shared_FDR_lt_005_both_n":sum(r["both_FDR_lt_0_05"] for r in diffs),
      "core_three_same_sign_n":sum(bysymbol[g]["same_sign"] for g in CORE),
      "core_three_bh_fdr_lt_005_both_n":sum(bysymbol[g]["both_FDR_lt_0_05"] for g in CORE),
      "P9_16_panel_same_sign_n":sum(bysymbol[g]["same_sign"] for g in PANEL),
      "baseline_high_effect_abs_logFC_ge_3":len(hs75),
      "reanalysis_high_effect_abs_logFC_ge_3":len(hs50),
      "internal_high_effect_shared_loci":len(set(hs75)&set(hs50)),
      "high_effect_strata":cutoff_strata,
      "source_score050_joined_paired_donors":57,
      "source_score050_total_matrix_matched_cells":sum(v["matched_cells"] for v in qc["per_arm"].values()),
      "all_18_matrix_crc_qc":"PASS_IN_P4",
      "original_author_reported_high_effect_80_exact_genes":"HOLD_NO_AUTHORIZED_SOURCE_GENE_LIST",
      "not_independent_biological_replication":True,
      "not_clinical_steroid_delta_iop":True
    }
    if summary["baseline_bh_fdr_lt_005"]!=13717:
        raise ValueError("REFERENCE_13717_GENES_DRIFT")
    (args.out/"GATE_B_score050_075_actual_refit_concordance.json").write_text(
        json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("GATE_B_COMPARISON",json.dumps({k:summary[k] for k in (
        "status","reanalysis_tested_genes","pearson_gene_level_log2FC_shared",
        "core_three_same_sign_n","core_three_bh_fdr_lt_005_both_n",
        "reanalysis_high_effect_abs_logFC_ge_3")},ensure_ascii=False),flush=True)
if __name__=="__main__":main()
