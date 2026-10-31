#!/usr/bin/env python3
"""Fixed common 55 donor BeamA vs TM_Fibro gene effect contrast, aggregate-only."""
import argparse,csv,gzip,json,math
from pathlib import Path
import numpy as np
MARKERS=("MYOC","TGFBI","MMP2","ACTA2","SDC1","EGR3","LPAR1","TFPI","LUM","IGFBP3","VCAM1","NTM","FKBP5","CHAC1","ERFE","PARM1")
CORE=("MYOC","TGFBI","MMP2")
def load(path):
    with gzip.open(path,"rt",encoding="utf8",newline="") as f:
        rd=csv.DictReader(f,delimiter="\t")
        if not {"gene_id","gene_symbol","logFC","FDR","logCPM"}.issubset(rd.fieldnames or []):
            raise ValueError("MISSING_GENE_STAT_FIELDS")
        out={}
        for r in rd:
            id=r["gene_id"].split(".")[0]
            if id in out:raise ValueError("DUP_GENE_ID")
            vals={n:float(r[n]) for n in ("logFC","FDR","logCPM")}
            if not all(math.isfinite(x) for x in vals.values()):raise ValueError("NONFINITE_GENE_VALUE")
            out[id]={"gene_symbol":r["gene_symbol"],**vals}
    return out
def table(path,rows):
    if not rows:raise ValueError("EMPTY_TABLE")
    with path.open("w",newline="",encoding="utf8") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        writer.writeheader();writer.writerows(rows)
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--beam",type=Path,required=True)
    p.add_argument("--fibro",type=Path,required=True)
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    q=json.loads(args.source.read_text(encoding="utf8"))
    if q.get("status")!="PASS_55_DUAL_ANNOTATION_RAW_PSEUDOBULK" or q.get("matrix_joined_common_donors_ge30_every_arm_and_annotation")!=55:
        raise ValueError("NOT_MATCHED_55")
    a=load(args.beam);b=load(args.fibro)
    common=sorted(set(a)&set(b))
    if len(common)<12000:raise ValueError("LOW_COMMON_GENE_COUNT")
    effects=np.array([a[g]["logFC"] for g in common]),np.array([b[g]["logFC"] for g in common])
    x,y=effects
    corr=float(np.corrcoef(x,y)[0,1])
    if not math.isfinite(corr):raise ValueError("CORRELATION_NONFINITE")
    rows=[]
    for id in common:
        z=a[id];t=b[id]
        if z["gene_symbol"]!=t["gene_symbol"]:raise ValueError("GENE_ID_NAME_DRIFT")
        if z["gene_symbol"] not in MARKERS:continue
        rows.append({"gene_id":id,"gene_symbol":z["gene_symbol"],
          "BeamA_55_raw_edgeR_log2FC":round(z["logFC"],7),
          "TM_Fibro_55_raw_edgeR_log2FC":round(t["logFC"],7),
          "BeamA_55_FDR":z["FDR"],"TM_Fibro_55_FDR":t["FDR"],
          "sign_agrees":z["logFC"]*t["logFC"]>0,
          "FDR_both_lt_005":z["FDR"]<.05 and t["FDR"]<.05,
          "same_55_donors_verified":True,
          "cells_can_overlap_across_annotations":True,
          "independent_study_replication":False,
          "batch_effect_independent_treatment_identified":False})
    # An annotation with fewer selected cells may filter genuinely low-CPM genes.
    # This is NOT a failure of the fitted models. Report all preselected targets
    # and mark untestable genes without inventing an FDR or forcing sign agreement.
    target_rows={r["gene_symbol"]:r for r in rows}
    raw_a={r["gene_symbol"]:r for r in a.values()}
    raw_b={r["gene_symbol"]:r for r in b.values()}
    if not all(g in target_rows for g in CORE):
        raise ValueError("THREE_CORE_UNTESTABLE_IN_ONE_ANNOTATION")
    complete=[]
    for gene in MARKERS:
        r=target_rows.get(g)
        if r is not None:
            complete.append({**r,"tested_in_both_annotations":True,
              "filter_status":"TESTED_IN_BOTH"})
        else:
            complete.append({"gene_id":"","gene_symbol":gene,
              "BeamA_55_raw_edgeR_log2FC":round(raw_a[gene]["logFC"],7) if gene in raw_a else "",
              "TM_Fibro_55_raw_edgeR_log2FC":round(raw_b[gene]["logFC"],7) if gene in raw_b else "",
              "BeamA_55_FDR":raw_a[gene]["FDR"] if gene in raw_a else "",
              "TM_Fibro_55_FDR":raw_b[gene]["FDR"] if gene in raw_b else "",
              "sign_agrees":"","FDR_both_lt_005":"",
              "same_55_donors_verified":True,
              "cells_can_overlap_across_annotations":True,
              "independent_study_replication":False,
              "batch_effect_independent_treatment_identified":False,
              "tested_in_both_annotations":False,
              "filter_status":"FILTERED_FROM_TM_FIBRO" if gene not in raw_b else "FILTERED_FROM_BEAMA"})
    if len(complete)!=16 or len({r["gene_symbol"] for r in complete})!=16:
        raise ValueError("PANEL_16_UNIQUE_EXPECTED")
    table(args.out/"GATE_B_55_matched_BeamA_vs_scPred_16_genes.tsv",complete)
    note=target_rows
    summary={"status":"PASS_MATCHED_55_DONOR_DUAL_ANNOTATION_EDGER_SENSITIVITY",
     "biological_pairs_same_donors":55,"source_libraries":18,
     "beam_tested_genes":len(a),"fibro_tested_genes":len(b),
     "common_Ensembl_ids":len(common),
     "shared_gene_effect_pearson":round(corr,8),
     "shared_genes_same_logFC_direction":int(np.count_nonzero(x*y>0)),
     "shared_genes_both_FDR_lt005":sum(a[g]["FDR"]<.05 and b[g]["FDR"]<.05 for g in common),
     "beam_FDR_lt005":sum(z["FDR"]<.05 for z in a.values()),
     "fibro_FDR_lt005":sum(z["FDR"]<.05 for z in b.values()),
     "marker_16_testable_both_n":len(rows),
     "marker_16_filtered_from_one_annotation":[z["gene_symbol"] for z in complete if not z["tested_in_both_annotations"]],
     "marker_16_same_sign":sum(z["sign_agrees"] for z in rows),
     "marker_16_FDR_both_lt005":sum(z["FDR_both_lt_005"] for z in rows),
     "core_3_same_sign":sum(note[s]["sign_agrees"] for s in CORE),
     "core_3_FDR_both_lt005":sum(note[s]["FDR_both_lt_005"] for s in CORE),
     "core_3":list(CORE),
     "not_independent_replication":True,
     "not_randomized_steroid_dose_or_clinical_IOP":True,
     "treatment_library_prep_aliased":True,
     "cell_label_post_exposure":True,
     "annotation_groups_may_overlap_cells":True,
     "RNA_assay_only_no_protein_function":True}
    (args.out/"GATE_B_55_matched_annotation_concordance.json").write_text(
        json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("MATCHED_55_GATE",json.dumps(summary,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
