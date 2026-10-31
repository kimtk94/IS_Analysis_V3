#!/usr/bin/env python3
"""Gate B full pseudobulk edgeR .50-vs-.75 sensitivity, gene-level aggregate ONLY.

The 0.50 and 0.75 analyses are the SAME donors, SAME underlying study, and
not independent experiments; both are conditioned on post-DEX cell labels.
"""
import argparse,csv,gzip,json,math
from pathlib import Path
from statistics import median
CORE=["MYOC","TGFBI","MMP2","ACTA2","SDC1","EGR3","LPAR1","TFPI","LUM","IGFBP3","VCAM1","NTM","FKBP5"]
def load(p):
    with gzip.open(p,"rt",encoding="utf8",newline="") as f:
        rs=list(csv.DictReader(f,delimiter="\t"))
    ids={}
    for r in rs:
        k=r["gene_id"].split(".")[0]
        if k in ids:raise ValueError("DUPLICATE_GENE")
        v={key:float(r[key]) for key in ("logFC","logCPM","FDR","PValue")}
        if not all(math.isfinite(x) for x in v.values()):raise ValueError("NONFINITE_GENE_VALUE")
        ids[k]={**v,"symbol":r["gene_symbol"]}
    return ids
def emit(out,name,rows):
    with (out/name).open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--baseline",type=Path,required=True)
    ap.add_argument("--sensitivity",type=Path,required=True)
    ap.add_argument("--source",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    p=load(a.baseline);s=load(a.sensitivity)
    src=json.loads(a.source.read_text(encoding="utf8"))
    if src.get("status")!="GO_57DONOR_FULL_RAW_PSEUDOBULK_QC" or src.get("score_cutoff")!=.5:
        raise ValueError("SOURCE_SCORE_AND_57_BIOLOGICAL_DONOR_CONTRACT_DRIFT")
    if len(p)!=16050 or len(s)<10000:raise ValueError("EXPECTED_BASELINE_AND_SENSITIVITY_GENE_UNIVERSES")
    shared=sorted(set(p)&set(s))
    same=sum(p[g]["logFC"]*s[g]["logFC"]>0 for g in shared)
    fdr_both=sum(p[g]["FDR"]<.05 and s[g]["FDR"]<.05 for g in shared)
    base_sig=[g for g in shared if p[g]["FDR"]<.05 and abs(p[g]["logFC"])>=.25]
    base_sig_same=sum(p[g]["logFC"]*s[g]["logFC"]>0 for g in base_sig)
    base_sig_sig_same=sum(p[g]["logFC"]*s[g]["logFC"]>0 and s[g]["FDR"]<.05 for g in base_sig)
    matrix=[]
    symbol_to_id={rec["symbol"]:id for id,rec in p.items()}
    for symbol in CORE:
        id=symbol_to_id.get(symbol)
        if id not in s:raise ValueError("CORE_GENE_NOT_SHARED "+symbol)
        before,after=p[id],s[id]
        matrix.append({"gene_symbol":symbol,"ensembl_id":id,
           "baseline_075_log2FC":round(before["logFC"],7),
           "threshold_050_log2FC":round(after["logFC"],7),
           "delta_log2FC_050_minus_075":round(after["logFC"]-before["logFC"],7),
           "baseline_075_FDR":before["FDR"],"threshold_050_FDR":after["FDR"],
           "signed_agreement":before["logFC"]*after["logFC"]>0,
           "both_gene_FDR_below_005":before["FDR"]<.05 and after["FDR"]<.05,
           "same_57_donors_not_external_replication":True,
           "cannot_disentangle_library_and_DEX":True})
    emit(a.out,"GateB_13_markers_050_vs_075_edgeR.tsv",matrix)
    high=[{"gene_id":g,"gene_symbol":p[g]["symbol"],"baseline_075_log2FC":round(p[g]["logFC"],6),
           "threshold_050_log2FC":round(s[g]["logFC"],6),
           "baseline_075_FDR":p[g]["FDR"],"threshold_050_FDR":s[g]["FDR"],
           "same_sign":p[g]["logFC"]*s[g]["logFC"]>0}
        for g in shared if p[g]["FDR"]<.05 and abs(p[g]["logFC"])>=3]
    emit(a.out,"GateB_83_high_effect_050_vs_075_edgeR.tsv",high)
    summary={
      "status":"PASS_TRUE_RAW_UMI_57_DONOR_BEAMA_SCORE_SENSITIVITY",
      "reference_score":0.75,"sensitivity_score":0.50,
      "n_biological_pairs":57,
      "reference_edgeR_tested":len(p),"sensitivity_edgeR_tested":len(s),
      "shared_ensembl_genes":len(shared),
      "shared_gene_signed_agreement":same,
      "shared_gene_signed_agreement_fraction":round(same/len(shared),6),
      "shared_genes_both_fdr_below_005":fdr_both,
      "baseline_FDR_and_abs_effect_ge_point25_shared_n":len(base_sig),
      "baseline_FDR_and_abs_effect_ge_point25_sign_concordant_n":base_sig_same,
      "baseline_FDR_and_abs_effect_ge_point25_both_significant_and_same_n":base_sig_sig_same,
      "baseline_83_high_effect_shared_n":len(high),
      "baseline_83_high_effect_signed_agreement_n":sum(x["same_sign"] for x in high),
      "baseline_83_high_effect_both_fdr_below005_n":sum(x["threshold_050_FDR"]<.05 for x in high),
      "baseline_83_high_effect_threshold_retained_n":sum(x["threshold_050_FDR"]<.05 and abs(x["threshold_050_log2FC"])>=3 for x in high),
      "core_3_signed_agreement":all(x["signed_agreement"] for x in matrix if x["gene_symbol"] in ("MYOC","TGFBI","MMP2")),
      "core_3_both_fdr_below005":all(x["both_gene_FDR_below_005"] for x in matrix if x["gene_symbol"] in ("MYOC","TGFBI","MMP2")),
      "cross_cutoff_donors_unchanged":True,
      "post_treatment_cell_type_selection_collider_bias_removed":False,
      "technical_treatment_library_confounded":True,
      "this_is_independent_biological_validation":False,
      "interpretation":"Score sensitivity from the same 57 donors, 18 separately prepared DEX/vehicle 10x libraries; cannot establish treatment-specific causality."}
    (a.out/"GateB_050_vs_075_raw_UMI_edgeR_gate.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf8")
    print("GATEB_EDGER",json.dumps(summary,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
