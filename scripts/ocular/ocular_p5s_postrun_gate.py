#!/usr/bin/env python3
"""P5s strict aggregate-only verification of successful P5q and P5r artifacts.

Recompute every 12 x 2 x 9 leave-pool-out summary from independent public rows.
No raw counts, donor identifiers, or new biological replication.
"""
import argparse
import csv
import json
import math
from collections import Counter
from itertools import product
from pathlib import Path

LOW={"ERFE","CHAC1","PARM1"}
ADEQUATE={"ADAMTS15","EGR3","LPAR1","ABHD3","DDIT4L","TFPI","STK17B","RAB32","SDC1"}
TARGETS=LOW|ADEQUATE
MODELS={"paired_donor","paired_pool_library"}
POOLS=set("ABCDEFGHI")

def getrows(file):
    with open(file,encoding="utf8",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))
def val(v):
    try:
        x=float(v)
    except (ValueError,TypeError):
        raise ValueError("UNPARSEABLE_NUMBER "+str(v))
    if not math.isfinite(x):raise ValueError("NONFINITE_NUMBER "+str(v))
    return x
def truth(v):
    if v in ("TRUE","True","true","1"):return True
    if v in ("FALSE","False","false","0"):return False
    raise ValueError("UNPARSEABLE_BOOL "+str(v))
def close(a,b):
    return math.isclose(a,b,rel_tol=1e-6,abs_tol=1e-9)
def save_tsv(path,rows):
    if not rows:raise ValueError("EMPTY_OUTPUT")
    with open(path,"w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        w.writeheader();w.writerows(rows)
def main():
    ap=argparse.ArgumentParser()
    for arg in ("p5q-detail","p5q-summary","p5q-json","p5r-consensus","p5r-timecourse","p5r-threshold-grid","output"):
        ap.add_argument("--"+arg,required=True,type=Path)
    a=ap.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    d=getrows(a.p5q_detail)
    s=getrows(a.p5q_summary)
    r=getrows(a.p5r_consensus)
    t=getrows(a.p5r_timecourse)
    grid=getrows(a.p5r_threshold_grid)
    j=json.loads(a.p5q_json.read_text(encoding="utf8"))
    if len(d)!=216 or len(s)!=24 or len(r)!=12 or len(t)!=12 or len(grid)!=24:
        raise ValueError("CONTRACT_DIMENSIONS_CHANGED")
    if j.get("status")!="PASS_18_FORMAL_EDGE_R_QL_REFITS" or j.get("refits")!=18:
        raise ValueError("P5Q_NOT_ACTUALLY_COMPLETE")
    keys={(x["gene_symbol"],x["model"],x["omitted_pool"]) for x in d}
    expected=set(product(TARGETS,MODELS,POOLS))
    if keys!=expected or len(keys)!=216:raise ValueError("MISSING_OR_DUPLICATED_TARGET_MODEL_POOL")
    keys2={(x["gene_symbol"],x["model"]) for x in s}
    if keys2!=set(product(TARGETS,MODELS)):raise ValueError("MISSING_MODEL_SUMMARIES")
    if {x["gene_symbol"] for x in r}!=TARGETS or {x["gene_symbol"] for x in t}!=TARGETS:
        raise ValueError("P5R_SHORTLIST_DRIFT")
    ref={x["gene_symbol"]:x for x in r}
    ref_time={x["gene_symbol"]:x for x in t}
    if len({x["gene_id"] for x in r})!=12 or len({x["gene_id"] for x in t})!=12:
        raise ValueError("ENSMBL_TARGET_IDS_NOT_UNIQUE")
    gridbase=[x for x in grid if close(val(x["logCPM_cutoff"]),2) and
                                   close(val(x["min_abs_log2fc_all_three"]),.5)]
    if len(gridbase)!=1 or int(gridbase[0]["n_genes"])!=9 or
       set(gridbase[0]["genes"].split(","))!=ADEQUATE:
        raise ValueError("ORIGINAL_NINE_SCREEN_NOT_REPRODUCED")
    bysummary={(x["gene_symbol"],x["model"]):x for x in s}
    output=[]
    for symbol in sorted(TARGETS):
        xref=ref[symbol]
        if xref["gene_id"]!=ref_time[symbol]["gene_id"]:
            raise ValueError("P5R_ENSMBL_JOIN_CONFLICT "+symbol)
        if not close(val(xref["primary_log2fc"]),val(ref_time[symbol]["primary_500nM_48h_log2fc"])):
            raise ValueError("CROSSWALK_PRIMARY_LFC_DRIFT "+symbol)
        for model in sorted(MODELS):
            items=[x for x in d if x["gene_symbol"]==symbol and x["model"]==model]
            if len(items)!=9 or {x["omitted_pool"] for x in items}!=POOLS:
                raise ValueError("INVALID_NINE_OMISSIONS "+symbol+" "+model)
            for y in items:
                if y["gene_id"]!=xref["gene_id"]:
                    raise ValueError("P5Q_P5R_ENSMBL_MISMATCH "+symbol)
                if not close(val(y["primary_log2fc"]),val(xref["primary_log2fc"])):
                    raise ValueError("PRIMARY_LFC_CROSS_SOURCE_MISMATCH "+symbol)
                if not 0<=val(y["FDR"])<=1 or not 0<=val(y["PValue"])<=1:
                    raise ValueError("INVALID_QL_FDR_PVALUE "+symbol)
            valid=[x for x in items if truth(x["tested"])]
            agreement=sum(truth(x["primary_sign_agree"]) for x in valid)
            sig=sum(val(x["FDR"])<.05 for x in valid)
            together=sum(truth(x["primary_sign_agree"]) and val(x["FDR"])<.05 for x in valid)
            worst=max(val(x["FDR"]) for x in valid) if valid else None
            minabs=min(abs(val(x["log2fc"])) for x in valid) if valid else None
            summary=bysummary[(symbol,model)]
            checks=[
                int(summary["n_refits"])==9,
                int(summary["n_tested"])==len(valid),
                int(summary["same_sign"])==agreement,
                int(summary["fdr_below_0_05"])==sig,
                int(summary["significant_and_same_sign"])==together,
                truth(summary["all_nine_sign_stable"])==(len(valid)==9 and agreement==9),
                truth(summary["all_nine_fdr_and_sign"])==(len(valid)==9 and together==9),
                (worst is None or close(val(summary["max_fdr"]),worst)),
                (minabs is None or close(val(summary["min_abs_log2fc"]),minabs))]
            if not all(checks):raise ValueError("P5Q_QL_SUMMARY_RECOMPUTE_FAILED "+symbol+" "+model)
            output.append({"gene_id":xref["gene_id"],"gene_symbol":symbol,
                "expression_tier":"LOW_EXPRESSION" if symbol in LOW else "ADEQUATE_EXPRESSION",
                "primary_mean_logCPM":val(xref["logCPM"]),
                "model":model,"tested_n":len(valid),"sign_stable_n":agreement,
                "BH_FDR_sig_n":sig,"BH_FDR_sig_and_signed_n":together,
                "max_refit_FDR":worst,"min_refit_abs_log2fc":minabs,
                "all_nine_edgeR_QL_signed_and_sig":len(valid)==9 and together==9,
                "descriptive_pool_sign_count":int(xref["pool_signed_n"]),
                "bootstrap_is_descriptive":True,
                "independent_replication":False,
                "dose_specific_mechanism_validated":False,
                "clinical_delta_IOP_validated":False})
    nall=sum(x["all_nine_edgeR_QL_signed_and_sig"] for x in output)
    if nall!=24:raise ValueError("P5Q_ANY_NINEFOLD_MODEL_FAILED "+str(nall))
    if j.get("donor_model_full_9_sig_sign")!=12 or j.get("pooled_model_full_9_sig_sign")!=12:
        raise ValueError("P5Q_JSON_VS_ROWS_MISMATCH")
    save_tsv(a.output/"P5s_12_gene_two_model_external_reconciliation.tsv",output)
    verdict={
        "status":"PASS_P5Q_P5R_AGGREGATE_CROSS_ARTIFACT_RECOMPUTATION",
        "p5q_run":41170641379,"p5r_run":41172214095,
        "p5q_detail_rows_verified":216,"p5q_gene_model_rows_verified":24,
        "formal_ql_refit_models":18,"formal_model_gene_combinations":216,
        "all_12_have_9_of_9_BH_FDR_and_primary_sign_in_both_models":True,
        "adequate_expression_nine":sorted(ADEQUATE),"low_expression_three":sorted(LOW),
        "prior_96_gene_pool_descriptive_recon_still_not_independent":True,
        "drug_vs_library_batch_causality":"HOLD",
        "matched_100_500nM_dose_time_interaction":"HOLD",
        "protein_orthogonal_mechanistic_novelty":"HOLD",
        "clinical_steroid_delta_IOP_prediction":"STOP",
        "note":"Post-hoc selected genes and 18 leave-out refits are correlated re-analyses of one study. 12/12 significant in both models is internal numerical robustness, not 12 externally replicated causal genes."
    }
    (a.output/"P5s_gate_verdict.json").write_text(
        json.dumps(verdict,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("P5S_PASS 216/216 rows, 24/24 summaries and 12/12 strict dual-model nine-LOO signed FDR criteria",flush=True)
if __name__=="__main__":
    main()
