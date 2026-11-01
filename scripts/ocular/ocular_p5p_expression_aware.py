#!/usr/bin/env python3
"""P5p expression-aware triage for externally opposite-sign TM genes.

Uses only verified public aggregate edgeR gene-level summary and P5m timecourse
manifest; no donor-level data; exploratory gene priority NOT causal discovery.
"""
import argparse,csv,gzip,json,math
from collections import Counter
from pathlib import Path

def read_gz(p):
    with gzip.open(p,"rt",encoding="utf8",newline="") as f:
        x=csv.DictReader(f,delimiter="\t")
        return list(x)

def read_tsv(p):
    with open(p,encoding="utf8",newline="") as f:
        x=csv.DictReader(f,delimiter="\t")
        return list(x)

def num(x):
    try:
        y=float(x)
        if math.isfinite(y):return y
    except (TypeError,ValueError):pass
    return None

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--primary",required=True)
    a.add_argument("--p5m",required=True)
    a.add_argument("--output-dir",required=True)
    args=a.parse_args()
    primary=read_gz(args.primary);time=read_tsv(args.p5m)
    if len(primary)!=16050 or len(time)!=96:raise ValueError("EXPECTED_SOURCE_GENE_UNIVERSE")
    byid={r["gene_id"]:r for r in primary}
    if len(byid)!=len(primary):raise ValueError("PRIMARY_GENE_ID_DUPLICATED")
    if len({r["gene_id"] for r in time})!=96:raise ValueError("CONFLICT_GENE_ID_DUPLICATED")
    enriched=[]
    for x in time:
        ref=byid.get(x["gene_id"])
        if ref is None:raise ValueError("P5M_GENE_NOT_FOUND_IN_PRIMARY")
        p=num(x["primary_BeamA_500nM_48h_log2fc"])
        r=num(ref["logFC"])
        if p is None or r is None or abs(p-r)>1e-4:
            raise ValueError("PRIMARY_EFFECT_DIVERGED")
        logcpm=num(ref["logCPM"])
        external=num(x["external_TM_100nM_48h_log2fc"])
        pool=num(x["conservative_poolunit_log2fc"])
        pool_fdr=num(x["conservative_poolunit_FDR"])
        if None in (logcpm,external,pool,pool_fdr):
            raise ValueError("P5M_MISSING_NUMERIC_VALUE")
        if not p*external<0:raise ValueError("NOT_EXTERNAL_OPPOSITE_EFFECT")
        common=(pool_fdr<.05 and pool*p>0 and
                abs(p)>=.5 and abs(external)>=.5 and abs(pool)>=.5)
        adequately_expressed=common and logcpm>=2
        strong_low=common and abs(p)>=1 and logcpm<2
        group=("ADEQUATE_EXPRESSION_CONFIRMED_DISCORDANCE" if adequately_expressed
               else "LARGE_EFFECT_LOW_EXPRESSION" if strong_low
               else "OTHER_DESCRIPTIVE_CONFLICT")
        xout=dict(x)
        xout["primary_mean_logCPM"]=round(logcpm,6)
        xout["abundance_criterion"]="PASS_GE2" if logcpm>=2 else "BELOW_LOGCPM2"
        xout["expression_aware_priority"]=group
        xout["clinical_marker_validated"]="NO"
        xout["mechanistic_novelty_confirmed"]="NO"
        enriched.append(xout)
    enriched.sort(key=lambda r:(
        {"ADEQUATE_EXPRESSION_CONFIRMED_DISCORDANCE":0,
         "LARGE_EFFECT_LOW_EXPRESSION":1,
         "OTHER_DESCRIPTIVE_CONFLICT":2}[r["expression_aware_priority"]],
        -abs(float(r["primary_BeamA_500nM_48h_log2fc"])),
        r["gene_id"]))
    counts=Counter(r["expression_aware_priority"] for r in enriched)
    target=[r for r in enriched if r["expression_aware_priority"]=="ADEQUATE_EXPRESSION_CONFIRMED_DISCORDANCE"]
    strong=[r for r in enriched if r["expression_aware_priority"]=="LARGE_EFFECT_LOW_EXPRESSION"]
    expected={"ADAMTS15","EGR3","LPAR1","ABHD3","DDIT4L","TFPI","STK17B","RAB32","SDC1"}
    if {r["gene_symbol"] for r in target}!=expected:
        raise ValueError("EXPRESSION_ADEQUATE_CANDIDATE_SET_DRIFT")
    if {r["gene_symbol"] for r in strong}!={"ERFE","PARM1","CHAC1"}:
        raise ValueError("STRONG_LOW_ABUNDANCE_SET_DRIFT")
    result={
      "status":"PASS_EXPRESSION_AWARE_PUBLIC_GENE_TRIAGE",
      "evidence":"GSE313736 P5b primary edgeR aggregate plus P5m official 100nM time course and pool-unit estimates",
      "n_original_discordant":96,
      "n_adequate_expression_nine":len(target),
      "adequately_expressed_candidates":[r["gene_symbol"] for r in target],
      "n_high_effect_low_expression":len(strong),
      "high_effect_low_abundance":[
        {"gene":r["gene_symbol"],"primary_logCPM":r["primary_mean_logCPM"],
         "primary_log2fc":num(r["primary_BeamA_500nM_48h_log2fc"])}
        for r in strong],
      "screen_criteria":{
          "primary_mean_logCPM_ge":2,
          "primary_abs_log2fc_ge":.5,
          "external_abs_log2fc_ge":.5,
          "conservative_poolunit_abs_log2fc_ge":.5,
          "conservative_poolunit_BH_FDR_lt":.05,
          "signed_stability_primary_vs_poolunit":True},
      "priority_counts":dict(counts),
      "caveats":[
          "Post hoc hypothesis filter; no prospective registration",
          "Mean logCPM is a broad abundance metric, not protein expression or cell-identity validation",
          "High-expression opposing effects are still potentially confounded with dose, cell state and treatment-linked sequencing libraries",
          "None of 96 have yet been subjected to full 9 edgeR leave-one-pool-out refits",
          "No clinical steroid-induced intraocular-pressure outcome data"]}
    out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True)
    with (out/"P5p_96_expression_aware_discordances.tsv").open("w",newline="",encoding="utf8") as f:
        writer=csv.DictWriter(f,fieldnames=list(enriched[0]),delimiter="\t")
        writer.writeheader();writer.writerows(enriched)
    with (out/"P5p_9_expression_adequate_candidates.tsv").open("w",newline="",encoding="utf8") as f:
        writer=csv.DictWriter(f,fieldnames=list(target[0]),delimiter="\t")
        writer.writeheader();writer.writerows(target)
    (out/"P5p_expression_aware_summary.json").write_text(
          json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("P5P_PASS",counts,"NINE",result["adequately_expressed_candidates"],flush=True)
if __name__=="__main__":main()
