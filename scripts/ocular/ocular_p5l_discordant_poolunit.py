#!/usr/bin/env python3
"""P5l: cross-check 96 sign-discordant published TM genes vs conservative pool model.

Only aggregate per-gene published results and model statistics. No donor data.
"""
import argparse,csv,gzip,io,json,math
from collections import Counter
from pathlib import Path

def read_tsv(path):
    with open(path,encoding="utf8",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))

def read_gz(path):
    with gzip.open(path,"rt",encoding="utf8") as f:
        return list(csv.DictReader(f,delimiter="\t"))

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--discordant",required=True)
    p.add_argument("--poolunit",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    disc=read_tsv(a.discordant)
    pool=read_gz(a.poolunit)
    if len(disc)!=96:raise ValueError("EXPECTED_96_SIGN_DISCORDANT_GENES")
    poolmap={r["gene_id"]:r for r in pool}
    if len(poolmap)<20000:raise ValueError("EXPECTED_COMPLETE_POOL_LIBRARY_MODEL")
    rows=[]
    for d in disc:
        gid=d["gene_id"];r=poolmap.get(gid)
        if r is None:raise ValueError("DISCORDANT_GENE_ABSENT_POOL_MODEL "+gid)
        p_lfc=float(d["primary_log2fc"])
        e_lfc=float(d["external_TM_log2fc"])
        u_lfc=float(r["pool_logFC"])
        p_fdr=float(d["primary_fdr"])
        e_fdr=float(d["external_TM_fdr"])
        u_fdr=float(r["pool_FDR"])
        if not p_fdr<.05 or not e_fdr<.05 or p_lfc*e_lfc>=0:
            raise ValueError("NOT_BOTH_SIGNIFICANT_OPPOSITE "+gid)
        rows.append({"gene_id":gid,"gene_symbol":d["gene_symbol"],
          "primary_log2fc":p_lfc,"primary_FDR":p_fdr,
          "external_2026_TM_log2fc":e_lfc,"external_2026_TM_FDR":e_fdr,
          "poolunit_log2fc":u_lfc,"poolunit_FDR":u_fdr,
          "poolunit_same_as_primary":p_lfc*u_lfc>0,
          "poolunit_FDR_lt_0_05":u_fdr<.05,
          "primary_abs_log2fc_ge1":abs(p_lfc)>=1,
          "triangulated_conflict":u_fdr<.05 and p_lfc*u_lfc>0})
    confirmed=[x for x in rows if x["triangulated_conflict"]]
    report={"status":"PASS_P5L_GENE_LEVEL_POOLUNIT_CONFLICT_AUDIT",
     "source_P5j_run":41158343447,"source_P5d_run":41090827791,
     "published_both_FDR_sign_opposite":len(rows),
     "pool_model_gene_coverage":len(rows),
     "pool_model_same_direction_as_primary":sum(x["poolunit_same_as_primary"] for x in rows),
     "pool_model_FDR_significant":sum(x["poolunit_FDR_lt_0_05"] for x in rows),
     "pool_model_significant_same_as_primary":len(confirmed),
     "primary_abs_logFC_ge1":sum(x["primary_abs_log2fc_ge1"] for x in rows),
     "primary_abslfc_lt0_25":sum(abs(x["primary_log2fc"])<.25 for x in rows),
     "poolunit_direction_discordance_genes":[x["gene_symbol"] for x in rows if not x["poolunit_same_as_primary"]],
     "top_ten_robust_conflicting_genes":[x["gene_symbol"] for x in confirmed[:10]],
     "warnings":["Neither primary nor pool-level model separates DEX from condition-specific library preparation",
       "External experiment is 100nM bulk TM, primary nominal 500nM BeamA pseudobulk",
       "Cross-study sign conflict is not evidence of novelty, causal biology or clinical steroid-induced IOP"]}
    out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
    with (out/"P5l_96_discordant_poolunit_crosscheck.tsv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        writer.writeheader();writer.writerows(rows)
    (out/"P5l_discordant_poolunit_summary.json").write_text(
       json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    print("P5L_PASS",len(rows),"pool_same",report["pool_model_same_direction_as_primary"],
      "pool_FDR",report["pool_model_FDR_significant"],"both",len(confirmed))
if __name__=="__main__":main()
