#!/usr/bin/env python3
"""P5g rank statistically robust gene-level aggregates from verified P5b/P5d ZIPs.

Does not access donor-level expression. Unreviewed DOES NOT mean novel.
Usage: python scripts/ocular/ocular_p5g_candidate_screen.py \
  --p5b <downloaded-p5b-artifact.zip> --p5d <downloaded-p5d-artifact.zip> \
  --output <folder>
"""
import argparse,csv,gzip,io,json,math,zipfile
from collections import defaultdict
from pathlib import Path
KNOWN={
 "HIPK2":"Singh et al 2026 chromatin preprint",
 "PTK2B":"Kathirvel et al 2022 DEX RNA-seq/qPCR",
 "FKBP5":"Known glucocorticoid response",
 "FAM107A":"Mehrotra et al 2026 TM",
 "ZBTB16":"Singh 2026 and Kathirvel 2022",
 "ATP10A":"Singh 2026",
 "EDIL3":"Singh 2026",
 "MAOA":"2026 steroid TM fibrosis MAOA-MIEF1-DRP1",
 "MIEF1":"2026 steroid TM fibrosis pathway",
 "DNM1L":"2026 steroid TM fibrosis pathway",
 "MYOC":"Known TM DEX/glaucoma marker",
 "ANGPTL7":"Known TM DEX marker",
 "TSC22D3":"Glucocorticoid response marker",
 "SAA1":"Kathirvel et al Scientific Reports 2022",
 "SAA2":"Kathirvel et al Scientific Reports 2022",
 "OCA2":"Kathirvel et al Scientific Reports 2022",
 "BDKRB1":"2022 Int J Ophthalmol DEX TM qPCR",
}
def zrows(zpath,filename):
    with zipfile.ZipFile(zpath) as z:
        with z.open(filename) as h:
            if filename.endswith(".gz"):
                with gzip.GzipFile(fileobj=h) as f:
                    return list(csv.DictReader(io.TextIOWrapper(f,encoding="utf8"),delimiter="\t"))
            return list(csv.DictReader(io.TextIOWrapper(h,encoding="utf8"),delimiter="\t"))
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--p5b",required=True)
    p.add_argument("--p5d",required=True)
    p.add_argument("--output",required=True)
    a=p.parse_args()
    base=zrows(a.p5b,"p5_all_gene_results.tsv.gz")
    pool=zrows(a.p5d,"p5d_poolunit_all_gene_results.tsv.gz")
    loo=zrows(a.p5b,"p5b_top_gene_leave_pool.tsv")
    sd=zrows(a.p5d,"p5d_top100_pool_consistency.tsv")
    bypool={x["gene_id"]:x for x in pool}
    bysd={x["gene_id"]:x for x in sd}
    byloo=defaultdict(list)
    for x in loo:byloo[x["gene_id"]].append(x)
    results=[]
    for x in base:
        gene=x["gene_id"];y=bypool.get(gene);z=bysd.get(gene);folds=byloo.get(gene,[])
        if y is None or z is None or len(folds)!=9:continue
        effect=float(x["logFC"])
        pool_effect=float(y["pool_logFC"])
        lfc=[float(v["sensitivity_logFC"]) for v in folds]
        lfdr=[float(v["sensitivity_FDR"]) for v in folds]
        per_pool_cv=float(z["pool_delta_sd"])/max(abs(float(z["pooled_mean_delta"])),0.01)
        if not (float(x["FDR"])<0.05 and float(y["pool_FDR"])<0.05
          and abs(effect)>=1.0 and float(x["logCPM"])>=2.0
          and min(map(abs,lfc))>=0.8 and max(lfdr)<0.05
          and sum(math.copysign(1,v)==math.copysign(1,effect) for v in lfc)==9
          and int(z["pools_agree"])==9
          and math.copysign(1,pool_effect)==math.copysign(1,effect)
          and per_pool_cv<=.35):continue
        symbol=x["gene_symbol"]
        results.append(dict(gene_symbol=symbol,gene_id=gene,
           primary_log2fc=round(effect,5),primary_fdr=float(x["FDR"]),
           primary_logcpm=round(float(x["logCPM"]),4),
           poolunit_log2fc=round(pool_effect,5),
           poolunit_fdr=float(y["pool_FDR"]),
           loo_abs_log2fc_min=round(min(map(abs,lfc)),5),
           loo_fdr_max=max(lfdr),loo_sign_agree=9,
           poolunit_9pool_sign_agree=9,
           pool_delta_sd=round(float(z["pool_delta_sd"]),5),
           pool_cv=round(per_pool_cv,5),
           literature_flag="KNOWN_RESPONSE" if symbol in KNOWN else "UNREVIEWED_NOT_NOVELTY_CONFIRMED",
           known_reference=KNOWN.get(symbol,"")))
    results.sort(key=lambda x:(x["literature_flag"]=="KNOWN_RESPONSE",
                                  -abs(x["primary_log2fc"]),x["gene_id"]))
    out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
    if not results:raise ValueError("NO_CANDIDATES_PASS_PRESPECIFIED_AGGREGATE_FILTERS")
    with (out/"P5g_robust_candidates.tsv").open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=results[0].keys(),delimiter="\t")
        w.writeheader();w.writerows(results)
    report={"status":"PASS_AGGREGATE_SCREEN_ONLY",
       "input_primary_run":41085982185,"input_poolunit_run":41090827791,
       "primary_tested_genes":len(base),"poolunit_tested_genes":len(pool),
       "screen_scope":"top100 primary genes with complete 9-pool summary",
       "screen_n":len(results),
       "known_count":sum(x["literature_flag"]=="KNOWN_RESPONSE" for x in results),
       "unreviewed_not_novelty_count":sum(x["literature_flag"]!="KNOWN_RESPONSE" for x in results),
       "criteria":{"both_FDR_lt":0.05,"abs_primary_log2fc_ge":1,"primary_logcpm_ge":2,
         "min_LOO_abs_log2fc_ge":0.8,"max_LOO_FDR_lt":0.05,
         "all_9_sign_match":True,"max_pool_sd_div_abs_mean":0.35},
       "warnings":["No causal novelty claim","Only top100 screened",
                   "Known lookup incomplete; unreviewed != novel",
                   "Technical library differences cannot be separated from DEX treatment"]}
    (out/"P5g_candidate_screen_summary.json").write_text(json.dumps(report,indent=2)+"\n")
    print("P5G_PASS",report["screen_n"],"known",report["known_count"],
          "unreviewed",report["unreviewed_not_novelty_count"])
if __name__=="__main__":main()
