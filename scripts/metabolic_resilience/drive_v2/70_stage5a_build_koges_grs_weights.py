#!/usr/bin/env python3
from pathlib import Path
import csv
import json

ROOT=Path("/srv/is-analysis")
INST=ROOT/"results/metabolic_resilience/stage3_full_pgwas/instruments/STAGE3C2_FROZEN_INSTRUMENTS_ALL.tsv"
OUT=ROOT/"results/metabolic_resilience/stage5_koges/grs_weights"
OUT.mkdir(parents=True,exist_ok=True)

rows=list(csv.DictReader(INST.open("r",encoding="utf-8"),delimiter="\t"))
primary=[r for r in rows if r["mode"]=="primary"]

summary=[]
for gene in sorted(set(r["gene_symbol"] for r in primary)):
    x=[r for r in primary if r["gene_symbol"]==gene]
    dst=OUT/f"{gene}.protein_cis_pqtl_grs_weights.tsv"
    fields=["gene","rsid","chrom_hg19","pos_hg19","effect_allele","other_allele","weight_beta_protein","F"]
    outrows=[]
    for r in x:
        outrows.append({
            "gene":gene,
            "rsid":r.get("rsid",""),
            "chrom_hg19":r["chrom_hg19"],
            "pos_hg19":r["pos_hg19"],
            "effect_allele":r["exposure_effect_allele"],
            "other_allele":r["exposure_other_allele"],
            "weight_beta_protein":r["exposure_beta"],
            "F":r["exposure_F"],
        })
    with dst.open("w",encoding="utf-8",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n")
        wr.writeheader();wr.writerows(outrows)
    summary.append({"gene":gene,"k":len(outrows),"file":str(dst)})

(OUT.parent/"STAGE5A_GRS_WEIGHT_SUMMARY.json").write_text(
    json.dumps({"genes":len(summary),"sets":summary},indent=2)+"\n"
)
print(json.dumps({"genes":len(summary),"sets":summary},indent=2))
