#!/usr/bin/env python3
"""Check four CKD MAF-outlier rsIDs in existing analysis-ready aggregate tables."""
import csv,gzip,json
from pathlib import Path
ROOT=Path("/srv/is-analysis/data/ckd/analysis_ready")
TARGET=("rs35646972","rs6663222","rs11834807","rs11860566")
FILES={"eur_egfrcrea":ROOT/"EUR/eGFRcrea.tsv.gz","ukb_st16":ROOT/"UKBPPP_ST16_cis_independent.tsv.gz"}
def lookup(path):
 with gzip.open(path,"rt",newline="") as f:
  reader=csv.DictReader(f,delimiter="\t")
  rows={r["rsid"]:r for r in reader if r.get("rsid") in TARGET}
  valid=("rsid","gene_symbol","variant_id","effect_allele","other_allele","ref_allele","alt_allele","eaf","alt_freq","beta","se","effect_scale","source_schema")
  return [{k:r.get(k) for k in valid if k in r} if (r:=rows.get(rsid)) else {"rsid":rsid,"status":"ABSENT_IN_ANALYSIS_READY"} for rsid in TARGET]
print(json.dumps({"type":"CKD_FOUR_OUTLIER_AGGREGATE_SOURCE_LOOKUP","results":{k:lookup(p) for k,p in FILES.items()},"not_original_full_locus_consortrium":True},indent=2))
