#!/usr/bin/env python3
"""Report only column names of already-known CKD aggregate source tables."""
import csv,gzip,json
from pathlib import Path
root=Path("/srv/is-analysis")
paths={
 "pqtl_manifest":root/"results/ckd/stage2/stage2_pqtl_download_manifest.tsv",
 "pqtl_schema":root/"results/ckd/stage2/stage2_pqtl_schema.tsv",
 "eur_egfrcrea":root/"data/ckd/analysis_ready/EUR/eGFRcrea.tsv.gz",
 "ukb_st16":root/"data/ckd/analysis_ready/UKBPPP_ST16_cis_independent.tsv.gz",
}
def header(path):
 if not path.is_file():return {"present":False}
 with (gzip.open(path,"rt") if path.suffix==".gz" else path.open()) as f:
  return {"present":True,"columns":next(csv.reader(f,delimiter="\t"),[])[:100]}
print(json.dumps({k:header(v) for k,v in paths.items()},indent=2))
