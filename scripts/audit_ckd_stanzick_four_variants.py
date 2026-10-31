#!/usr/bin/env python3
"""Read-only Stanzick eGFRcrea original GWAS four-target lookup."""
import gzip,json
from pathlib import Path
p=Path("/srv/is-analysis/data/ckd/rawdata/outcome/EUR/ckdgen_stanzick2021_egfr_eur/metal_eGFR_meta_ea1.TBL.map.annot.gc.gz")
target=("rs35646972","rs6663222","rs11834807","rs11860566")
coords=("1:242690143","1:242990635","12:58433974","16:20208930")
rows=0
matches=[]
with gzip.open(p,"rt",errors="replace") as fh:
    headers=fh.readline().split()
    for line in fh:
        rows+=1
        if not any(t in line for t in target) and not any(t in line for t in coords):
            continue
        cols=line.split()
        if len(cols)==len(headers):
            matches.append(dict(zip(headers,cols)))
print(json.dumps({"audit":"ORIGINAL_STANZICK_EGFRCREA_TARGETS","columns":headers,"rows_scanned":rows,"matches":matches[:30]},indent=2))
