#!/usr/bin/env python3
"""Read-only targeted Stanzick 2021 eGFRcrea original GWAS lookup."""
import gzip,json,re
from pathlib import Path
path=Path("/srv/is-analysis/data/ckd/rawdata/outcome/EUR/ckdgen_stanzick2021_egfr_eur/metal_eGFR_meta_ea1.TBL.map.annot.gc.gz")
target={"rs35646972","rs6663222","rs11834807","rs11860566"}
positions={"1:242690143","1:242990635","12:58433974","16:20208930"}
found=[];total=0;columns=[]
with gzip.open(path,"rt",errors="replace") as f:
 columns=re.split(r"\\s+",f.readline().strip())
 for line in f:
  if not any(rs in line for rs in target) and not any(pos in line for pos in positions):\n   total+=1\n   continue\n  parts=re.split(r"\\s+",line.strip())
  if len(parts)!=len(columns):continue
  total+=1
  values=dict(zip(columns,parts))
  joined=" ".join(parts)
  if any(re.search(r"(?<![A-Za-z0-9])"+rs+r"(?![A-Za-z0-9])",joined) for rs in target) or any(pos in joined for pos in positions):
   found.append(values)
   if len(found)>20:break
print(json.dumps({"audit":"ORIGINAL_STANZICK_EGFRCREA_TARGETS","columns":columns,"rows_scanned":total,"matches":found},indent=2))
