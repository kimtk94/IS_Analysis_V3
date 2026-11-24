#!/usr/bin/env python3
"""Read-only audit of original pQTL coordinate conventions and GWAS provenance."""
import csv,gzip,json,os
from pathlib import Path
root=Path("/srv/is-analysis")
p=root/"data/ckd/analysis_ready/EUR/eGFRcrea.summary.json"
result={"egfrcrea_summary":None,"source_hints":[]}
if p.is_file():
 try:result["egfrcrea_summary"]=json.loads(p.read_text())
 except Exception as e:result["egfrcrea_summary_error"]=type(e).__name__
for base in (root/"results/ckd/stage2",root/"data/ckd/rawdata",root/"results/ckd/stage1"):
 if not base.exists():continue
 for d,dirs,names in os.walk(base):
  depth=len(Path(d).relative_to(base).parts)
  dirs[:]=dirs if depth<5 else []
  for name in names:
   if any(t in name.lower() for t in ("manifest","source","provenance")) and name.endswith((".json",".tsv",".txt")):
    if len(result["source_hints"])<60:result["source_hints"].append(str(Path(d)/name))
print(json.dumps(result,indent=2))
