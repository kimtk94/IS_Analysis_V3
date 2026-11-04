#!/usr/bin/env python3
"""Read-only three-locus source archive inventory; prints metadata only."""
import csv,json,os
from pathlib import Path
R=Path("/srv/is-analysis")
genes=("UMOD","INHBC","SDCCAG8")
q=R/"results/ckd/stage2/stage2_pqtl_schema.tsv"
with q.open(newline="") as f:
 rows=[r for r in csv.DictReader(f,delimiter="\t") if r.get("gene_symbol","").upper() in genes]
short=[{k:r.get(k) for k in ("gene_symbol","ancestry","archive","member","header")} for r in rows]
archive_files=[]
for root in (R/"data/ckd",R/"results/ckd"):
 for current,ds,names in os.walk(root):
  depth=len(Path(current).relative_to(root).parts)
  ds[:]=[d for d in ds if not d.startswith(".")] if depth<7 else []
  for name in names:
   if not any(g.lower() in name.lower() for g in genes):continue
   if not (name.endswith(".tar") or name.endswith(".tar.gz") or name.endswith(".tsv.gz")):continue
   p=Path(current)/name
   try: size=p.stat().st_size
   except OSError:continue
   if len(archive_files)<80:archive_files.append({"file":str(p),"size_bytes":size})
print(json.dumps({"schema_records":short,"local_gene_files":archive_files},indent=2))
