#!/usr/bin/env python3
"""Read-only inventory of stage2 pQTL source manifest and eGFRcrea source candidates.

Prints only study-level metadata, no raw genotypes or per-participant records.
"""
import csv,json,os
from pathlib import Path
B=Path("/srv/is-analysis")
manifest=B/"results/ckd/stage2/stage2_pqtl_download_manifest.tsv"
genes={"SDCCAG8","INHBC","UMOD"}
with manifest.open(newline="") as f:
    rec=[r for r in csv.DictReader(f,delimiter="\t") if r.get("gene_symbol","").upper() in genes]
safe=("gene_symbol","ancestry","stage2_role","source_file","synapse_id","expected_size_bytes","md5","sha256")
out={"pqtl_manifest_records":[{k:r.get(k) for k in safe} for r in rec],"gwas_candidates":[]}
for root in (B/"data/ckd",B/"results/ckd"):
    for current,dirs,files in os.walk(root):
        depth=len(Path(current).relative_to(root).parts)
        dirs[:]=[d for d in dirs if not d.startswith(".")] if depth<7 else []
        for n in files:
            l=n.lower()
            if not (("egfr" in l or "ckdgen" in l) and ("crea" in l or "gwas" in l)):continue
            p=Path(current)/n
            try:size=p.stat().st_size
            except OSError:continue
            if len(out["gwas_candidates"])<70:
                out["gwas_candidates"].append({"path":str(p),"bytes":size})
print(json.dumps(out,indent=2))
