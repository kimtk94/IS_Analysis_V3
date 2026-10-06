#!/usr/bin/env python3
"""
Stage 4-B scaffold.

Purpose:
- consume confirmed EUR candidates after Stage 3 MR + coloc;
- map to CKB/KoGES/BBJ/TWB EAS resources;
- keep platform/source provenance explicit;
- never count EAS replication as a new metabolic domain.

This script intentionally stops unless a reviewed EAS registry is supplied.
"""
from pathlib import Path
import csv
import json
import sys

ROOT=Path("/srv/is-analysis")
REG=ROOT/"results/metabolic_resilience/stage4_eas/STAGE4B_EAS_REGISTRY.tsv"
OUT=ROOT/"results/metabolic_resilience/stage4_eas"
OUT.mkdir(parents=True,exist_ok=True)

if not REG.exists():
    template=OUT/"STAGE4B_EAS_REGISTRY_TEMPLATE.tsv"
    fields=[
        "resource","ancestry","platform","trait","path","build",
        "protein_key_column","variant_key_column","beta_column","se_column",
        "effect_allele_column","other_allele_column","eaf_column","status"
    ]
    with template.open("w",encoding="utf-8",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n")
        wr.writeheader()
        wr.writerows([
            {"resource":"CKB","ancestry":"EAS","platform":"SomaScan v4.1","trait":"REVIEW",
             "path":str(ROOT/"data/metabolic_resilience/stage1_pqtl/ckb/CKB_SomaScan_MR_coloc.gz"),
             "status":"REVIEW_SCHEMA"},
            {"resource":"KoGES","ancestry":"Korean/EAS","platform":"GWAS","trait":"REVIEW",
             "path":str(ROOT/"data/metabolic_resilience/stage1_pqtl/koges"),"status":"REVIEW_SCHEMA"},
        ])
    print("[HOLD] Review registry:",template)
    sys.exit(2)

rows=list(csv.DictReader(REG.open("r",encoding="utf-8"),delimiter="\t"))
bad=[r for r in rows if r["status"]!="PASS"]
if bad:
    print("[HOLD] EAS registry contains non-PASS rows:",len(bad))
    sys.exit(2)

print(json.dumps({
    "registry_rows":len(rows),
    "status":"READY_FOR_RESOURCE_SPECIFIC_REPLICATION_CODE",
},indent=2))
