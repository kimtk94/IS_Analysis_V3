#!/usr/bin/env python3
from pathlib import Path
import json

ROOT=Path("/srv/is-analysis")

checks=[
    ("Stage2 final shortlist", ROOT/"results/metabolic_resilience/stage2_gwas/final_shortlist/STAGE2_FINAL_SCREENING_SHORTLIST.tsv"),
    ("Stage3A downloads", ROOT/"results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3A_DOWNLOAD_MANIFEST.tsv"),
    ("Stage3B1 coordinate lock", ROOT/"results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3B1_COORDINATE_SUMMARY.json"),
    ("Stage3B2 cis windows", ROOT/"results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3B2_GENE_CIS_WINDOWS.tsv"),
    ("Stage3B3 cis marginal", ROOT/"results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3B3_CIS_EXTRACTION_SUMMARY.json"),
    ("Stage3C0R LD reference", ROOT/"results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3C0R_LD_MAPPING_SUMMARY.json"),
    ("Stage3C1 clumping", ROOT/"results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3C1_LD_CLUMP_SUMMARY.json"),
    ("Stage3C2 frozen instruments", ROOT/"results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3C2_INSTRUMENT_FREEZE.json"),
    ("Stage3D0 outcome registry", ROOT/"results/metabolic_resilience/stage3_confirmatory_mr/STAGE3D0_OUTCOME_REGISTRY.tsv"),
    ("Stage3D3 MR", ROOT/"results/metabolic_resilience/stage3_confirmatory_mr/STAGE3D3_MR_RESULTS.tsv"),
    ("Stage3E2 coloc", ROOT/"results/metabolic_resilience/stage3_coloc/STAGE3E2_COLOC_ABF.tsv"),
    ("Stage4A CKB coverage", ROOT/"results/metabolic_resilience/stage4_eas/STAGE4A_CKB_COVERAGE.json"),
    ("Stage5A KoGES GRS weights", ROOT/"results/metabolic_resilience/stage5_koges/STAGE5A_GRS_WEIGHT_SUMMARY.json"),
]

print("="*86)
print("METABOLIC RESILIENCE PIPELINE STATUS")
print("="*86)
for name,p in checks:
    state="PASS/EXISTS" if p.exists() and p.stat().st_size>0 else "PENDING"
    print(f"{state:12s}  {name:34s}  {p}")
