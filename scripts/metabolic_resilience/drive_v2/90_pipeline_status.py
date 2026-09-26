#!/usr/bin/env python3
from pathlib import Path
import os

ROOT = Path(os.environ.get("IS_ANALYSIS_ROOT", "/srv/is-analysis"))

checks = [
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
    ("Stage3D5 MR summary", ROOT/"results/metabolic_resilience/stage3_confirmatory_mr/STAGE3D5_PRIMARY_MR_SUMMARY.tsv"),
    ("Stage3E2 coloc", ROOT/"results/metabolic_resilience/stage3_coloc/STAGE3E2_COLOC_ABF.tsv"),
    ("Stage3E3 coloc summary", ROOT/"results/metabolic_resilience/stage3_coloc/STAGE3E3_COLOC_SUMMARY.tsv"),
    ("Stage4A CKB coverage", ROOT/"results/metabolic_resilience/stage4_eas/STAGE4A_CKB_COVERAGE.json"),
    ("Stage4C EUR candidate gate", ROOT/"results/metabolic_resilience/stage4_eas/STAGE4C_EUR_CANDIDATE_GATE.json"),
    ("Stage4D EAS registry", ROOT/"results/metabolic_resilience/stage4_eas/STAGE4D_EAS_RESOURCE_REGISTRY_SUMMARY.json"),
    ("Stage4E cross-ancestry", ROOT/"results/metabolic_resilience/stage4_eas/STAGE4E_CROSS_ANCESTRY_SUMMARY.json"),
    ("Stage5A KoGES GRS weights", ROOT/"results/metabolic_resilience/stage5_koges/STAGE5A_GRS_WEIGHT_SUMMARY.json"),
    ("Stage5B KoGES feasibility", ROOT/"results/metabolic_resilience/stage5_koges/STAGE5B_KOGES_FEASIBILITY_SUMMARY.json"),
    ("Stage5C GRS coverage", ROOT/"results/metabolic_resilience/stage5_koges/STAGE5C_GRS_VARIANT_COVERAGE.json"),
    ("Stage5D analysis contract", ROOT/"results/metabolic_resilience/stage5_koges/STAGE5D_ANALYSIS_CONTRACT.json"),
    ("Stage5E incident MetS Cox", ROOT/"results/metabolic_resilience/stage5_koges/models/STAGE5E_INCIDENT_METS_COX_SUMMARY.json"),
    ("Stage5F repeated MBI LMM", ROOT/"results/metabolic_resilience/stage5_koges/models/STAGE5F_REPEATED_MBI_LMM_SUMMARY.json"),
    ("Stage5G PA interaction", ROOT/"results/metabolic_resilience/stage5_koges/models/STAGE5G_PA_INTERACTION_SUMMARY.json"),
]

print("="*92)
print("METABOLIC RESILIENCE PIPELINE STATUS")
print("="*92)
for name, p in checks:
    state = "PASS/EXISTS" if p.exists() and p.stat().st_size > 0 else "PENDING"
    print(f"{state:12s}  {name:36s}  {p}")
