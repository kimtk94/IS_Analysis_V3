#!/usr/bin/env python3
"""Build and validate the East-Asian replication resource registry.

Resource-specific schemas are never guessed. The pipeline remains HOLD until
all rows selected for analysis are explicitly reviewed and marked PASS.
Sample overlap with the pQTL source must be recorded.
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

ROOT = Path(os.environ.get("IS_ANALYSIS_ROOT", "/srv/is-analysis"))
OUTDIR = ROOT / "results/metabolic_resilience/stage4_eas"
REG = OUTDIR / "STAGE4D_EAS_RESOURCE_REGISTRY.tsv"
TEMPLATE = OUTDIR / "STAGE4D_EAS_RESOURCE_REGISTRY_TEMPLATE.tsv"
SUMMARY = OUTDIR / "STAGE4D_EAS_RESOURCE_REGISTRY_SUMMARY.json"

OUTDIR.mkdir(parents=True, exist_ok=True)

FIELDS = [
    "resource",
    "ancestry",
    "role",
    "platform",
    "trait",
    "domain",
    "path",
    "member",
    "build",
    "sample_n",
    "case_n",
    "control_n",
    "protein_key_column",
    "variant_key_column",
    "chr_column",
    "pos_column",
    "effect_allele_column",
    "other_allele_column",
    "beta_column",
    "se_column",
    "p_column",
    "eaf_column",
    "sample_overlap_with_pqtl",
    "sample_overlap_note",
    "status",
]

if not REG.exists():
    template_rows = [
        {
            "resource": "CKB",
            "ancestry": "EAS",
            "role": "pQTL",
            "platform": "SomaScan",
            "trait": "PROTEIN",
            "domain": "protein",
            "path": str(ROOT / "data/metabolic_resilience/stage1_pqtl/ckb/CKB_SomaScan_MR_coloc.gz"),
            "build": "REVIEW",
            "sample_overlap_with_pqtl": "SELF",
            "sample_overlap_note": "CKB pQTL source; review before combining with any CKB outcome.",
            "status": "REVIEW_SCHEMA",
        },
        {
            "resource": "KoGES",
            "ancestry": "Korean/EAS",
            "role": "outcome_GWAS",
            "platform": "KoreanChip+imputation",
            "trait": "REVIEW",
            "domain": "REVIEW",
            "path": str(ROOT / "data/metabolic_resilience/stage2_gwas/eas/koges"),
            "build": "GRCh37",
            "sample_overlap_with_pqtl": "REVIEW",
            "status": "REVIEW_SCHEMA",
        },
        {
            "resource": "BBJ",
            "ancestry": "Japanese/EAS",
            "role": "outcome_GWAS",
            "platform": "GWAS",
            "trait": "REVIEW",
            "domain": "REVIEW",
            "path": str(ROOT / "data/metabolic_resilience/stage2_gwas/eas/bbj"),
            "build": "REVIEW",
            "sample_overlap_with_pqtl": "REVIEW",
            "status": "REVIEW_SCHEMA",
        },
        {
            "resource": "TWB",
            "ancestry": "Taiwanese/EAS",
            "role": "outcome_GWAS",
            "platform": "GWAS",
            "trait": "REVIEW",
            "domain": "REVIEW",
            "path": str(ROOT / "data/metabolic_resilience/stage2_gwas/eas/twb"),
            "build": "REVIEW",
            "sample_overlap_with_pqtl": "REVIEW",
            "status": "REVIEW_SCHEMA",
        },
    ]
    with TEMPLATE.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(template_rows)

    payload = {
        "status": "HOLD_REGISTRY_MISSING",
        "template": str(TEMPLATE),
        "instruction": "Review exact file schemas and sample overlap, then save as STAGE4D_EAS_RESOURCE_REGISTRY.tsv.",
    }
    SUMMARY.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    raise SystemExit(2)

with REG.open("r", encoding="utf-8", newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

required_common = [
    "resource", "ancestry", "role", "trait", "path", "build",
    "sample_overlap_with_pqtl", "sample_overlap_note", "status",
]

problems = []
for i, r in enumerate(rows, start=2):
    for col in required_common:
        if not str(r.get(col, "")).strip():
            problems.append(f"line {i}: missing {col}")

    if r.get("status") == "PASS":
        role = r.get("role", "")
        if role == "pQTL":
            for col in ["protein_key_column", "variant_key_column", "effect_allele_column", "other_allele_column", "beta_column", "se_column"]:
                if not str(r.get(col, "")).strip():
                    problems.append(f"line {i}: PASS pQTL row missing {col}")
        elif role == "outcome_GWAS":
            for col in ["variant_key_column", "effect_allele_column", "other_allele_column", "beta_column", "se_column", "p_column"]:
                if not str(r.get(col, "")).strip():
                    problems.append(f"line {i}: PASS outcome row missing {col}")
        else:
            problems.append(f"line {i}: unknown role={role!r}")

        p = Path(r["path"])
        if not p.exists():
            problems.append(f"line {i}: PASS row path does not exist: {p}")

nonpass = [r for r in rows if r.get("status") != "PASS"]

if problems or nonpass:
    payload = {
        "status": "HOLD_REVIEW_REQUIRED",
        "rows": len(rows),
        "pass_rows": len(rows) - len(nonpass),
        "nonpass_rows": len(nonpass),
        "problems": problems,
    }
    SUMMARY.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    raise SystemExit(2)

payload = {
    "status": "PASS",
    "rows": len(rows),
    "resources": sorted({r["resource"] for r in rows}),
    "traits": sorted({r["trait"] for r in rows}),
    "sample_overlap_explicit": all(r.get("sample_overlap_with_pqtl") for r in rows),
}
SUMMARY.write_text(json.dumps(payload, indent=2) + "\n")
print(json.dumps(payload, indent=2))
