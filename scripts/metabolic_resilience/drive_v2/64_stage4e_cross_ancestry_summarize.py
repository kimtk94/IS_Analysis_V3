#!/usr/bin/env python3
"""Summarize EUR-to-EAS replication without over-interpreting non-significance.

Expected input is a reviewed, resource-specific EAS result table assembled
after allele harmonization. This script classifies each tested pair as:
- CONCORDANT_SIGNIFICANT
- CONCORDANT_NONSIGNIFICANT
- OPPOSITE_DIRECTION
- NOT_TESTABLE

The default nominal EAS significance threshold is configurable with
EAS_REPLICATION_ALPHA and is reported in the output.
"""

from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path

ROOT = Path(os.environ.get("IS_ANALYSIS_ROOT", "/srv/is-analysis"))
OUTDIR = ROOT / "results/metabolic_resilience/stage4_eas"
SRC = OUTDIR / "STAGE4E_EAS_HARMONIZED_RESULTS.tsv"
OUT = OUTDIR / "STAGE4E_CROSS_ANCESTRY_SUMMARY.tsv"
SUMMARY = OUTDIR / "STAGE4E_CROSS_ANCESTRY_SUMMARY.json"
TEMPLATE = OUTDIR / "STAGE4E_EAS_HARMONIZED_RESULTS_TEMPLATE.tsv"

OUTDIR.mkdir(parents=True, exist_ok=True)
ALPHA = float(os.environ.get("EAS_REPLICATION_ALPHA", "0.05"))

FIELDS = [
    "gene", "protein", "trait", "domain", "resource", "platform",
    "pqtl_snp_eur", "pqtl_snp_eas", "eaf_eur", "eaf_eas",
    "beta_protein_eur", "beta_protein_eas",
    "mr_beta_eur", "mr_beta_eas", "mr_p_eas",
    "coloc_eur_pph4", "coloc_eas_pph4",
    "sample_overlap_with_pqtl", "status", "notes",
]

if not SRC.exists():
    with TEMPLATE.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
    payload = {
        "status": "HOLD_EAS_HARMONIZED_RESULTS_MISSING",
        "template": str(TEMPLATE),
    }
    SUMMARY.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    raise SystemExit(2)


def fnum(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except Exception:
        return None


with SRC.open("r", encoding="utf-8", newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

out = []
for r in rows:
    eur = fnum(r.get("mr_beta_eur"))
    eas = fnum(r.get("mr_beta_eas"))
    p = fnum(r.get("mr_p_eas"))
    source_status = (r.get("status") or "").strip().upper()

    if source_status in {"NOT_TESTABLE", "ABSENT", "RARE", "NO_PQTL", "NO_OUTCOME"} or eur is None or eas is None:
        cls = "NOT_TESTABLE"
        concordant = ""
    else:
        concordant_bool = (eur > 0 and eas > 0) or (eur < 0 and eas < 0)
        concordant = int(concordant_bool)
        if concordant_bool and p is not None and p < ALPHA:
            cls = "CONCORDANT_SIGNIFICANT"
        elif concordant_bool:
            cls = "CONCORDANT_NONSIGNIFICANT"
        else:
            cls = "OPPOSITE_DIRECTION"

    rr = dict(r)
    rr["direction_concordance"] = concordant
    rr["replication_class"] = cls
    rr["eas_alpha"] = ALPHA
    out.append(rr)

fields = list(FIELDS)
for extra in ["direction_concordance", "replication_class", "eas_alpha"]:
    if extra not in fields:
        fields.append(extra)

with OUT.open("w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    w.writeheader()
    w.writerows(out)

payload = {
    "status": "PASS",
    "rows": len(out),
    "eas_alpha": ALPHA,
    "concordant_significant": sum(r["replication_class"] == "CONCORDANT_SIGNIFICANT" for r in out),
    "concordant_nonsignificant": sum(r["replication_class"] == "CONCORDANT_NONSIGNIFICANT" for r in out),
    "opposite_direction": sum(r["replication_class"] == "OPPOSITE_DIRECTION" for r in out),
    "not_testable": sum(r["replication_class"] == "NOT_TESTABLE" for r in out),
    "interpretation_rule": "EAS non-significance alone is not treated as refutation.",
}
SUMMARY.write_text(json.dumps(payload, indent=2) + "\n")
print(json.dumps(payload, indent=2))
