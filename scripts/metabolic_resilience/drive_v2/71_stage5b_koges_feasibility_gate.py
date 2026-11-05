#!/usr/bin/env python3
"""KoGES longitudinal feasibility gate for Stage 5.

The gate follows the master plan and deliberately requires reviewed metadata
for criteria whose numeric threshold was not prespecified.
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

ROOT = Path(os.environ.get("IS_ANALYSIS_ROOT", "/srv/is-analysis"))
OUTDIR = ROOT / "results/metabolic_resilience/stage5_koges"
ELIG = OUTDIR / "KOGES_LONGITUDINAL_ELIGIBILITY.tsv"
META = OUTDIR / "STAGE5B_KOGES_FEASIBILITY_INPUT.json"
SUMMARY = OUTDIR / "STAGE5B_KOGES_FEASIBILITY_SUMMARY.json"
TEMPLATE = OUTDIR / "STAGE5B_KOGES_FEASIBILITY_INPUT_TEMPLATE.json"

OUTDIR.mkdir(parents=True, exist_ok=True)

if not META.exists():
    template = {
        "metabolic_variables_usable_most_waves": None,
        "physical_activity_usable_waves": None,
        "baseline_healthy_genetic_n_sufficient": None,
        "genotype_pqtl_allele_harmonization_possible": None,
        "medication_variables_complete": None,
        "notes": "Fill after participant-level/raw-data audit. Do not infer unspecified sufficiency thresholds.",
    }
    TEMPLATE.write_text(json.dumps(template, indent=2) + "\n")
    payload = {
        "status": "HOLD_REVIEWED_FEASIBILITY_METADATA_MISSING",
        "template": str(TEMPLATE),
    }
    SUMMARY.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    raise SystemExit(2)

meta = json.loads(META.read_text())

required_meta = [
    "metabolic_variables_usable_most_waves",
    "physical_activity_usable_waves",
    "baseline_healthy_genetic_n_sufficient",
    "genotype_pqtl_allele_harmonization_possible",
    "medication_variables_complete",
]
missing = [k for k in required_meta if meta.get(k) is None]
if missing:
    payload = {"status": "HOLD_METADATA_INCOMPLETE", "missing": missing}
    SUMMARY.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    raise SystemExit(2)

if not ELIG.exists():
    payload = {
        "status": "HOLD_ELIGIBILITY_TABLE_MISSING",
        "expected": str(ELIG),
        "required_columns": [
            "participant_id", "n_total_visits", "n_valid_metabolic_visits",
            "baseline_healthy", "incident_mets", "incident_t2d",
            "followup_years", "genotype_available",
        ],
    }
    SUMMARY.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    raise SystemExit(2)

with ELIG.open("r", encoding="utf-8", newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

required_cols = {
    "participant_id", "n_total_visits", "n_valid_metabolic_visits",
    "baseline_healthy", "incident_mets", "incident_t2d",
    "followup_years", "genotype_available",
}
if not rows:
    raise SystemExit("Eligibility table is empty")
missing_cols = required_cols - set(rows[0])
if missing_cols:
    raise SystemExit(f"Missing eligibility columns: {sorted(missing_cols)}")


def truth(x):
    return str(x).strip().lower() in {"1", "true", "yes", "y"}


def integer(x):
    try:
        return int(float(x))
    except Exception:
        return 0


longitudinal = [r for r in rows if integer(r["n_valid_metabolic_visits"]) >= 2]
genetic_longitudinal = [r for r in longitudinal if truth(r["genotype_available"])]
baseline_healthy_genetic = [
    r for r in genetic_longitudinal if truth(r["baseline_healthy"])
]

n_long = len(longitudinal)
n_gen_long = len(genetic_longitudinal)
n_base_healthy_gen = len(baseline_healthy_genetic)

go_core = (
    int(meta["metabolic_variables_usable_most_waves"]) >= 4
    and n_gen_long >= 4000
    and bool(meta["baseline_healthy_genetic_n_sufficient"])
    and bool(meta["genotype_pqtl_allele_harmonization_possible"])
)

modifications = []
if int(meta["physical_activity_usable_waves"]) < 3:
    modifications.append("USE_BASELINE_PA_STRATIFICATION")
if not bool(meta["medication_variables_complete"]):
    modifications.append("STRENGTHEN_LAB_ONLY_CONTINUOUS_MBI")

if not go_core:
    decision = "NO_GO_OR_ROLE_SPLIT"
elif modifications:
    decision = "MODIFY"
else:
    decision = "GO"

payload = {
    "status": decision,
    "participants_total": len(rows),
    "longitudinal_n_ge2_metabolic_visits": n_long,
    "genotype_longitudinal_n": n_gen_long,
    "baseline_healthy_genetic_n": n_base_healthy_gen,
    "incident_mets_genetic_longitudinal_n": sum(
        truth(r["incident_mets"]) for r in genetic_longitudinal
    ),
    "incident_t2d_genetic_longitudinal_n": sum(
        truth(r["incident_t2d"]) for r in genetic_longitudinal
    ),
    "metabolic_variables_usable_most_waves": meta["metabolic_variables_usable_most_waves"],
    "physical_activity_usable_waves": meta["physical_activity_usable_waves"],
    "baseline_healthy_genetic_n_sufficient_reviewed": meta["baseline_healthy_genetic_n_sufficient"],
    "genotype_pqtl_allele_harmonization_possible": meta["genotype_pqtl_allele_harmonization_possible"],
    "medication_variables_complete": meta["medication_variables_complete"],
    "modifications": modifications,
    "alternative_if_no_go": "KoGES longitudinal observational validation + EAS summary-MR role split",
}
SUMMARY.write_text(json.dumps(payload, indent=2) + "\n")
print(json.dumps(payload, indent=2))
