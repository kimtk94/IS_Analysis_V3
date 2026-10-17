#!/usr/bin/env python3
"""Write the locked KoGES longitudinal analysis contract.

This stage records, rather than invents, the prespecified analysis hierarchy.
It is intentionally data-independent and can be reviewed before model fitting.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

ROOT = Path(os.environ.get("IS_ANALYSIS_ROOT", "/srv/is-analysis"))
OUTDIR = ROOT / "results/metabolic_resilience/stage5_koges"
OUT = OUTDIR / "STAGE5D_ANALYSIS_CONTRACT.json"
OUTDIR.mkdir(parents=True, exist_ok=True)

contract = {
    "primary": {
        "name": "Baseline healthy to incident MetS",
        "baseline_definition": "MetS component count <= 1",
        "event_definition": "incident MetS: component count >= 3",
        "model": "Cox proportional hazards",
        "genetic_predictor": "external cis-pQTL ProteinGRS",
        "base_covariates": ["age", "sex", "genetic_PCs"],
        "note": "Do not over-adjust the primary genetic model with lifestyle covariates.",
    },
    "secondary": {
        "name": "Repeated continuous Metabolic Burden Index",
        "model": "linear mixed model",
        "formula": "MBI ~ Time + ProteinGRS + ProteinGRS:Time + (1 + Time | participant)",
        "target_term": "ProteinGRS:Time",
        "mbi_components": ["waist", "SBP", "fasting_glucose", "log_TG", "-HDL"],
        "sex_stratified_standardization": "consider for variables with strong sex differences",
    },
    "supportive": {
        "name": "EUR MR + coloc to EAS replication",
        "interpretation": [
            "concordant + significant = strong replication",
            "concordant + non-significant = potentially underpowered",
            "opposite direction = potential heterogeneity",
            "EUR pQTL absent/rare in EAS = not testable",
        ],
    },
    "exploratory": {
        "trajectory_model": "multinomial logistic regression: TrajectoryClass ~ ProteinGRS + Age + Sex + PCs",
        "primary_modifier": "physical_activity",
        "interaction": "MBI ~ ... + ProteinGRS:PA + ProteinGRS:PA:Time",
        "secondary_modifiers": ["sleep"],
        "other_exploratory_modifiers": ["diet", "alcohol", "smoking"],
    },
    "execution_gate": "Do not advance automatically when an upstream stage is REVIEW/HOLD/WARN/FAIL.",
}

OUT.write_text(json.dumps(contract, indent=2) + "\n")
print(json.dumps(contract, indent=2))
