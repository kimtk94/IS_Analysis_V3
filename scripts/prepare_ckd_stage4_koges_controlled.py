#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path("/srv/is-analysis")
RES = ROOT / "results/ckd"

ANCHOR = RES / "stage4_koges/STAGE4_KOGES_ANCHOR_PANEL.tsv"
AVAIL = (
    RES
    / "stage4_koges/genotype_availability/"
    "KOGES_GENOTYPE_AVAILABILITY.json"
)

OUT = RES / "stage4_koges/controlled_validation"
OUT.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# 1. Anchor panel
# ------------------------------------------------------------

with ANCHOR.open() as f:
    anchors = list(csv.DictReader(f, delimiter="\t"))

if len(anchors) != 9:
    raise RuntimeError(
        f"Expected 9 anchors, found {len(anchors)}"
    )


# ------------------------------------------------------------
# 2. Current genotype availability
# ------------------------------------------------------------

if AVAIL.exists():
    availability = json.loads(AVAIL.read_text())
else:
    availability = {
        "status": "UNKNOWN",
    }


# ------------------------------------------------------------
# 3. Controlled-data analysis contract
# ------------------------------------------------------------

contract = {
    "project": "CKD proteogenomic thesis",
    "stage": "KoGES controlled-access longitudinal validation",

    "current_status": (
        "READY_PENDING_CONTROLLED_DATA"
        if availability.get("status")
        == "CONTROLLED_KOGES_GENOTYPE_NOT_PRESENT_LOCALLY"
        else "CONTROLLED_DATA_REQUIRES_REAUDIT"
    ),

    "required_genetic_inputs": {
        "accepted_formats": [
            "PLINK1 bed/bim/fam",
            "PLINK2 pgen/pvar/psam",
            "BGEN/sample",
            "VCF/VCF.GZ",
        ],
        "required_build_metadata": True,
        "preferred_build": "GRCh37/hg19",
        "required_variant_fields": [
            "chromosome",
            "position",
            "variant ID or rsID",
            "reference allele",
            "alternate allele",
            "dosage or genotype",
        ],
    },

    "required_phenotype_inputs": {
        "participant_id": True,
        "sex": True,
        "age": True,
        "creatinine": True,
        "visit_time_or_date": True,
        "baseline_definition": True,
    },

    "preferred_covariates": [
        "genetic principal components",
        "genotyping batch",
        "study site",
        "age",
        "sex",
    ],

    "primary_endpoints": {
        "longitudinal_egfr": {
            "model": "linear mixed-effects model",
            "primary_term": "genotype × time",
            "minimum_measurements": 3,
            "minimum_followup_years": 2,
        },

        "incident_ckd": {
            "baseline_requirement": "baseline eGFR >= 60",
            "event_definition": (
                "two consecutive later eGFR values <60 "
                "separated by >=90 days"
            ),
        },
    },

    "anchor_strategy": {
        "n_primary_anchors": 9,
        "direct_anchor_first": True,
        "missing_anchor_action": (
            "search ancestry-matched EAS/Korean LD proxy"
        ),
        "proxy_documentation_required": [
            "target rsID",
            "proxy rsID",
            "r2",
            "LD population",
            "allele orientation",
            "genome build",
        ],
    },

    "effect_orientation": (
        "All KoGES effects must ultimately be aligned to "
        "the UKB-PPP protein-increasing/effect allele."
    ),

    "public_training_policy": (
        "Public KoGES training phenotype results are workflow QA only "
        "and must not be reported as manuscript-level genetic inference."
    ),
}


(OUT / "KOGES_CONTROLLED_ANALYSIS_CONTRACT.json").write_text(
    json.dumps(contract, indent=2)
)


# ------------------------------------------------------------
# 4. Anchor extraction specification
# ------------------------------------------------------------

fields = [
    "gene_symbol",
    "rsid",
    "variant_id",
    "chrom_hg19",
    "pos_hg19",
    "ref_allele",
    "effect_allele_alt",
    "beta_cond_alt",
    "se_cond",
    "alt_freq_eur",
    "pip",
    "f_stat",
    "stage2_class",
]

with (
    OUT / "KOGES_CONTROLLED_ANCHOR_SPEC.tsv"
).open("w", newline="") as f:

    w = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
        extrasaction="ignore",
    )

    w.writeheader()
    w.writerows(anchors)


# ------------------------------------------------------------
# 5. Analysis state
# ------------------------------------------------------------

state = {
    "n_anchors": len(anchors),
    "genotype_availability":
        availability.get("status", "UNKNOWN"),
    "phenotype_prototype": "PASS",
    "integrated_evidence": "FROZEN",
    "controlled_validation":
        contract["current_status"],
    "blocker": (
        "Controlled-access KoGES individual genotype data"
        if contract["current_status"]
        == "READY_PENDING_CONTROLLED_DATA"
        else None
    ),
}

(OUT / "KOGES_CONTROLLED_ANALYSIS_STATE.json").write_text(
    json.dumps(state, indent=2)
)


# ------------------------------------------------------------
# 6. Human-readable status
# ------------------------------------------------------------

print("=" * 64)
print("KOGES CONTROLLED VALIDATION PREPARATION")
print("=" * 64)

print("anchors =", len(anchors))
print(
    "genotype =",
    state["genotype_availability"]
)
print(
    "controlled_validation =",
    state["controlled_validation"]
)
print(
    "blocker =",
    state["blocker"]
)

print()

for r in anchors:
    print(
        r["gene_symbol"],
        r["rsid"],
        r["chrom_hg19"],
        r["pos_hg19"],
        r["ref_allele"],
        r["effect_allele_alt"],
        sep="\t",
    )
