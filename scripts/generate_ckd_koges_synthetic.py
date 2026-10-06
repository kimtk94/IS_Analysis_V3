#!/usr/bin/env python3

from __future__ import annotations

import csv
import math
import random
from pathlib import Path


SEED = 20260927
random.seed(SEED)

N = 1200

OUT = Path(
    "/srv/is-analysis/results/ckd/stage4_koges/"
    "controlled_validation/synthetic_test"
)

OUT.mkdir(parents=True, exist_ok=True)

ANCHOR = Path(
    "/srv/is-analysis/results/ckd/stage4_koges/"
    "controlled_validation/KOGES_CONTROLLED_ANCHOR_SPEC.tsv"
)


with ANCHOR.open() as f:
    anchors = list(
        csv.DictReader(
            f,
            delimiter="\t",
        )
    )

assert len(anchors) == 9


# ------------------------------------------------------------
# Synthetic truths
#
# Only selected genes have non-zero effects.
# These values are simulation parameters, NOT CKD findings.
# ------------------------------------------------------------

truth = {
    "ACP1": {
        "slope": +0.45,
        "event_log_hr": -0.25,
    },

    "SDCCAG8": {
        "slope": +0.30,
        "event_log_hr": -0.15,
    },

    "UMOD": {
        "slope": -0.65,
        "event_log_hr": +0.35,
    },
}

for a in anchors:
    truth.setdefault(
        a["gene_symbol"],
        {
            "slope": 0.0,
            "event_log_hr": 0.0,
        },
    )


# ------------------------------------------------------------
# Participants
# ------------------------------------------------------------

subjects = []

for i in range(1, N + 1):

    pid = f"SYN{i:05d}"

    sex = random.choice(["M", "F"])

    age0 = random.uniform(40, 69)

    baseline = random.gauss(
        92 - 0.35 * (age0 - 50),
        10,
    )

    baseline = max(
        45,
        min(130, baseline),
    )

    site = random.choice(
        ["ANSAN", "ANSUNG"]
    )

    batch = random.choice(
        ["B1", "B2", "B3"]
    )

    pcs = {
        f"PC{k}": random.gauss(0, 1)
        for k in range(1, 11)
    }

    subjects.append({
        "participant_id": pid,
        "sex": sex,
        "age0": age0,
        "baseline_egfr": baseline,
        "site": site,
        "batch": batch,
        **pcs,
    })


# ------------------------------------------------------------
# Genotype dosage
# ------------------------------------------------------------

geno_rows = []

geno_by_subject = {}

for a in anchors:

    gene = a["gene_symbol"]

    # synthetic MAF
    maf = random.uniform(
        0.12,
        0.42,
    )

    for s in subjects:

        # Hardy-Weinberg binomial genotype
        dosage = (
            int(random.random() < maf)
            + int(random.random() < maf)
        )

        geno_by_subject[
            (s["participant_id"], gene)
        ] = dosage

        geno_rows.append({
            "participant_id":
                s["participant_id"],

            "gene_symbol":
                gene,

            "target_rsid":
                a["rsid"],

            "analysis_rsid":
                a["rsid"],

            "is_proxy":
                0,

            "proxy_r2":
                "",

            "effect_allele":
                a["effect_allele_alt"],

            "other_allele":
                a["ref_allele"],

            "dosage":
                dosage,

            "info_score":
                0.99,

            "maf":
                maf,
        })


# ------------------------------------------------------------
# Longitudinal phenotype
#
# 5 visits:
# 0, 2, 4, 6, 8 years
# ------------------------------------------------------------

times = [0, 2, 4, 6, 8]

pheno_rows = []


for s in subjects:

    pid = s["participant_id"]

    random_intercept = random.gauss(
        0,
        4,
    )

    random_slope = random.gauss(
        0,
        0.35,
    )

    base_slope = (
        -1.25
        - 0.012 * (s["age0"] - 50)
        + random_slope
    )

    # Sum simulated genetic slope contributions.
    genetic_slope = 0.0

    event_lp = -4.3

    for gene, pars in truth.items():

        g = geno_by_subject[
            (pid, gene)
        ]

        genetic_slope += (
            pars["slope"] * g
        )

        event_lp += (
            pars["event_log_hr"] * g
        )

    # Subject-specific synthetic CKD risk.
    event_lp += (
        0.035 * (s["age0"] - 50)
        - 0.025 * (
            s["baseline_egfr"] - 90
        )
    )

    hazard = math.exp(event_lp)

    event_occurred = False

    for j, t in enumerate(times):

        age = s["age0"] + t

        egfr = (
            s["baseline_egfr"]
            + random_intercept
            + (
                base_slope
                + genetic_slope
            ) * t
            + random.gauss(0, 3.0)
        )

        egfr = max(
            10,
            min(150, egfr),
        )

        # approximate creatinine solely for interface testing
        creatinine = max(
            0.4,
            90 / egfr,
        )

        incident = 0

        if (
            j > 0
            and not event_occurred
        ):

            interval_prob = (
                1
                - math.exp(
                    -hazard * 2
                )
            )

            if random.random() < interval_prob:
                incident = 1
                event_occurred = True

        row = {
            "participant_id": pid,
            "visit_id": f"V{j}",
            "visit_date":
                f"{2010 + int(t):04d}-01-01",
            "time_years": t,
            "age": age,
            "sex": s["sex"],
            "creatinine_mg_dl":
                creatinine,
            "egfr": egfr,
            "baseline_egfr":
                s["baseline_egfr"],
            "baseline_ckd":
                int(
                    s["baseline_egfr"] < 60
                ),
            "incident_ckd":
                incident,
            "interval_index": j,
            "site": s["site"],
            "batch": s["batch"],
        }

        for k in range(1, 11):
            row[f"PC{k}"] = s[f"PC{k}"]

        pheno_rows.append(row)


# ------------------------------------------------------------
# Write
# ------------------------------------------------------------

geno_path = OUT / "synthetic_genotype.tsv"

with geno_path.open(
    "w",
    newline="",
) as f:

    fields = list(
        geno_rows[0].keys()
    )

    w = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
    )

    w.writeheader()
    w.writerows(geno_rows)


pheno_path = OUT / "synthetic_phenotype.tsv"

with pheno_path.open(
    "w",
    newline="",
) as f:

    fields = list(
        pheno_rows[0].keys()
    )

    w = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
    )

    w.writeheader()
    w.writerows(pheno_rows)


truth_path = OUT / "SYNTHETIC_TRUTH.tsv"

with truth_path.open(
    "w",
    newline="",
) as f:

    w = csv.writer(
        f,
        delimiter="\t",
    )

    w.writerow([
        "gene_symbol",
        "true_beta_genotype_time",
        "true_event_log_hr",
    ])

    for gene in sorted(truth):

        w.writerow([
            gene,
            truth[gene]["slope"],
            truth[gene]["event_log_hr"],
        ])


print("seed =", SEED)
print("participants =", N)
print("phenotype_rows =", len(pheno_rows))
print("genotype_rows =", len(geno_rows))
print("genes =", len(anchors))
print("PHENO =", pheno_path)
print("GENO =", geno_path)
print("TRUTH =", truth_path)
