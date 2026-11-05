#!/usr/bin/env python3
"""Calculate participant-level protein GRS from reviewed KoGES dosages.

PGS_ik = sum_j G_ij * beta_j^pQTL

The dosage input must already use the dosage allele recorded in Stage5C.
Primary score completeness is explicit; incomplete scores are retained only
for audit and must not silently replace the intended full external GRS.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
from collections import defaultdict
from pathlib import Path

ROOT = Path(os.environ.get("IS_ANALYSIS_ROOT", "/srv/is-analysis"))
STAGE5 = ROOT / "results/metabolic_resilience/stage5_koges"
WEIGHTDIR = STAGE5 / "harmonized_grs_weights"
DOSAGE = STAGE5 / "STAGE5I_KOGES_DOSAGE_LONG.tsv.gz"
TEMPLATE = STAGE5 / "STAGE5I_KOGES_DOSAGE_LONG_TEMPLATE.tsv"
OUT = STAGE5 / "STAGE5I_PROTEIN_GRS_LONG.tsv.gz"
SUMMARY = STAGE5 / "STAGE5I_PROTEIN_GRS_SUMMARY.json"

if not DOSAGE.exists():
    with TEMPLATE.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["participant_id","koges_variant_id","dosage"],
            delimiter="\t",
            lineterminator="\n",
        )
        w.writeheader()

    payload = {
        "status": "HOLD_DOSAGE_LONG_MISSING",
        "template": str(TEMPLATE),
        "dosage_definition": "effect of dosage_allele as reviewed in Stage5C",
    }
    SUMMARY.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    raise SystemExit(2)

weights = {}
gene_expected = {}

for path in sorted(WEIGHTDIR.glob("*.koges_harmonized_grs_weights.tsv")):
    gene = path.name.split(".")[0]
    with path.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    usable = []
    for r in rows:
        vid = r.get("koges_variant_id", "")
        beta = r.get("score_weight_per_koges_dosage", "")
        if not vid or beta in ("", "NA", None):
            continue
        try:
            b = float(beta)
        except Exception:
            continue
        usable.append((vid, b))

    gene_expected[gene] = len(usable)
    for vid, beta in usable:
        weights[(gene, vid)] = beta

if not weights:
    payload = {"status": "HOLD_NO_HARMONIZED_GRS_WEIGHTS"}
    SUMMARY.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    raise SystemExit(2)

variant_to_genes = defaultdict(list)
for (gene, vid), beta in weights.items():
    variant_to_genes[vid].append((gene, beta))

score = defaultdict(float)
seen = defaultdict(set)
participants = set()

with gzip.open(DOSAGE, "rt", encoding="utf-8", newline="") as f:
    rd = csv.DictReader(f, delimiter="\t")
    required = {"participant_id","koges_variant_id","dosage"}
    missing = required - set(rd.fieldnames or [])
    if missing:
        raise SystemExit(f"Missing dosage columns: {sorted(missing)}")

    for r in rd:
        pid = r["participant_id"]
        vid = r["koges_variant_id"]
        participants.add(pid)

        if vid not in variant_to_genes:
            continue
        try:
            dosage = float(r["dosage"])
        except Exception:
            continue
        if not math.isfinite(dosage):
            continue

        for gene, beta in variant_to_genes[vid]:
            key = (pid, gene)
            if vid in seen[key]:
                continue
            score[key] += dosage * beta
            seen[key].add(vid)

rows = []
for pid in sorted(participants):
    for gene in sorted(gene_expected):
        key = (pid, gene)
        expected = gene_expected[gene]
        observed = len(seen[key])
        rows.append({
            "participant_id": pid,
            "gene": gene,
            "protein_grs": score[key] if observed else "",
            "n_variants_expected": expected,
            "n_variants_observed": observed,
            "score_complete": int(expected > 0 and observed == expected),
            "score_status": (
                "PASS_COMPLETE"
                if expected > 0 and observed == expected
                else ("REVIEW_PARTIAL" if observed > 0 else "FAIL_ZERO")
            ),
        })

fields = [
    "participant_id","gene","protein_grs",
    "n_variants_expected","n_variants_observed",
    "score_complete","score_status",
]
with gzip.open(OUT, "wt", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)

payload = {
    "status": "PASS" if rows else "FAIL_ZERO",
    "participants": len(participants),
    "genes": len(gene_expected),
    "participant_gene_scores": len(rows),
    "complete_scores": sum(r["score_complete"] == 1 for r in rows),
    "partial_scores": sum(r["score_status"] == "REVIEW_PARTIAL" for r in rows),
    "zero_scores": sum(r["score_status"] == "FAIL_ZERO" for r in rows),
    "formula": "sum(dosage_j * external_pQTL_beta_j)",
    "primary_rule": "use complete intended external GRS unless a sensitivity analysis is explicitly prespecified",
}
SUMMARY.write_text(json.dumps(payload, indent=2) + "\n")
print(json.dumps(payload, indent=2))
