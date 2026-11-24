#!/usr/bin/env python3
"""Validate KoGES genotype coverage for external cis-pQTL GRS weights.

This stage does not calculate participant scores. It verifies that the exact
external-weight variants can be represented in the reviewed KoGES genotype
variant manifest and harmonizes the effect allele orientation.
"""

from __future__ import annotations

import csv
import json
import os
from collections import defaultdict
from pathlib import Path

ROOT = Path(os.environ.get("IS_ANALYSIS_ROOT", "/srv/is-analysis"))
STAGE5 = ROOT / "results/metabolic_resilience/stage5_koges"
WEIGHTS = STAGE5 / "grs_weights"
MANIFEST = STAGE5 / "STAGE5C_KOGES_GENOTYPE_VARIANT_MANIFEST.tsv"
TEMPLATE = STAGE5 / "STAGE5C_KOGES_GENOTYPE_VARIANT_MANIFEST_TEMPLATE.tsv"
OUTDIR = STAGE5 / "harmonized_grs_weights"
SUMMARY = STAGE5 / "STAGE5C_GRS_VARIANT_COVERAGE.json"
TABLE = STAGE5 / "STAGE5C_GRS_VARIANT_COVERAGE.tsv"

OUTDIR.mkdir(parents=True, exist_ok=True)

MFIELDS = [
    "variant_id", "rsid", "chrom_hg19", "pos_hg19",
    "allele_a", "allele_b", "dosage_allele", "info", "status",
]

if not MANIFEST.exists():
    with TEMPLATE.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=MFIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
    payload = {
        "status": "HOLD_GENOTYPE_VARIANT_MANIFEST_MISSING",
        "template": str(TEMPLATE),
        "instruction": "Populate from reviewed KoGES genotype/imputation metadata; do not infer alleles.",
    }
    SUMMARY.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    raise SystemExit(2)

with MANIFEST.open("r", encoding="utf-8", newline="") as f:
    gvars = list(csv.DictReader(f, delimiter="\t"))

by_rsid = defaultdict(list)
by_pos = defaultdict(list)
for r in gvars:
    if r.get("status") not in ("", "PASS"):
        continue
    if r.get("rsid") not in ("", ".", "NA"):
        by_rsid[r["rsid"]].append(r)
    try:
        by_pos[(str(r["chrom_hg19"]).replace("chr", ""), int(float(r["pos_hg19"])))] .append(r)
    except Exception:
        pass

summary = []

for wf in sorted(WEIGHTS.glob("*.protein_cis_pqtl_grs_weights.tsv")):
    gene = wf.name.split(".")[0]
    with wf.open("r", encoding="utf-8", newline="") as f:
        weights = list(csv.DictReader(f, delimiter="\t"))

    harmonized = []
    for w in weights:
        cand = []
        rsid = w.get("rsid", "")
        if rsid not in ("", ".", "NA"):
            cand = by_rsid.get(rsid, [])
        if not cand:
            try:
                key = (str(w["chrom_hg19"]).replace("chr", ""), int(float(w["pos_hg19"])))
                cand = by_pos.get(key, [])
            except Exception:
                cand = []

        ea = w["effect_allele"].upper()
        oa = w["other_allele"].upper()
        hit = None
        orientation = ""

        for g in cand:
            ga = g.get("allele_a", "").upper()
            gb = g.get("allele_b", "").upper()
            if {ga, gb} != {ea, oa}:
                continue

            dosage = g.get("dosage_allele", "").upper()
            if dosage == ea:
                orientation = "DIRECT_DOSAGE_EFFECT"
            elif dosage == oa:
                orientation = "FLIP_DOSAGE_TO_EFFECT"
            else:
                orientation = "REVIEW_DOSAGE_ALLELE"
            hit = g
            break

        if hit:
            rr = dict(w)
            rr["koges_variant_id"] = hit.get("variant_id", "")
            rr["koges_dosage_allele"] = hit.get("dosage_allele", "")
            rr["orientation"] = orientation
            rr["score_weight_per_koges_dosage"] = (
                w["weight_beta_protein"]
                if orientation == "DIRECT_DOSAGE_EFFECT"
                else (
                    str(-float(w["weight_beta_protein"]))
                    if orientation == "FLIP_DOSAGE_TO_EFFECT"
                    else ""
                )
            )
            harmonized.append(rr)

    dst = OUTDIR / f"{gene}.koges_harmonized_grs_weights.tsv"
    fields = [
        "gene", "rsid", "chrom_hg19", "pos_hg19",
        "effect_allele", "other_allele", "weight_beta_protein", "F",
        "koges_variant_id", "koges_dosage_allele", "orientation",
        "score_weight_per_koges_dosage",
    ]
    with dst.open("w", encoding="utf-8", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        wr.writeheader()
        wr.writerows(harmonized)

    k = len(weights)
    kh = len(harmonized)
    if k == 0 or kh == 0:
        status = "FAIL_ZERO"
    elif kh == k and all(r["orientation"] != "REVIEW_DOSAGE_ALLELE" for r in harmonized):
        status = "PASS_COMPLETE"
    else:
        status = "REVIEW_PARTIAL"

    summary.append({
        "gene": gene,
        "weights_total": k,
        "harmonized": kh,
        "coverage_fraction": (kh / k) if k else 0,
        "status": status,
        "output": str(dst),
    })

fields = ["gene", "weights_total", "harmonized", "coverage_fraction", "status", "output"]
with TABLE.open("w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(summary)

payload = {
    "status": "PASS" if summary and all(r["status"] == "PASS_COMPLETE" for r in summary) else "REVIEW",
    "genes": len(summary),
    "complete": sum(r["status"] == "PASS_COMPLETE" for r in summary),
    "partial": sum(r["status"] == "REVIEW_PARTIAL" for r in summary),
    "zero": sum(r["status"] == "FAIL_ZERO" for r in summary),
    "note": "Partial coverage changes the intended external GRS and requires explicit review.",
}
SUMMARY.write_text(json.dumps(payload, indent=2) + "\n")
print(json.dumps(payload, indent=2))
