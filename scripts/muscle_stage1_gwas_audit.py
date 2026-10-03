#!/usr/bin/env python3
"""Stage 1 audit for resistance-training trainability GWAS loci.

The script intentionally separates published claims from recalculated thresholds.
It uses only the Python standard library so it can run in CI and Google Colab.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List

GWAS_THRESHOLD = 5e-8
SUGGESTIVE_THRESHOLD = 1e-5
CADD_SUGGESTED_THRESHOLD = 12.37

REQUIRED_COLUMNS = {
    "study_id", "year", "population", "phenotype", "n", "rsid", "chr", "pos",
    "ea", "nea", "maf", "beta", "pve_pct", "cadd", "regulomedb", "gene",
    "function", "p", "source_threshold", "source_claimed_genomewide",
    "skeletal_muscle_eqtl_reported", "source_url",
}


def as_float(value: str, field: str) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid numeric value for {field}: {value!r}") from exc
    if not math.isfinite(out):
        raise ValueError(f"Non-finite numeric value for {field}: {value!r}")
    return out


def tier_for_p(p: float) -> str:
    if p < GWAS_THRESHOLD:
        return "A_genomewide"
    if p < SUGGESTIVE_THRESHOLD:
        return "B_suggestive"
    return "C_subthreshold"


def regulomedb_priority(value: str) -> bool:
    value = (value or "").strip().lower()
    return value.startswith("1") or value.startswith("2")


def audit_row(row: Dict[str, str]) -> Dict[str, str]:
    p = as_float(row["p"], "p")
    cadd = as_float(row["cadd"], "cadd")
    source_threshold = as_float(row["source_threshold"], "source_threshold")
    claimed = row["source_claimed_genomewide"].strip() in {"1", "true", "True", "YES", "yes"}
    muscle_eqtl = row["skeletal_muscle_eqtl_reported"].strip() in {"1", "true", "True", "YES", "yes"}

    tier = tier_for_p(p)
    recalculated_genomewide = p < GWAS_THRESHOLD
    source_threshold_pass = p < source_threshold
    claim_discordant = claimed and not recalculated_genomewide

    evidence_points = 0
    evidence_points += 2 if recalculated_genomewide else 1 if p < SUGGESTIVE_THRESHOLD else 0
    evidence_points += 2 if muscle_eqtl else 0
    evidence_points += 1 if cadd >= CADD_SUGGESTED_THRESHOLD else 0
    evidence_points += 1 if regulomedb_priority(row["regulomedb"]) else 0

    if evidence_points >= 4:
        priority = "HIGH_FUNCTIONAL_PRIORITY"
    elif evidence_points >= 2:
        priority = "MEDIUM_FUNCTIONAL_PRIORITY"
    else:
        priority = "EXPLORATORY"

    out = dict(row)
    out.update({
        "recalculated_tier": tier,
        "recalculated_genomewide": "1" if recalculated_genomewide else "0",
        "source_threshold_pass": "1" if source_threshold_pass else "0",
        "claim_threshold_discordance": "1" if claim_discordant else "0",
        "cadd_ge_12_37": "1" if cadd >= CADD_SUGGESTED_THRESHOLD else "0",
        "regulomedb_1_or_2": "1" if regulomedb_priority(row["regulomedb"]) else "0",
        "functional_evidence_points": str(evidence_points),
        "stage1_priority": priority,
    })
    return out


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        missing = sorted(REQUIRED_COLUMNS - fields)
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(missing)}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise ValueError("Input contains no loci")
    rsids = [r["rsid"] for r in rows]
    if len(rsids) != len(set(rsids)):
        raise ValueError("Duplicate rsid values detected")
    return rows


def write_tsv(path: Path, rows: Iterable[Dict[str, str]]) -> None:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows: List[Dict[str, str]]) -> Dict[str, object]:
    tiers: Dict[str, int] = {}
    priorities: Dict[str, int] = {}
    for row in rows:
        tiers[row["recalculated_tier"]] = tiers.get(row["recalculated_tier"], 0) + 1
        priorities[row["stage1_priority"]] = priorities.get(row["stage1_priority"], 0) + 1

    discordant = [r["rsid"] for r in rows if r["claim_threshold_discordance"] == "1"]
    high = [r["rsid"] for r in rows if r["stage1_priority"] == "HIGH_FUNCTIONAL_PRIORITY"]
    return {
        "n_loci": len(rows),
        "gwas_threshold": GWAS_THRESHOLD,
        "suggestive_threshold": SUGGESTIVE_THRESHOLD,
        "tier_counts": tiers,
        "priority_counts": priorities,
        "n_claim_threshold_discordant": len(discordant),
        "claim_threshold_discordant_rsids": discordant,
        "high_functional_priority_rsids": high,
        "interpretation": (
            "Tier A requires recalculated p < 5e-8. Tier B is suggestive (5e-8 <= p < 1e-5). "
            "Functional priority is for downstream validation only and is not evidence of causality."
        ),
    }


def write_report(path: Path, summary: Dict[str, object], rows: List[Dict[str, str]]) -> None:
    ranked = sorted(rows, key=lambda r: (-int(r["functional_evidence_points"]), float(r["p"])))
    lines = [
        "# MUSCLE Stage 1 GWAS Audit",
        "",
        f"- Loci audited: **{summary['n_loci']}**",
        f"- Recalculated genome-wide threshold: **P < {GWAS_THRESHOLD:g}**",
        f"- Suggestive threshold: **P < {SUGGESTIVE_THRESHOLD:g}**",
        f"- Published-claim/threshold discordances: **{summary['n_claim_threshold_discordant']}**",
        "",
        "## Recalculated tier counts",
        "",
    ]
    for key, value in sorted(summary["tier_counts"].items()):
        lines.append(f"- {key}: {value}")
    lines += ["", "## Functional follow-up priority", "", "| rsID | Gene | P | Tier | Points | Priority |", "|---|---|---:|---|---:|---|"]
    for r in ranked:
        lines.append(
            f"| {r['rsid']} | {r['gene']} | {float(r['p']):.3g} | {r['recalculated_tier']} | "
            f"{r['functional_evidence_points']} | {r['stage1_priority']} |"
        )
    lines += [
        "",
        "## Interpretation",
        "",
        "A functional-priority label is a triage device for GTEx/MoTrPAC/MetaMEx/GSE277819 follow-up. "
        "It must not be interpreted as causal evidence or as replication of resistance-training response.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(input_path: Path, outdir: Path) -> Dict[str, object]:
    raw = read_tsv(input_path)
    audited = [audit_row(r) for r in raw]
    outdir.mkdir(parents=True, exist_ok=True)

    write_tsv(outdir / "MUSCLE_STAGE1_GWAS_AUDIT.tsv", audited)
    summary = summarize(audited)
    (outdir / "MUSCLE_STAGE1_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    write_report(outdir / "MUSCLE_STAGE1_REPORT.md", summary, audited)
    return summary


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, type=Path, help="Curated/obtained GWAS locus TSV")
    p.add_argument("--outdir", required=True, type=Path, help="Output directory")
    return p


def main() -> None:
    args = build_parser().parse_args()
    summary = run(args.input, args.outdir)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
