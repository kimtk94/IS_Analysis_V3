#!/usr/bin/env python3
"""MUSCLE Stage 1 v1.1: audit resistance-training trainability GWAS loci.

The audit deliberately separates:
1) conventional GWAS significance (P < 5e-8),
2) each paper's analysis threshold,
3) thresholds/claims stated in the abstract,
4) functional evidence used only for downstream prioritization.

Only the Python standard library is required.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Optional

VERSION = "1.1"
CONVENTIONAL_GWAS_THRESHOLD = 5e-8
SUGGESTIVE_THRESHOLD = 1e-5
CADD_SUGGESTED_THRESHOLD = 12.37

REQUIRED_COLUMNS = {
    "study_id", "year", "cohort_type", "population", "phenotype", "n",
    "rsid", "chr", "pos", "ea", "nea", "maf", "beta", "pve_pct",
    "cadd", "regulomedb", "gene", "function", "p", "methods_threshold",
    "abstract_claimed_threshold", "paper_calls_genomewide",
    "skeletal_muscle_eqtl_reported", "source_url", "notes",
}

MISSING = {"", "NA", "N/A", ".", "None", "null"}


def as_float(value: str, field: str) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid numeric value for {field}: {value!r}") from exc
    if not math.isfinite(out):
        raise ValueError(f"Non-finite numeric value for {field}: {value!r}")
    return out


def as_optional_float(value: str, field: str) -> Optional[float]:
    if value is None or str(value).strip() in MISSING:
        return None
    return as_float(str(value).strip(), field)


def as_bool(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def tier_for_p(p: float) -> str:
    if p < CONVENTIONAL_GWAS_THRESHOLD:
        return "A_genomewide"
    if p < SUGGESTIVE_THRESHOLD:
        return "B_suggestive"
    return "C_subthreshold"


def regulomedb_priority(value: str) -> bool:
    value = (value or "").strip().lower()
    if value in MISSING:
        return False
    return value.startswith("1") or value.startswith("2")


def audit_row(row: Dict[str, str]) -> Dict[str, str]:
    p = as_float(row["p"], "p")
    methods_threshold = as_float(row["methods_threshold"], "methods_threshold")
    abstract_threshold = as_optional_float(
        row["abstract_claimed_threshold"], "abstract_claimed_threshold"
    )
    cadd = as_optional_float(row["cadd"], "cadd")
    paper_calls_genomewide = as_bool(row["paper_calls_genomewide"])
    muscle_eqtl = as_bool(row["skeletal_muscle_eqtl_reported"])

    recalculated_genomewide = p < CONVENTIONAL_GWAS_THRESHOLD
    methods_threshold_pass = p < methods_threshold
    abstract_threshold_pass = (
        "" if abstract_threshold is None else ("1" if p < abstract_threshold else "0")
    )
    internal_threshold_discordance = (
        abstract_threshold is not None
        and not math.isclose(methods_threshold, abstract_threshold, rel_tol=0.0, abs_tol=0.0)
    )
    nonstandard_genomewide_threshold = (
        paper_calls_genomewide
        and not math.isclose(
            methods_threshold, CONVENTIONAL_GWAS_THRESHOLD, rel_tol=0.0, abs_tol=0.0
        )
    )

    evidence_points = 0
    evidence_points += 2 if recalculated_genomewide else 1 if p < SUGGESTIVE_THRESHOLD else 0
    evidence_points += 2 if muscle_eqtl else 0
    evidence_points += 1 if (cadd is not None and cadd >= CADD_SUGGESTED_THRESHOLD) else 0
    evidence_points += 1 if regulomedb_priority(row["regulomedb"]) else 0

    if evidence_points >= 4:
        priority = "HIGH_FUNCTIONAL_PRIORITY"
    elif evidence_points >= 2:
        priority = "MEDIUM_FUNCTIONAL_PRIORITY"
    else:
        priority = "EXPLORATORY"

    out = dict(row)
    out.update({
        "recalculated_tier": tier_for_p(p),
        "recalculated_genomewide": "1" if recalculated_genomewide else "0",
        "methods_threshold_pass": "1" if methods_threshold_pass else "0",
        "abstract_threshold_pass": abstract_threshold_pass,
        "paper_internal_threshold_discordance": "1" if internal_threshold_discordance else "0",
        "nonstandard_genomewide_threshold": "1" if nonstandard_genomewide_threshold else "0",
        "cadd_ge_12_37": (
            "" if cadd is None else ("1" if cadd >= CADD_SUGGESTED_THRESHOLD else "0")
        ),
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
        raise ValueError(f"Input contains no loci: {path}")

    keys = [(r["study_id"], r["cohort_type"], r["rsid"]) for r in rows]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate study/cohort/rsid rows detected")

    return rows


def write_tsv(path: Path, rows: Iterable[Dict[str, str]]) -> None:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows: List[Dict[str, str]], dataset_role: str) -> Dict[str, object]:
    tier_counts = Counter(r["recalculated_tier"] for r in rows)
    priority_counts = Counter(r["stage1_priority"] for r in rows)
    study_counts = Counter(r["study_id"] for r in rows)

    internal = [
        r["rsid"] for r in rows if r["paper_internal_threshold_discordance"] == "1"
    ]
    nonstandard = [
        r["rsid"] for r in rows if r["nonstandard_genomewide_threshold"] == "1"
    ]
    abstract_fail = [
        r["rsid"]
        for r in rows
        if r["abstract_threshold_pass"] == "0"
    ]
    medium_or_high = [
        r["rsid"]
        for r in rows
        if r["stage1_priority"] in {"MEDIUM_FUNCTIONAL_PRIORITY", "HIGH_FUNCTIONAL_PRIORITY"}
    ]

    return {
        "version": VERSION,
        "dataset_role": dataset_role,
        "n_loci": len(rows),
        "study_counts": dict(sorted(study_counts.items())),
        "conventional_gwas_threshold": CONVENTIONAL_GWAS_THRESHOLD,
        "suggestive_threshold": SUGGESTIVE_THRESHOLD,
        "tier_counts": dict(sorted(tier_counts.items())),
        "priority_counts": dict(sorted(priority_counts.items())),
        "n_methods_threshold_pass": sum(r["methods_threshold_pass"] == "1" for r in rows),
        "n_paper_internal_threshold_discordant": len(internal),
        "paper_internal_threshold_discordant_rsids": internal,
        "n_fail_abstract_claimed_threshold": len(abstract_fail),
        "fail_abstract_claimed_threshold_rsids": abstract_fail,
        "n_nonstandard_genomewide_threshold": len(nonstandard),
        "medium_or_high_functional_priority_rsids": medium_or_high,
        "interpretation": (
            "Tier A uses the conventional P<5e-8 threshold. Tier B is suggestive "
            "(5e-8<=P<1e-5). Passing a paper-specific threshold does not make a locus "
            "conventionally genome-wide significant. Functional priority is triage only."
        ),
    }


def write_report(
    path: Path,
    primary_summary: Dict[str, object],
    primary_rows: List[Dict[str, str]],
    comparator_summary: Optional[Dict[str, object]] = None,
) -> None:
    ranked = sorted(
        primary_rows,
        key=lambda r: (-int(r["functional_evidence_points"]), float(r["p"]))
    )
    lines = [
        "# MUSCLE Stage 1 v1.1 GWAS Audit",
        "",
        f"- Primary RT loci audited: **{primary_summary['n_loci']}**",
        f"- Study counts: **{primary_summary['study_counts']}**",
        f"- Conventional genome-wide threshold: **P < {CONVENTIONAL_GWAS_THRESHOLD:g}**",
        f"- Suggestive threshold: **P < {SUGGESTIVE_THRESHOLD:g}**",
        f"- Loci passing their paper's Methods threshold: **{primary_summary['n_methods_threshold_pass']}**",
        f"- Abstract↔Methods threshold discordances: **{primary_summary['n_paper_internal_threshold_discordant']}**",
        f"- Loci failing an explicitly stated abstract threshold: **{primary_summary['n_fail_abstract_claimed_threshold']}**",
        "",
        "## Conventional reclassification",
        "",
    ]
    for key, value in sorted(primary_summary["tier_counts"].items()):
        lines.append(f"- {key}: {value}")

    lines += [
        "",
        "## Functional follow-up priority",
        "",
        "| Study | rsID | Gene | P | Tier | Methods pass | Abstract pass | Points | Priority |",
        "|---|---|---|---:|---|---:|---:|---:|---|",
    ]
    for r in ranked:
        lines.append(
            f"| {r['study_id']} | {r['rsid']} | {r['gene']} | {float(r['p']):.3g} | "
            f"{r['recalculated_tier']} | {r['methods_threshold_pass']} | "
            f"{r['abstract_threshold_pass'] or 'NA'} | {r['functional_evidence_points']} | "
            f"{r['stage1_priority']} |"
        )

    if comparator_summary is not None:
        lines += [
            "",
            "## HIIT comparator",
            "",
            f"- Comparator loci: **{comparator_summary['n_loci']}**",
            f"- Conventional tier counts: **{comparator_summary['tier_counts']}**",
            "- HIIT loci are retained only for exercise-mode comparison and are not part of the primary RT discovery set.",
        ]

    lines += [
        "",
        "## Interpretation guardrails",
        "",
        "The 2024 and 2026 studies use P<1e-5 in their analysis methods. "
        "That threshold is treated as exploratory/suggestive rather than the conventional "
        "GWAS threshold of P<5e-8. The 2026 abstract additionally states P<5e-8, creating "
        "an internal threshold inconsistency because its reported lead SNP P values are larger.",
        "",
        "Functional-priority labels only decide which loci advance first to skeletal-muscle "
        "eQTL/sQTL and exercise multi-omics validation. They are not causal claims.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(
    input_path: Path,
    outdir: Path,
    comparator_path: Optional[Path] = None,
) -> Dict[str, object]:
    primary_raw = read_tsv(input_path)
    primary = [audit_row(r) for r in primary_raw]
    outdir.mkdir(parents=True, exist_ok=True)

    write_tsv(outdir / "MUSCLE_STAGE1_GWAS_AUDIT.tsv", primary)
    primary_summary = summarize(primary, "primary_RT")

    comparator_summary = None
    if comparator_path is not None:
        comparator_raw = read_tsv(comparator_path)
        comparator = [audit_row(r) for r in comparator_raw]
        write_tsv(outdir / "MUSCLE_STAGE1_HIIT_COMPARATOR_AUDIT.tsv", comparator)
        comparator_summary = summarize(comparator, "HIIT_comparator")
        (outdir / "MUSCLE_STAGE1_HIIT_COMPARATOR_SUMMARY.json").write_text(
            json.dumps(comparator_summary, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )

    combined_summary = dict(primary_summary)
    if comparator_summary is not None:
        combined_summary["comparator"] = comparator_summary

    (outdir / "MUSCLE_STAGE1_SUMMARY.json").write_text(
        json.dumps(combined_summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    write_report(
        outdir / "MUSCLE_STAGE1_REPORT.md",
        primary_summary,
        primary,
        comparator_summary,
    )
    return combined_summary


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, type=Path, help="Primary RT locus TSV")
    p.add_argument("--outdir", required=True, type=Path, help="Output directory")
    p.add_argument(
        "--comparator",
        type=Path,
        default=None,
        help="Optional HIIT comparator locus TSV",
    )
    return p


def main() -> None:
    args = build_parser().parse_args()
    summary = run(args.input, args.outdir, args.comparator)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
