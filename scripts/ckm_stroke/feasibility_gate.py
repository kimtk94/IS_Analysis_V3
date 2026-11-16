#!/usr/bin/env python3
"""Fail-closed CKM→IS research feasibility inventory (metadata only)."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

EXPOSURES = ("BMI", "LDL_C", "HDL_C", "TG", "SBP", "T2D", "EGFR", "UACR")
OUTCOMES = ("LAA", "CE", "SVD")
# Independent 2026 supplementary table, not treated as the authoritative source.
EXTERNAL_SUBTYPE_MAPPING = {
    "LAA": "GCST90104541",
    "CE": "GCST90104542",
    "SVD": "GCST90104543",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"empty TSV header: {path}")
        return list(reader)


def assess(
    phenotypes: list[dict[str, str]],
    sources: list[dict[str, str]],
    tasks: dict,
) -> dict:
    for keys, rows in (
        (("phenotype", "role", "ancestry", "accession", "build", "data_file"), phenotypes),
        (("trait", "ancestry", "source", "kind", "url"), sources),
    ):
        for row in rows:
            if not set(keys).issubset(row):
                raise ValueError(f"required columns missing: {keys}")
    phenotype_lookup: dict[tuple[str, str], dict[str, str]] = {}
    for row in phenotypes:
        key = (row["phenotype"], row["ancestry"])
        if key in phenotype_lookup:
            raise ValueError(f"duplicate phenotype/ancestry: {key}")
        phenotype_lookup[key] = row
    source_lookup: dict[tuple[str, str], dict[str, str]] = {}
    for row in sources:
        key = (row["trait"], row["ancestry"])
        if key in source_lookup:
            raise ValueError(f"duplicate source/ancestry: {key}")
        source_lookup[key] = row

    missing_eur_sources = [
        name for name in EXPOSURES
        if not source_lookup.get((name, "EUR"), {}).get("url", "").strip()
    ]
    exposure_candidates = [
        name for name in EXPOSURES
        if source_lookup.get((name, "EUR"), {}).get("url", "").strip()
    ]
    missing_eur_normalized = [
        name for name in EXPOSURES
        if not phenotype_lookup.get((name, "EUR"), {}).get("data_file", "").strip()
        or phenotype_lookup.get((name, "EUR"), {}).get("data_file", "").startswith("REVIEW_")
    ]
    source_mismatches = [
        f"{name}:{ancestry}" for (name, ancestry), row in source_lookup.items()
        if row.get("url", "").strip() and (name, ancestry) not in phenotype_lookup
    ]
    disputed = []
    for name in OUTCOMES:
        actual = phenotype_lookup.get((name, "EUR"), {}).get("accession", "")
        other = EXTERNAL_SUBTYPE_MAPPING[name]
        if actual and actual != other:
            disputed.append({
                "phenotype": name,
                "manifest_accession": actual,
                "independent_supplement_accession": other,
                "decision": "HOLD_ORIGINAL_CONSORTIUM_METADATA_REQUIRED",
            })
    missing_eas = [
        name for name in OUTCOMES
        if not source_lookup.get((name, "EAS"), {}).get("url", "").strip()
        or phenotype_lookup.get((name, "EAS"), {}).get("accession", "").startswith("REVIEW_")
    ]
    unconfigured_tasks = [
        key for key in ("p2", "p3", "p4", "p5")
        if not isinstance(tasks.get(key), dict)
        or not all(tasks[key].get(field) for field in ("argv", "inputs", "outputs"))
    ]
    missing_outcome_build = [
        name for name in OUTCOMES
        if phenotype_lookup.get((name, "EUR"), {}).get("build", "").startswith("REVIEW_")
    ]
    problems = []
    if missing_eur_sources:
        problems.append("MISSING_EUR_EXPOSURE_SOURCE")
    if missing_eur_normalized:
        problems.append("NO_APPROVED_NORMALIZED_INPUT_PATHS")
    if disputed:
        problems.append("OUTCOME_ACCESSION_LABEL_CONFLICT")
    if missing_outcome_build:
        problems.append("OUTCOME_BUILD_UNVERIFIED")
    if unconfigured_tasks:
        problems.append("DOWNSTREAM_COMMANDS_UNCONFIGURED")
    if missing_eas:
        problems.append("EAS_SUBTYPE_SOURCES_UNSELECTED")
    if source_mismatches:
        problems.append("ORPHAN_SOURCE_ROWS")
    return {
        "technical_status": "PASS_METADATA_AUDITED",
        "scientific_status": "HOLD" if problems else "REVIEW_REQUIRED",
        "can_run_ldsc": False,  # An indexed source URL alone cannot establish P1 completeness.
        "evidence_scope": "repository manifests only, not full-GWAS or cohort QC",
        "eur_exposure_source_urls_present": len(exposure_candidates),
        "eur_exposure_source_urls_total": len(EXPOSURES),
        "eur_source_url_present_traits": exposure_candidates,
        "eur_missing_source_traits": missing_eur_sources,
        "eur_missing_approved_normalized_input_paths": missing_eur_normalized,
        "stroke_eur_accession_conflicts": disputed,
        "stroke_eur_unverified_builds": missing_outcome_build,
        "stroke_eas_missing_source_traits": missing_eas,
        "downstream_unconfigured": unconfigured_tasks,
        "orphan_source_rows": source_mismatches,
        "blocking_reasons": problems,
        "external_comparison": {
            "publication": "2026 Journal of Stroke and Cerebrovascular Diseases supplement Table S2",
            "source_url": "https://media.springernature.com/original/springer-static/esm/art%3A10.1007%2Fs12975-026-01509-7/MediaObjects/12975_2026_1509_MOESM2_ESM.pdf",
            "authority": "independent_secondary_source_not_authoritative",
            "action": "Verify GCST90104541/42 against original GIGASTROKE study-level metadata and file headers before analysis",
        },
        "novelty_prior_art": [{
            "doi": "10.1177/0271678X261462718",
            "published_online": "2026-06-15",
            "overlap": "CKM metabolic signature, incident stroke, MR and colocalization",
            "differentiation_needed": "ischemic stroke subtype heterogeneity and EUR/EAS evidence",
        }],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phenotypes", default="config/ckm_stroke/phenotypes.tsv")
    parser.add_argument("--sources", default="config/ckm_stroke/public_sources.tsv")
    parser.add_argument("--tasks", default="config/ckm_stroke/tasks.json")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    result = assess(read_tsv(Path(args.phenotypes)), read_tsv(Path(args.sources)),
                    json.loads(Path(args.tasks).read_text(encoding="utf-8")))
    target = Path(args.out)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
