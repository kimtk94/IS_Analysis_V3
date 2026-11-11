#!/usr/bin/env python3
"""Validate a stroke GWAS phenotype register without pretending to adjudicate cases.

Output schema consistency PASS is NOT cohort-level phenotype verification.
"""
import argparse
import json
import os
from pathlib import Path

EXPECTED = {
    "GCST90104544": ("any_stroke", 27413, False),
    "GCST90104545": ("ischemic_stroke", 19032, True),
    "GCST90104546": ("cardioembolic_stroke", 926, False),
    "GCST90104547": ("large_artery_stroke", 1735, False),
    "GCST90104548": ("small_vessel_stroke", 5532, False),
}
REQUIRED_UNKNOWN = ("imaging_confirmation", "toast_adjudication", "icd_rules", "incident_or_prevalent")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default="data/metadata/stroke_eas_phenotype_register_2026-11-11.json")
    parser.add_argument("--output-dir", default="artifacts/is_eas_phenotypes")
    args = parser.parse_args()
    payload = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    assert payload["evidence_type"] == "published_aggregate_gwas_metadata_only"
    assert payload["record_source_status"] == "COHORT_CASE_DEFINITION_UNVERIFIED"
    assert payload["research_primary_outcome"] == "ischemic_stroke"
    records = payload["eas_studies"]
    assert len(records) == len(EXPECTED)
    assert {x["accession"] for x in records} == set(EXPECTED)
    assert len([x for x in records if x["primary"]]) == 1
    for r in records:
        phenotype, case_count, primary = EXPECTED[r["accession"]]
        assert r["phenotype"] == phenotype
        assert r["published_cases"] == case_count
        assert r["published_controls"] == 237242
        assert r["primary"] is primary
        assert r["case_ascertainment"] == "NOT_VERIFIED"
        for field in REQUIRED_UNKNOWN:
            assert r[field] is None, f"{r['accession']}: avoid ungrounded {field}"
    lines = [
        "# Stroke EAS case-definition provenance audit",
        "",
        "**SCHEMA_PASS / PHENOTYPE_COHORT_HOLD** – clinical cohort diagnosis details not verified.",
        "",
        "| ID | Published trait | Cases | Cohort diagnostic basis |",
        "|---|---|---:|---|",
    ]
    lines.extend(
        f"| {r['accession']} | {r['phenotype']} | {r['published_cases']:,} | {r['case_ascertainment']} |"
        for r in records
    )
    lines.extend([
        "",
        "Do not interpret a GWAS Catalog trait name as proof of imaging-confirmed TOAST case ascertainment.",
        "Outcome and sample counts are published/registry metadata, not validation of source GWAS bytes.",
    ])
    result = "\n".join(lines) + "\n"
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "phenotype_gate.md").write_text(result, encoding="utf-8")
    (out / "phenotype_gate.json").write_text(
        json.dumps({"schema_validation": "PASS",
                    "cohort_case_definition": "HOLD",
                    "imaging_status": "UNVERIFIED",
                    "primary_outcome": "ischemic_stroke"}, indent=2) + "\n",
        encoding="utf-8"
    )
    print(result)
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a", encoding="utf-8") as stream:
            stream.write(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
