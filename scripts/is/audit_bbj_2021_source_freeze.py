#!/usr/bin/env python3
"""Validate a *historical Actions log transcription*. Never claim new source recomputation.

No source downloads, no private sample data, no GWAS N adoption and no coloc.
"""
import argparse
import json
from pathlib import Path

SOURCE_RUN = "https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39124421779"
SOURCE_ZIP_SHA256 = "a8032286cf2206a02c746b138a263930dc7d371c57b94352820005f8a32c8d2b"
REQUIRED_KEYS = (
    "archive_candidate_variant_ids", "official_source_rows",
    "same_coordinate_allele_pair", "comparable_effect_p",
    "same_p_with_log_tolerance", "same_beta_sign_and_magnitude",
    "opposite_beta_sign_same_magnitude", "p_and_beta_match_either_sign",
    "hypothesis_effect_A2_agree", "hypothesis_effect_A1_agree",
)


def audit_evidence(evidence):
    if evidence.get("record_type") != "HISTORICAL_ACTIONS_LOG_RECONSTRUCTION_NOT_NEW_GWAS_RUN":
        raise ValueError("Must be explicitly designated historical evidence")
    if evidence.get("original_run") != SOURCE_RUN or evidence.get("original_job_id") != 117666449564:
        raise ValueError("Unrecognized originating Actions job")
    if evidence.get("source_zip_sha256") != SOURCE_ZIP_SHA256 or evidence.get("source_zip_bytes") != 1611973884:
        raise ValueError("Source ZIP historical integrity pin differs")
    if evidence.get("archive_sha256") != "NOT_EXPORTED_BY_COMPARISON_V1":
        raise ValueError("Archive provenance must not be invented")

    raw = evidence.get("original_log_payload", {})
    if raw.get("schema") != "is-bbj-2021-source-effect-comparison-v1":
        raise ValueError("Unknown comparison schema")
    if raw.get("study_release") != "Sakaue2021_BBJ_IS_candidate":
        raise ValueError("Unexpected provisional study release")
    for key in ("genomic_build_verified_independently", "allele_effect_orientation_verified_independently",
                "gwas_source_accepted", "sample_size_adopted", "coloc_executed"):
        if raw.get(key) is not False:
            raise ValueError("Historical HOLD/FAIL-CLOSED flag changed: " + key)
    counts = raw.get("summary", {})
    if not all(isinstance(counts.get(k), int) and not isinstance(counts[k], bool) and counts[k] >= 0 for k in REQUIRED_KEYS):
        raise ValueError("Missing or non-integer count")
    n = counts["comparable_effect_p"]
    if not (
        n == 25441
        == counts["archive_candidate_variant_ids"]
        == counts["same_coordinate_allele_pair"]
        == counts["same_p_with_log_tolerance"]
        == counts["same_beta_sign_and_magnitude"]
        == counts["p_and_beta_match_either_sign"]
        == counts["hypothesis_effect_A2_agree"]
    ):
        raise ValueError("Historical match counts diverged; do not promote")
    if counts["official_source_rows"] != 13435541:
        raise ValueError("Historical total source rows differ")
    if not (
        counts["hypothesis_effect_A1_agree"] == counts["opposite_beta_sign_same_magnitude"] == 20
        and raw.get("p_beta_concordance_fraction_either_sign") == 1.0
    ):
        raise ValueError("Opposite-sign ambiguity or reported fraction changed")
    if any(counts.get(k, 0) for k in (
        "archive_ambiguous_duplicate_rows", "ambiguous_archived_beta_p",
        "duplicate_source_overlap_rows", "different_p", "invalid_source_beta_p"
    )):
        raise ValueError("Unresolved duplicate or invalid source comparison")
    # With A2 matching every observed variant, the 20 A1 matches overlap A2:
    # small absolute-tolerance beta values cannot determine signed orientation.
    dual_compatible = counts["hypothesis_effect_A1_agree"]
    return {
        "schema": "is-bbj2021-historical-evidence-gate-v1",
        "source": SOURCE_RUN,
        "source_zip_sha256": SOURCE_ZIP_SHA256,
        "scientific_status": "HOLD",
        "historical_numeric_concordance_screen": "PASS",
        "new_original_data_comparison_performed": False,
        "source_archive_sha256_recorded_in_v1": False,
        "source_variants_scanned": counts["official_source_rows"],
        "compared_variant_keys": n,
        "p_and_beta_match_either_sign": counts["p_and_beta_match_either_sign"],
        "effect_A2_hypothesis_compatible": counts["hypothesis_effect_A2_agree"],
        "effect_A2_only_under_tolerances": n - dual_compatible,
        "effect_A1_and_A2_both_compatible_under_tolerances": dual_compatible,
        "effect_allele_official_definition_verified": False,
        "genomic_build_independently_verified": False,
        "approved_gwas_source": False,
        "sample_size_adopted": False,
        "coloc_executed": False,
        "limitations": [
            "Aggregate re-audit of existing workflow logs; not a new run over primary files.",
            "20 near-zero beta comparisons are sign-ambiguous under original abs_tol=0.00005.",
            "100 percent concordance applies to a 25,441-variant EAS5 overlap subset, not all GWAS SNPs.",
            "Original v1 result lacks EAS5 archive SHA-256 and independently certified build/effect-allele metadata.",
            "BBJ may overlap GIGASTROKE EAS discovery; not an independent replication cohort."
        ]
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    result = audit_evidence(json.loads(args.input.read_text(encoding="utf-8")))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("IS_BBJ_HISTORICAL_EVIDENCE_GATE=" + json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
