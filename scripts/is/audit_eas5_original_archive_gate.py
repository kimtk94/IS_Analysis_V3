#!/usr/bin/env python3
"""Fail-closed archive-byte provenance gate for historical IS EAS5 GWAS overlap.

No downloads, no SSH, no sample identifiers exported, no accepted coloc or N.
The default GitHub-only run has no source archive: emit scientific HOLD.
If a verified, authorized local archive is mounted, calculate SHA-256 and
independently check the gzip TSV structure before any source-promotion review.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path

REQUIRED = {
    "mapped_genomic_variant_GRCh37", "BBJ_beta_ALT",
    "BBJ_p", "BBJ_allele_exact_Fasta37",
}
VARIANT = re.compile(r"(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y|MT):[1-9][0-9]*:[ACGT]+:[ACGT]+$", re.I)
TARGET_VARIANTS = 25441
MAX_ROWS = 3_000_000

def valid_numeric(value, *, p=False):
    try:
        f = float(value)
    except (ValueError, TypeError, OverflowError):
        return False
    return math.isfinite(f) and (not p or 0 <= f <= 1)

def sha256_file(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for part in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(part)
    return h.hexdigest()

def audit(path=None, *, expected_sha256=None, max_rows=MAX_ROWS):
    result = {
        "schema": "is-eas5-archive-provenance-v1",
        "historical_source_run": "https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39124421779",
        "historical_overlap_target": TARGET_VARIANTS,
        "archive_bytes_available": False,
        "archive_sha256": None,
        "provided_expected_sha256_matched": False,
        "expected_sha256_independently_authenticated": False,
        "original_build_and_ref_alt_independently_verified": False,
        "BBJ_effect_allele_independently_verified": False,
        "EAS_AIS_source_approved": False,
        "BBJ_source_approved": False,
        "study_sample_size_adopted": False,
        "coloc_executed": False,
        "scientific_status": "HOLD",
    }
    if path is None or not path.is_file():
        result["file_gate"] = "HOLD_ARCHIVE_NOT_AVAILABLE_ON_CI_RUNNER"
        return result
    result["archive_bytes_available"] = True
    result["file_bytes"] = path.stat().st_size
    result["archive_sha256"] = sha256_file(path)
    counts = Counter()
    candidate_keys = {}
    try:
        with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if not REQUIRED.issubset(reader.fieldnames or []):
                result["file_gate"] = "HOLD_REQUIRED_COLUMNS_MISSING"
                return result
            for row in reader:
                counts["rows"] += 1
                if counts["rows"] > max_rows:
                    result["file_gate"] = "HOLD_ROW_SCAN_BOUND"
                    return result
                if str(row["BBJ_allele_exact_Fasta37"]).lower() not in ("true", "1", "yes"):
                    continue
                variant = str(row["mapped_genomic_variant_GRCh37"]).strip()
                if not VARIANT.fullmatch(variant):
                    counts["invalid_variant_key"] += 1
                    continue
                chrom, pos, ref, alt = variant.upper().removeprefix("CHR").split(":")
                if ref == alt:
                    counts["invalid_variant_key"] += 1
                    continue
                key = (str(int(chrom)) if chrom.isdigit() else chrom, int(pos), ref, alt)
                if not valid_numeric(row["BBJ_beta_ALT"]) or not valid_numeric(row["BBJ_p"], p=True):
                    counts["missing_or_invalid_source_effect"] += 1
                    continue
                value = (float(row["BBJ_beta_ALT"]), float(row["BBJ_p"]))
                if key not in candidate_keys:
                    candidate_keys[key] = value
                elif candidate_keys[key] != value:
                    counts["conflicting_duplicate_keys"] += 1
    except (OSError, UnicodeError, EOFError, csv.Error, ValueError) as exc:
        result["file_gate"] = "HOLD_ARCHIVE_PARSE_OR_GZIP_ERROR"
        result["failure_class"] = type(exc).__name__
        return result

    counts["distinct_source_eligible_variant_keys"] = len(candidate_keys)
    result["bounded_qc"] = dict(counts)
    if expected_sha256 is None:
        result["file_gate"] = "HOLD_NO_TRUSTED_EXPECTED_SHA256"
    elif not re.fullmatch("[a-fA-F0-9]{64}", expected_sha256):
        result["file_gate"] = "HOLD_INVALID_TRUSTED_SHA256"
    elif result["archive_sha256"].lower() != expected_sha256.lower():
        result["file_gate"] = "HOLD_SOURCE_SHA256_MISMATCH"
    elif counts["conflicting_duplicate_keys"]:
        result["file_gate"] = "HOLD_CONFLICTING_ARCHIVE_DUPLICATES"
    elif len(candidate_keys) != TARGET_VARIANTS:
        result["file_gate"] = "HOLD_SOURCE_ELIGIBLE_VARIANT_COUNT_MISMATCH"
    else:
        result["file_gate"] = "FILE_BYTES_AND_HISTORICAL_ELIGIBLE_KEY_COUNT_MATCH_REVIEW"
        result["provided_expected_sha256_matched"] = True
    return result

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path)
    parser.add_argument("--expected-sha256")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.archive, expected_sha256=args.expected_sha256)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    print("IS_EAS5_ORIGINAL_ARCHIVE_GATE=" + json.dumps(report, sort_keys=True))

if __name__ == "__main__":
    main()
