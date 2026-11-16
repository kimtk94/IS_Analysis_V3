#!/usr/bin/env python3
"""Fail-closed evidence gate for *locally supplied, authorized* SBP/T2D GWAS archives.

Do not fetch, extract, commit, or inspect private participant-level records.
Only examine archive envelope, bounded header and cryptographic digest.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import re
import zipfile
from pathlib import Path

STUDIES = {
    "SBP": {"ICBP_Evangelou2018_original_GRASP", "IEU_OpenGWAS_ieu-b-38"},
    "T2D": {"T2DGGI_Suzuki2024_EUR"},
}
PUBLISHED_T2D_MD5 = "46f13bf073c5560d9b068b623229ff01"
MAX_RAW_FILE_BYTES = 1_500_000_000
MAX_HEADER_BYTES = 131072
MAX_MEMBER_LIST = 200


def ledger_row(ledger: Path, trait: str, study: str) -> dict:
    if study not in STUDIES.get(trait, set()):
        raise ValueError("unsupported trait/study selection")
    with ledger.open(encoding="utf-8-sig", newline="") as handle:
        matches = [x for x in csv.DictReader(handle, delimiter="\t")
                   if x.get("trait") == trait and x.get("study") == study]
    if len(matches) != 1 or matches[0].get("approved_as_ldsc_input") != "NO":
        raise ValueError("ambiguous or already-approved ledger source")
    return matches[0]


def digest_file(path: Path, max_bytes: int = MAX_RAW_FILE_BYTES):
    size = path.stat().st_size
    if size < 32 or size > max_bytes:
        raise ValueError("input size outside bounded range")
    sha = hashlib.sha256()
    md5 = hashlib.md5(usedforsecurity=False)
    count = 0
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            count += len(chunk)
            if count > max_bytes:
                raise ValueError("file grew beyond byte limit")
            sha.update(chunk)
            md5.update(chunk)
    if count != size:
        raise ValueError("input changed during hash")
    return {"bytes": count, "sha256": sha.hexdigest(), "md5": md5.hexdigest()}


def examine_gzip_header(path: Path, is_vcf: bool):
    with gzip.open(path, "rb") as handle:
        consumed = 0
        metadata = []
        while consumed < MAX_HEADER_BYTES:
            line = handle.readline(8192)
            consumed += len(line)
            if len(line) >= 8192:
                raise ValueError("oversized gzip header line")
            if not line:
                raise ValueError("no tabular header in source")
            decoded = line.decode("utf-8-sig", "strict").rstrip("\r\n")
            if is_vcf and decoded.startswith("##"):
                metadata.append(decoded[:100])
                continue
            if is_vcf and not decoded.startswith("#CHROM"):
                raise ValueError("no VCF #CHROM header")
            header = decoded.lstrip("#").split("\t")
            if len(header) < 4 or len(header) > 500:
                raise ValueError("unreasonable header width")
            return {"type": "VCF" if is_vcf else "TSV",
                    "header_fields": header[:50], "metadata_lines_sampled": len(metadata)}
        raise ValueError("header exceeds byte cap")


def examine_zip(path: Path):
    with zipfile.ZipFile(path) as z:
        members = z.infolist()
        if not members or len(members) > MAX_MEMBER_LIST:
            raise ValueError("unexpected ZIP member count")
        names = []
        sum_raw = 0
        for info in members:
            if info.flag_bits & 1:
                raise ValueError("encrypted ZIP members cannot be audited")
            if info.file_size > 15_000_000_000:
                raise ValueError("ZIP member exceeds allowed size")
            if info.filename.startswith(("/", "\\")) or ".." in Path(info.filename).parts:
                raise ValueError("unsafe ZIP member name")
            sum_raw += info.file_size
            names.append({"name": info.filename[:220], "uncompressed_bytes": info.file_size,
                          "compressed_bytes": info.compress_size})
        return {"type": "ZIP_LIST_ONLY", "members": names, "sum_uncompressed_bytes": sum_raw,
                "zip_members_extracted": False, "crc_uncompressed_verified": False}


def check(trait: str, study: str, file: Path, ledger: Path, *, repo_root: Path):
    row = ledger_row(ledger, trait, study)
    path = file.expanduser().resolve(strict=True)
    root = repo_root.resolve()
    if path == root or root in path.parents:
        raise ValueError("supplied raw file must be kept outside Git repository")
    if not path.is_file():
        raise ValueError("input is not an ordinary file")
    is_zip = path.suffix.lower() == ".zip"
    is_vcf = path.name.lower().endswith(".vcf.gz")
    if trait == "T2D" and not is_zip:
        raise ValueError("selected T2DGGI release expects user-authorized ZIP")
    if trait == "SBP" and (is_zip or not path.name.lower().endswith(".gz")):
        raise ValueError("SBP file should be original gz or OpenGWAS VCF.gz")
    digest = digest_file(path)
    envelope = examine_zip(path) if is_zip else examine_gzip_header(path, is_vcf=is_vcf)
    observed_md5 = digest["md5"].lower()
    if trait == "T2D":
        match = observed_md5 == PUBLISHED_T2D_MD5
        checksum = {"published_md5": PUBLISHED_T2D_MD5,
                    "actual_md5": observed_md5, "matches_published_md5": match}
        verdict = "PASS_ARCHIVE_MD5_ONLY" if match else "HOLD_PUBLISHED_MD5_MISMATCH"
    else:
        checksum = {"published_md5": None, "actual_md5": observed_md5,
                    "matches_published_md5": None}
        verdict = "HOLD_NO_PUBLISHER_CHECKSUM"
    return {
        "schema": "ckm-sbp-t2d-user-supplied-input-v1",
        "checked_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "trait": trait, "study": study,
        "source_distribution": row.get("verified_landing_url"),
        "raw_local_path_recorded": False,
        "access_authorized_by_user": "USER_SUPPLIED_FILE_ASSUMED",
        "source_status": verdict,
        "archive_integrity": digest,
        "publisher_checksum": checksum,
        "container_schema": envelope,
        "original_allele_schema_approved": False,
        "outcome_overlap_checked": False,
        "ready_for_ldsc": False,
        "scientific_status": "HOLD",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trait", required=True, choices=["SBP", "T2D"])
    ap.add_argument("--study", required=True)
    ap.add_argument("--file", required=True, type=Path)
    ap.add_argument("--ledger", type=Path, default=Path("config/ckm_stroke/p1_exact_source_endpoints.tsv"))
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    try:
        result = check(args.trait, args.study, args.file, args.ledger, repo_root=args.repo_root)
    except Exception as exc:
        result = {"trait": args.trait, "study": args.study,
                  "source_status": "HOLD_INPUT_UNVERIFIED", "error_type": type(exc).__name__,
                  "error_detail": ("input file not found" if isinstance(exc, FileNotFoundError)
                                   else str(exc)[:200]), "ready_for_ldsc": False,
                  "scientific_status": "HOLD", "raw_local_path_recorded": False}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"trait": result["trait"], "study": result["study"],
                      "status": result["source_status"],
                      "ready_for_ldsc": result["ready_for_ldsc"]}))
    if result["source_status"] not in ("PASS_ARCHIVE_MD5_ONLY", "HOLD_NO_PUBLISHER_CHECKSUM"):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
