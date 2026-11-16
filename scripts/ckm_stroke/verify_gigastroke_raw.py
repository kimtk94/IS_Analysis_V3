#!/usr/bin/env python3
"""Bounded original GIGASTROKE EUR GWAS byte-integrity audit, no data retention.

This audit verifies original EBI metadata + published md5sum sidecar, downloads
one declared GRCh37 gzip stream directly into hashlib, and persists JSON only.
Passing MD5 is NOT allele harmonization, LDSC, MR, or scientific validation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlparse

BASE = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/"
MAX_BYTES = 400_000_000
FIELDS = ("phenotype", "ancestry", "accession", "reported_trait", "build",
          "sample_total", "published_raw_md5", "metadata_url", "verification")


def get_text(url: str, fetcher=urllib.request.urlopen) -> str:
    with fetcher(urllib.request.Request(url, headers={"User-Agent": "CKM-IS-integrity-audit/1.0"}),
                 timeout=30) as response:
        payload = response.read(64_000)
        if len(payload) >= 64_000:
            raise ValueError("metadata size limit exceeded")
        return payload.decode("utf-8-sig")


def load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not set(FIELDS).issubset(set(reader.fieldnames or ())):
            raise ValueError("manifest missing required columns")
        rows = list(reader)
    if len(rows) != 4 or set(row["phenotype"] for row in rows) != {"AIS", "LAA", "CE", "SVD"}:
        raise ValueError("expect exactly AIS/LAA/CE/SVD manifest entries")
    if len({row["accession"] for row in rows}) != 4:
        raise ValueError("duplicated accession")
    return rows


def source_urls(row: dict[str, str]) -> tuple[str, str, str]:
    acc = row["accession"]
    if re.fullmatch(r"GCST[0-9]{8}", acc) is None or row["ancestry"] != "EUR":
        raise ValueError("invalid accession or ancestry")
    folder = BASE + acc + "/"
    raw = folder + acc + "_buildGRCh37.tsv.gz"
    meta = raw + "-meta.yaml"
    sums = folder + "md5sum.txt"
    if row["metadata_url"] != meta or row["build"] != "GRCh37":
        raise ValueError("source URL or build is not official pinned publisher source")
    digest = row["published_raw_md5"].lower()
    if re.fullmatch(r"[0-9a-f]{32}", digest) is None:
        raise ValueError("invalid expected digest")
    return raw, meta, sums


def check_sidecars(row: dict[str, str], metadata: str, md5file: str) -> None:
    acc = row["accession"]
    required = {
        "gwas_id": acc,
        "genome_assembly": "GRCh37",
        "data_file_name": acc + "_buildGRCh37.tsv.gz",
        "data_file_md5sum": row["published_raw_md5"],
    }
    for key, expected in required.items():
        match = re.search(r"^" + re.escape(key) + r":\s*([^\s#]+)\s*$", metadata, re.M)
        if not match or match.group(1).strip('"') != expected:
            raise ValueError("official YAML mismatch for " + key)
    m = re.search(r"^trait_description:\s*\n\s*-\s*(.+)$", metadata, re.M)
    if not m or m.group(1).strip() != row["reported_trait"]:
        raise ValueError("official trait description differs")
    m = re.search(r"^\s*sample_size:\s*(\d+)\s*$", metadata, re.M)
    if not m or m.group(1) != row["sample_total"]:
        raise ValueError("publisher sample count differs")
    filename = acc + "_buildGRCh37.tsv.gz"
    found = [m.group(1).lower() for m in re.finditer(
        r"^([a-fA-F0-9]{32})\s+\*?" + re.escape(filename) + r"\s*$", md5file, re.M)]
    if len(found) != 1 or found[0] != row["published_raw_md5"].lower():
        raise ValueError("publisher md5sum.txt mismatch")


def digest_stream(response, *, max_bytes: int) -> dict:
    md5 = hashlib.md5(usedforsecurity=False)
    read_bytes = 0
    first_block = b""
    length = response.headers.get("Content-Length")
    if length is not None and (not length.isdecimal() or int(length) > max_bytes):
        raise ValueError("official download exceeds byte cap")
    while True:
        block = response.read(1024 * 1024)
        if not block:
            break
        read_bytes += len(block)
        if read_bytes > max_bytes:
            raise ValueError("download exceeded cap")
        md5.update(block)
        if len(first_block) < 1024:
            first_block += block[:1024 - len(first_block)]
    if length is not None and int(length) != read_bytes:
        raise ValueError("incomplete Content-Length")
    if first_block[:2] != b"\x1f\x8b":
        raise ValueError("not gzip")
    if not read_bytes:
        raise ValueError("empty source")
    return {"raw_bytes": read_bytes, "actual_raw_md5": md5.hexdigest(),
            "gzip_magic_verified": True, "http_content_length": int(length) if length else None}


def audit_one(row: dict[str, str], *, max_bytes: int = MAX_BYTES,
              fetcher=urllib.request.urlopen) -> dict:
    raw, meta, sums = source_urls(row)
    check_sidecars(row, get_text(meta, fetcher), get_text(sums, fetcher))
    with fetcher(urllib.request.Request(raw, headers={"User-Agent": "CKM-IS-integrity-audit/1.0"}),
                 timeout=30) as source:
        if getattr(source, "status", 200) != 200:
            raise ValueError("non-200 raw source")
        result = digest_stream(source, max_bytes=max_bytes)
    equal = result["actual_raw_md5"] == row["published_raw_md5"].lower()
    return {
        "phenotype": row["phenotype"], "accession": row["accession"],
        "reported_trait": row["reported_trait"], "build": row["build"],
        "official_metadata_verified": True, "official_md5_sidecar_verified": True,
        "published_md5": row["published_raw_md5"].lower(), "raw_source_url": raw,
        **result,
        "raw_md5_verified": equal,
        "status": "PASS_RAW_SOURCE_INTEGRITY" if equal else "FAIL_RAW_MD5_MISMATCH",
        "scientific_status": "HOLD",
        "limitations": "Byte MD5 does not verify alleles, variant schema, harmonization, LD or MR.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path("config/ckm_stroke/gigastroke_eur_official.tsv"))
    parser.add_argument("--phenotype", choices=("AIS", "LAA", "CE", "SVD"), required=True)
    parser.add_argument("--max-bytes", type=int, default=MAX_BYTES)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    records = load_manifest(args.manifest)
    row = next(r for r in records if r["phenotype"] == args.phenotype)
    report = {"audit_version": "ckm-gigastroke-source-integrity-v1",
              "checked_at_utc": datetime.now(timezone.utc).isoformat(),
              "scope": "publisher metadata + md5 sidecar + exact raw compressed bytes",
              "result": None}
    try:
        if args.max_bytes <= 0 or args.max_bytes > MAX_BYTES:
            raise ValueError("invalid byte cap")
        report["result"] = audit_one(row, max_bytes=args.max_bytes)
    except Exception as exc:
        report["result"] = {"phenotype": args.phenotype, "accession": row["accession"],
                            "status": "HOLD_SOURCE_UNAVAILABLE_OR_INVALID",
                            "scientific_status": "HOLD",
                            "error_type": type(exc).__name__, "error_detail": str(exc)[:240]}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report))
    if report["result"]["status"] != "PASS_RAW_SOURCE_INTEGRITY":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
