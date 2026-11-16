#!/usr/bin/env python3
"""Compare published GIGASTROKE YAML MD5 to exact original compressed file.

Streams one original file per action, without storing variant rows or requiring
credentials. MD5 is publisher-provided checksum, SHA256 is additional local
fingerprint. Not an LDSC or MAGMA analysis.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import urllib.error
import urllib.request
from urllib.parse import urlparse

FTP_ROOT = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000"
INFO = {
    "EUR": {"accession": "GCST90104540", "cases": 62100, "controls": 1234808},
    "EAS": {"accession": "GCST90104545", "cases": 19032, "controls": 237242},
}
MAX_BYTES = 700 * 1024 * 1024


def source_paths(accession):
    basename = accession + "_buildGRCh37.tsv.gz"
    root = f"{FTP_ROOT}/{accession}/"
    return root + basename, root + basename + "-meta.yaml"


def metadata_values(text):
    selected = {}
    for key in ("gwas_id", "data_file_name", "data_file_md5sum",
                "genome_assembly", "sample_size"):
        match = re.search(r"^\s*" + re.escape(key) + r":\s*(.*?)\s*$", text,
                          flags=re.I | re.M)
        if match:
            selected[key] = match.group(1).strip().strip("'\"")
    return selected


def urlopen_source(url):
    r = urllib.request.Request(url, headers={
        "User-Agent": "IS_Analysis_V3_SourceIndependentChecksumVerifier/1.0"
    })
    response = urllib.request.urlopen(r, timeout=80)
    parsed = urlparse(response.geturl())
    if parsed.scheme != "https" or parsed.netloc != "ftp.ebi.ac.uk":
        response.close()
        raise ValueError("Source redirected to different host")
    return response


def read_bounded_metadata(url):
    with urlopen_source(url) as stream:
        content = stream.read(128001)
    if len(content) > 128000:
        raise ValueError("Metadata too large")
    return metadata_values(content.decode("utf-8", errors="replace"))


def stream_hash_file(url):
    size, md5, sha = 0, hashlib.md5(), hashlib.sha256()
    with urlopen_source(url) as stream:
        status = stream.status
        declared = stream.headers.get("Content-Length")
        if declared and declared.isdigit() and int(declared) > MAX_BYTES:
            raise ValueError("Source exceeds max bytes")
        while chunk := stream.read(4 * 1024 * 1024):
            size += len(chunk)
            if size > MAX_BYTES:
                raise ValueError("Source exceeds max bytes")
            md5.update(chunk)
            sha.update(chunk)
    return {
        "http_status": status,
        "bytes_downloaded": size,
        "declared_bytes": int(declared) if declared and declared.isdigit() else None,
        "local_md5": md5.hexdigest(),
        "local_sha256": sha.hexdigest(),
        "transfer_size_consistent": (
            declared is None or not declared.isdigit() or size == int(declared)
        )
    }


def audit(ancestry):
    study = INFO[ancestry]
    url, metadata_url = source_paths(study["accession"])
    outcome = {
        "ancestry": ancestry,
        "accession": study["accession"],
        "source_url": url,
        "official_metadata_url": metadata_url,
        "independent_publisher_checksum_status": "HOLD",
        "full_row_scan_this_run": False,
        "ancestry_ld_matched": False,
        "magma_analysis": "NOT_RUN",
        "ldsc_seg_analysis": "NOT_RUN",
    }
    try:
        meta = read_bounded_metadata(metadata_url)
        outcome["publisher_metadata"] = meta
        expected = meta.get("data_file_md5sum", "")
        if (meta.get("gwas_id") != study["accession"] or
                meta.get("genome_assembly") != "GRCh37" or
                meta.get("data_file_name") != study["accession"] + "_buildGRCh37.tsv.gz" or
                not re.fullmatch("[0-9a-fA-F]{32}", expected)):
            outcome["independent_publisher_checksum_status"] = "METADATA_REVIEW"
            return outcome
        outcome.update(stream_hash_file(url))
        if outcome["http_status"] != 200 or not outcome["transfer_size_consistent"]:
            outcome["independent_publisher_checksum_status"] = "TRANSFER_REVIEW"
        elif outcome["local_md5"].lower() == expected.lower():
            outcome["independent_publisher_checksum_status"] = "PUBLISHER_MD5_MATCH"
        else:
            outcome["independent_publisher_checksum_status"] = "PUBLISHER_MD5_MISMATCH"
    except (urllib.error.URLError, OSError, TimeoutError, ValueError) as exc:
        outcome["independent_publisher_checksum_status"] = "SOURCE_UNAVAILABLE"
        outcome["error"] = f"{type(exc).__name__}: {str(exc)[:180]}"
    # Note: the literature-derived total sample size does not mean each SNP had
    # that N, especially in meta-GWAS. Keep the illustrative N separate.
    cases, controls = study["cases"], study["controls"]
    outcome["published_cases_not_per_variant"] = cases
    outcome["published_controls_not_per_variant"] = controls
    outcome["published_total_n_not_per_variant"] = cases + controls
    outcome["illustrative_balanced_effective_n_not_per_variant"] = round(
        4 * cases * controls / (cases + controls)
    )
    outcome["sample_size_model_ready"] = "HOLD"
    return outcome


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ancestry", required=True, choices=INFO)
    parser.add_argument("--out-dir", default="artifacts/is_gwas_publisher_checksum")
    args = parser.parse_args()
    result = audit(args.ancestry)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    report = {
        "time_utc": datetime.now(timezone.utc).isoformat(),
        "audit_kind": "PUBLISHER_YAML_MD5_VS_ORIGINAL_REMOTE_BYTES",
        "result": result,
        "biology_analysis": "NOT_RUN",
        "scope_note": "Checksums establish artifact integrity; NOT GWAS phenotype or LD/eQTL readiness.",
    }
    (out / "checksum.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    content = (
        f"## GIGASTROKE publisher MD5 comparison: {args.ancestry}\n\n"
        f"- Accession: {result['accession']}\n"
        f"- Status: **{result['independent_publisher_checksum_status']}**\n"
        f"- Downloaded bytes: {result.get('bytes_downloaded', 'not retrieved')}\n"
        f"- Source MD5: {result.get('publisher_metadata', {}).get('data_file_md5sum', 'unavailable')}\n"
        f"- Calculated MD5: {result.get('local_md5', 'unavailable')}\n"
        f"- Balanced Neff (illustrative only, not per-SNP): "
        f"{result.get('illustrative_balanced_effective_n_not_per_variant', 'unknown')}\n"
        f"- MAGMA/LDSC: NOT_RUN; source allele and LD mapping: HOLD\n"
    )
    (out / "checksum.md").write_text(content, encoding="utf-8")
    print(content)
    if os.getenv("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a", encoding="utf-8") as stream:
            stream.write(content)
    # Sources unavailable do not silently PASS, but retain evidence in artifact.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
