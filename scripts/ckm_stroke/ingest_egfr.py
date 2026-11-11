#!/usr/bin/env python3
"""Real CKDGen 2019 EUR eGFR source download + full gzip/row integrity audit.

The official Wuttke2019 source page lists no source MD5 for this file.
SHA-256 is locally computed, NOT described as consortium-verified.
Raw and standardized GWAS are deliberately ephemeral.
"""
import argparse
import datetime as dt
import json
import os
from pathlib import Path
import tempfile
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import Request, urlopen

from scripts.ckm_stroke.ingest_uacr import hash_file, stream_to_file
from scripts.ckm_stroke.normalize_sumstats import normalize

EGFR = {
    "url": "https://ckdgen.imbi.uni-freiburg.de/files/Wuttke2019/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt.gz",
    "filename": "20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt.gz",
    "source_landing": "https://ckdgen.imbi.uni-freiburg.de/datasets/wuttke-li-y-li-m-sieber-feitosa-gorski-et-al-2019/",
    "source_readme": "https://ckdgen.imbi.uni-freiburg.de/files/Wuttke2019/README.txt",
}

def ingest(root, max_mb=1400, opener=urlopen, retry_sleep=time.sleep):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    limit = max_mb * 1024 * 1024
    url = EGFR["url"]
    target = root / EGFR["filename"]
    normalized = root / "eGFR_EUR_GRCh37.normalized.tsv.gz"
    agent = {"User-Agent": "CKM-Stroke-eGFR-Source-Audit/1.0"}
    expected_length = None
    try:
        with opener(Request(url, method="HEAD", headers=agent), timeout=45) as head:
            if urlparse(head.geturl()).hostname != "ckdgen.imbi.uni-freiburg.de":
                return {"status": "HOLD", "scientific_status": "NOT_RUN",
                        "reason": "unapproved download host"}
            reported = head.headers.get("Content-Length")
            if reported:
                expected_length = int(reported)
                if expected_length > limit:
                    return {"status": "HOLD", "scientific_status": "NOT_RUN",
                            "reason": "official source exceeds configured size limit",
                            "declared_bytes": expected_length}
    except (HTTPError, URLError, OSError, ValueError):
        expected_length = None
    last_error = None
    size = 0
    for attempt in range(1, 4):
        try:
            with opener(Request(url, headers=agent), timeout=150) as resp:
                if resp.status != 200:
                    raise OSError("non-200 download response " + str(resp.status))
                if urlparse(resp.geturl()).hostname != "ckdgen.imbi.uni-freiburg.de":
                    raise OSError("unapproved final download host")
                size = stream_to_file(resp, target, limit)
            last_error = None
            break
        except (HTTPError, URLError, OSError, ValueError) as exc:
            last_error = str(exc)[:240]
            target.unlink(missing_ok=True)
            if attempt < 3:
                retry_sleep(attempt)
    if last_error:
        return {"status": "HOLD", "scientific_status": "NOT_RUN",
                "reason": "download incomplete", "error": last_error}
    if expected_length is not None and size != expected_length:
        return {"status": "FAIL", "scientific_status": "NOT_RUN",
                "reason": "download differs from declared length",
                "declared": expected_length, "actual": size}
    source_sha = hash_file(target, "sha256")
    try:
        result = normalize(target, normalized, "ckdgen2019", source_sha)
    except (OSError, ValueError, EOFError, AssertionError) as exc:
        return {"status": "FAIL", "scientific_status": "NOT_RUN",
                "reason": "full gzip/normalization failure", "details": str(exc)[:320]}
    return {
        "status": "PASS" if result["status"] in {"PASS", "HOLD"} else "FAIL",
        "scientific_status": "HOLD",
        "source_provenance": "CKDGen_Wuttke2019_EUR",
        "source_checksum_status": "UNVERIFIED_OFFICIAL_CHECKSUM",
        "local_sha256_verified_against_official": False,
        "gzip_stream_read_to_end": True,
        "download_size_checked": expected_length is not None,
        "declared_bytes": expected_length,
        "source_bytes": size,
        "source_sha256": source_sha,
        "normalization_status": result["status"],
        "input_rows": result["input_rows"],
        "normalized_rows": result["output_rows"],
        "excluded_total": result["excluded_total"],
        "excluded": result["excluded"],
        "excluded_examples": result["excluded_examples"],
        "normalized_sha256": result["output_sha256"],
        "note": "Full gzip/row parse and locally computed SHA-256, without an independently published source hash. NOT a statistical or allele-harmonization PASS.",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-mb", type=int, default=1400)
    ap.add_argument("--out", default="artifacts/ckm_stroke/egfr_full_integrity.json")
    args = ap.parse_args()
    try:
        if not 1 <= args.max_mb <= 1500:
            raise ValueError("max-mb must be 1..1500")
        with tempfile.TemporaryDirectory(prefix="ckm_egfr_") as directory:
            result = ingest(directory, args.max_mb)
    except (OSError, ValueError, OverflowError) as exc:
        result = {"status": "FAIL", "scientific_status": "NOT_RUN", "reason": str(exc)}
    result.update({
        "trait": "EGFR", "ancestry": "EUR", "url": EGFR["url"],
        "landing": EGFR["source_landing"], "readme": EGFR["source_readme"],
        "timestamp_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source_commit": os.getenv("GITHUB_SHA", ""),
    })
    dest = Path(args.out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False), flush=True)
    summ = os.getenv("GITHUB_STEP_SUMMARY")
    if summ:
        with open(summ, "a", encoding="utf-8") as fh:
            fh.write("\n### CKDGen eGFR EUR full source evidence\n")
            for key in ("status", "scientific_status", "source_checksum_status",
                        "source_bytes", "input_rows", "normalized_rows",
                        "excluded_total", "normalization_status"):
                fh.write("- " + key + ": " + str(result.get(key, "not_available")) + "\n")
    return int(result["status"] == "FAIL")


if __name__ == "__main__":
    raise SystemExit(main())
