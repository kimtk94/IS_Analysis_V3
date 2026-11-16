#!/usr/bin/env python3
"""GIANT 2018 corrected EUR BMI original-source streaming QC (no raw retention).

A matching gzip/Content-Length and complete structural scan is *not* LDSC-ready:
no consortium-provided reference checksum, duplicate-variant, imputation, cohort
overlap, or reference-LD validation has been established by this script.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import math
import re
import urllib.request
from collections import Counter
from pathlib import Path
from urllib.parse import urlsplit

SOURCE = ("https://giant-consortium.web.broadinstitute.org/images/c/c8/"
          "Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz")
EXPECTED = ("CHR", "POS", "SNP", "Tested_Allele", "Other_Allele",
            "Freq_Tested_Allele_in_HRS", "BETA", "SE", "P", "N")
MAX_BYTES = 75_000_000
MAX_ROWS = 12_000_000
SAMPLE_ERR_LIMIT = 8
ALLOWED_HOST = "giant-consortium.web.broadinstitute.org"


class BoundedHashReader:
    """Incrementally digest compressed bytes without saving the source."""
    def __init__(self, response, cap):
        self.response = response
        self.cap = cap
        self.count = 0
        self.sha256 = hashlib.sha256()
        self.md5 = hashlib.md5(usedforsecurity=False)

    def read(self, n=-1):
        # gzip.GzipFile requests bounded chunks, but enforce an upper bound.
        n = 1_048_576 if n is None or n < 0 else min(n, 1_048_576)
        block = self.response.read(n)
        self.count += len(block)
        if self.count > self.cap:
            raise ValueError("compressed transfer exceeds byte cap")
        self.sha256.update(block)
        self.md5.update(block)
        return block


def pinned_source(manifest):
    with Path(manifest).open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    hits = [r for r in rows if r.get("trait") == "BMI" and
            r.get("study") == "GIANT_Yengo2018_corrected"]
    if len(hits) != 1:
        raise ValueError("expected one pinned BMI source row")
    row = hits[0]
    if (row.get("exact_original_url") != SOURCE or
            row.get("approved_as_ldsc_input") != "NO" or
            row.get("genome_build") != "hg19_GRCh37"):
        raise ValueError("original source/build/approval contract mismatch")
    return row


def numeric_cell(value):
    v = float(value.strip())
    if not math.isfinite(v):
        raise ValueError("nonfinite numeric value")
    return v


def scan(response, *, max_bytes=MAX_BYTES, max_rows=MAX_ROWS):
    status = getattr(response, "status", 200)
    if status != 200:
        raise ValueError("full scan requires HTTP 200 (no partial content)")
    headers = response.headers
    size_raw = headers.get("Content-Length")
    if size_raw and (not size_raw.isdigit() or int(size_raw) > max_bytes):
        raise ValueError("publisher compressed size exceeds cap")
    reader = BoundedHashReader(response, max_bytes)
    result = {
        "source_status": "HOLD",
        "rows_total": 0,
        "valid_autosome_rows": 0,
        "invalid_rows": 0,
        "non_autosomal_rows": 0,
        "p_underflow_zero_rows": 0,
        "chr_counts": {},
        "errors_by_type": {},
        "error_examples": [],
        "min_N": None,
        "max_N": None,
        "no_duplicated_variant_check": True,
        "no_original_publisher_checksum": True,
        "ready_for_ldsc": False,
        "scientific_status": "HOLD",
    }
    counters = Counter()
    chr_counts = Counter()
    max_n = None
    min_n = None

    with gzip.GzipFile(fileobj=reader, mode="rb") as gz:
        header_line = gz.readline(8193)
        if len(header_line) > 8192 or not header_line.endswith(b"\n"):
            raise ValueError("header missing or too large")
        header = header_line.decode("utf-8-sig", "strict").rstrip("\r\n").split("\t")
        if tuple(header) != EXPECTED:
            raise ValueError("original field names/order mismatch")
        for rownum, rawline in enumerate(gz, start=2):
            if rownum > max_rows + 1:
                raise ValueError("row cap exceeded")
            if len(rawline) > 8192:
                raise ValueError("unexpectedly long variant row")
            result["rows_total"] += 1
            reason = None
            try:
                fields = rawline.decode("utf-8", "strict").rstrip("\r\n").split("\t")
                if len(fields) != len(EXPECTED):
                    raise ValueError("field_count")
                chrom, pos_str, snp, allele_a, allele_b, freq, beta, se, p, n = fields
                if not re.fullmatch(r"[0-9]+", chrom):
                    raise ValueError("chromosome_noninteger")
                chrom_i = int(chrom)
                if chrom_i < 1 or chrom_i > 22:
                    result["non_autosomal_rows"] += 1
                    continue
                if not pos_str.isdigit() or not (0 < int(pos_str) <= 300_000_000):
                    raise ValueError("position_invalid")
                if not snp or len(snp) > 128:
                    raise ValueError("variant_identifier_invalid")
                if (allele_a not in "ACGT" or allele_b not in "ACGT"
                        or len(allele_a) != 1 or len(allele_b) != 1 or allele_a == allele_b):
                    raise ValueError("alleles_invalid")
                f = numeric_cell(freq)
                b = numeric_cell(beta)
                s = numeric_cell(se)
                pv = numeric_cell(p)
                sample = numeric_cell(n)
                if not (0 <= f <= 1):
                    raise ValueError("HRS_frequency_out_of_bounds")
                if not (s > 0):
                    raise ValueError("SE_nonpositive")
                if not (0 <= pv <= 1):
                    raise ValueError("P_out_of_bounds")
                if not (sample > 0):
                    raise ValueError("N_nonpositive")
                if pv == 0.0:
                    result["p_underflow_zero_rows"] += 1
                result["valid_autosome_rows"] += 1
                chr_counts[str(chrom_i)] += 1
                min_n = sample if min_n is None else min(min_n, sample)
                max_n = sample if max_n is None else max(max_n, sample)
            except (ValueError, UnicodeError, OverflowError) as exc:
                reason = str(exc)
                counters[reason] += 1
                result["invalid_rows"] += 1
                if len(result["error_examples"]) < SAMPLE_ERR_LIMIT:
                    result["error_examples"].append({"row": rownum, "issue": reason})
    if size_raw is not None and reader.count != int(size_raw):
        raise ValueError("truncated compressed body vs Content-Length")
    if result["rows_total"] == 0 or result["valid_autosome_rows"] == 0:
        raise ValueError("no usable variant records")
    result.update({
        "full_compressed_bytes_read": reader.count,
        "http_content_length": int(size_raw) if size_raw else None,
        "source_sha256": reader.sha256.hexdigest(),
        "source_md5_unverified": reader.md5.hexdigest(),
        "gzip_stream_completed": True,
        "original_header_verified": True,
        "min_N": min_n,
        "max_N": max_n,
        "chr_counts": dict(sorted(chr_counts.items(), key=lambda k: int(k[0]))),
        "errors_by_type": dict(counters),
    })
    result["source_status"] = (
        "PASS_STREAM_STRUCTURE_ONLY" if result["invalid_rows"] == 0 else "HOLD_INVALID_RECORDS"
    )
    return result


def audit_one(manifest, *, opener=urllib.request.urlopen):
    pinned_source(manifest)
    report = {
        "schema": "ckm-bmi-giant-original-full-stream-v1",
        "audit_timestamp_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "trait": "BMI",
        "source": SOURCE,
        "build": "hg19/GRCh37 (README claim, not independent coordinate liftover)",
        "source_approved": False,
        "ready_for_ldsc": False,
        "scientific_status": "HOLD",
        "provenance": "Corrected Yengo2018 / GIANT official distribution",
        "warning": "HRS frequency is not participant-wide EAF; UKB sample overlap not evaluated",
    }
    try:
        request = urllib.request.Request(
            SOURCE, headers={"User-Agent": "Academic-CKM-Stroke-GWAS-Integrity-Audit/1.0"}
        )
        with opener(request, timeout=60) as resp:
            url = urlsplit(resp.geturl())
            if url.scheme != "https" or url.hostname != ALLOWED_HOST:
                raise ValueError("source redirected off official host")
            report["result"] = scan(resp)
    except Exception as exc:
        report["result"] = {
            "source_status": "HOLD_SOURCE_ERROR",
            "error_type": type(exc).__name__, "error_detail": str(exc)[:220],
            "ready_for_ldsc": False, "scientific_status": "HOLD",
        }
    return report


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", type=Path, default=Path("config/ckm_stroke/p1_exact_source_endpoints.tsv"))
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    report = audit_one(args.manifest)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    summary = report["result"]
    print(json.dumps({"source_status": summary["source_status"],
                      "rows_total": summary.get("rows_total"),
                      "valid_autosome_rows": summary.get("valid_autosome_rows"),
                      "invalid_rows": summary.get("invalid_rows"),
                      "bytes": summary.get("full_compressed_bytes_read"),
                      "source_sha256": summary.get("source_sha256"),
                      "scientific_status": report["scientific_status"]}, ensure_ascii=False))
    if summary["source_status"] != "PASS_STREAM_STRUCTURE_ONLY":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
