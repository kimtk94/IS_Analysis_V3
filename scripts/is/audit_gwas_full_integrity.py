#!/usr/bin/env python3
"""GIGASTROKE full-original-file transfer and gzip/row QC on ephemeral CI storage.

No patient-level data, no raw summary statistics in Git or Actions artifacts.
Records independently calculated file SHA-256, NOT an independently verified
publisher checksum. Does NOT imply phenotype cohort audit, MAGMA/LDSC, eQTL
or colocalization scientific readiness.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import tempfile
import urllib.error
import urllib.request

FILE_ROOT = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
    "GCST90104001-GCST90105000"
)
SOURCES = {
    "EUR": {
        "accession": "GCST90104540",
        "remote_size_at_prefix_audit_bytes": 283482604,
        "n_cases_published": 62100, "n_controls_published": 1234808
    },
    "EAS": {
        "accession": "GCST90104545",
        "remote_size_at_prefix_audit_bytes": 260689208,
        "n_cases_published": 19032, "n_controls_published": 237242
    },
}
REQUIRED = [
    "chromosome", "base_pair_location", "effect_allele_frequency",
    "beta", "standard_error", "p_value", "effect_allele", "other_allele",
]
MAX_FILE_BYTES = 650 * 1024 * 1024
MIN_VARIANTS = 100_000
CHUNK_BYTES = 4 * 1024 * 1024


def full_url(accession):
    return f"{FILE_ROOT}/{accession}/{accession}_buildGRCh37.tsv.gz"


def read_numeric(row, name):
    x = float(row[name])
    if not math.isfinite(x):
        raise ValueError(name + " nonfinite")
    return x


def validate_row(row):
    """Return a single primary failure class or None. No raw variant data logged."""
    try:
        chrom = row["chromosome"].strip().replace("chr", "").upper()
        if chrom not in {str(x) for x in range(1, 23)}:
            return "CHROMOSOME_OUTSIDE_AUTOSOMES"
        pos = int(float(row["base_pair_location"]))
        if pos < 1 or pos > 260_000_000:
            return "INVALID_POSITION"
        p = read_numeric(row, "p_value")
        if not 0 <= p <= 1:
            return "INVALID_P_VALUE"
        se = read_numeric(row, "standard_error")
        if se <= 0:
            return "INVALID_SE"
        read_numeric(row, "beta")
        eaf = read_numeric(row, "effect_allele_frequency")
        if not 0 <= eaf <= 1:
            return "INVALID_FREQUENCY"
        a1 = row["effect_allele"].strip().upper()
        a2 = row["other_allele"].strip().upper()
        if a1 not in {"A", "C", "G", "T"} or a2 not in {"A", "C", "G", "T"} or a1 == a2:
            return "INVALID_ALLELES"
    except (TypeError, ValueError, KeyError, OverflowError):
        return "INVALID_OR_MISSING_VALUE"
    return None


def stream_download(url, destination):
    request = urllib.request.Request(
        url, headers={"User-Agent": "IS-Analysis-V3-FullGwas-IntegrityAudit/1.0"}
    )
    result = {"url": url, "transfer_status": "SOURCE_UNAVAILABLE"}
    with urllib.request.urlopen(request, timeout=60) as resp:
        final = resp.geturl()
        if not final.startswith(FILE_ROOT + "/"):
            raise ValueError("Unexpected cross-host or prefix redirect")
        declared = resp.headers.get("content-length")
        digest = hashlib.sha256()
        size = 0
        with open(destination, "wb") as out:
            while True:
                chunk = resp.read(CHUNK_BYTES)
                if not chunk:
                    break
                size += len(chunk)
                if size > MAX_FILE_BYTES:
                    raise RuntimeError("Source exceeds bounded file cap")
                digest.update(chunk)
                out.write(chunk)
        result.update({
            "http_status": resp.status,
            "bytes_downloaded": size,
            "content_length_declared": int(declared) if declared and declared.isdigit() else None,
            "sha256_calculated_local_not_publisher_verified": digest.hexdigest(),
        })
        if resp.status != 200 or size < 1000:
            result["transfer_status"] = "HTTP_OR_EMPTY_FILE"
        elif declared and declared.isdigit() and size != int(declared):
            result["transfer_status"] = "TRUNCATED_TRANSFER"
        else:
            result["transfer_status"] = "COMPLETE_TRANSFER_SIZE_CHECK"
    return result


def scan_source_gzip(path):
    count = 0
    failures = Counter()
    autosome_count = Counter()
    below_5e8 = 0
    zero_p = 0
    palindromic = 0
    excluded_nonautosomal = 0
    chi2_sum = 0.0
    chi2_n = 0
    fields = []
    with gzip.open(path, "rt", newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = list(reader.fieldnames or [])
        missing = sorted(set(REQUIRED) - set(fields))
        if missing:
            return {
                "full_gzip_crc_integrity": "UNVERIFIED_EARLY_STOP",
                "row_qc_status": "MISSING_FIELDS",
                "missing_headers": missing,
                "header_columns": fields,
            }
        for row in reader:
            count += 1
            error = validate_row(row)
            if error == "CHROMOSOME_OUTSIDE_AUTOSOMES":
                excluded_nonautosomal += 1
                continue
            if error:
                failures[error] += 1
                continue
            chrom = row["chromosome"].strip().replace("chr", "")
            autosome_count[chrom] += 1
            p = float(row["p_value"])
            below_5e8 += p < 5e-8
            zero_p += p == 0
            a1 = row["effect_allele"].strip().upper()
            a2 = row["other_allele"].strip().upper()
            if frozenset((a1, a2)) in (frozenset(("A", "T")), frozenset(("C", "G"))):
                palindromic += 1
            chi2 = (float(row["beta"]) / float(row["standard_error"])) ** 2
            if math.isfinite(chi2):
                chi2_sum += chi2
                chi2_n += 1
    # Exhausting the gzip iterator verifies its EOF/CRC. Not an external hash match.
    valid = count - excluded_nonautosomal - sum(failures.values())
    bad_ratio = sum(failures.values()) / max(count - excluded_nonautosomal, 1)
    status = ("STRUCTURE_QC_PASS" if count >= MIN_VARIANTS and bad_ratio < 0.001
              else "STRUCTURE_QC_REVIEW")
    return {
        "full_gzip_crc_integrity": "PASS",
        "row_qc_status": status,
        "total_variant_rows": count,
        "valid_autosomal_snp_rows": valid,
        "excluded_nonautosomal_rows_not_invalid": excluded_nonautosomal,
        "invalid_ratio": round(bad_ratio, 8),
        "invalid_row_counts": dict(failures),
        "n_genome_wide_significant_p_below_5e8": below_5e8,
        "p_value_zero_count_underflow_possible": zero_p,
        "palindromic_variant_rows": palindromic,
        "mean_beta_se_squared_autosomal_not_ldsc_h2": round(chi2_sum / max(chi2_n, 1), 5),
        "chr_autosome_counts": dict(sorted(autosome_count.items(), key=lambda x: int(x[0]))),
        "header_columns": fields,
        "rsid_column_present": bool(set(fields) & {"rsid", "rs_id", "SNP"}),
        "per_variant_sample_size_column_present": bool(
            set(fields) & {"n", "N", "sample_size", "n_total", "n_cases", "n_controls"}
        ),
        "source_independent_sha256_verified": False,
        "snps_id_and_ld_reference_qc": "HOLD",
        "magma_analysis": "NOT_RUN",
        "ldsc_seg_analysis": "NOT_RUN",
    }


def audit_source(ancestry, source, workdir):
    accession = source["accession"]
    dest = workdir / f"{accession}.tsv.gz"
    result = {
        "ancestry": ancestry, "accession": accession,
        "source_url": full_url(accession),
        "as_of": datetime.now(timezone.utc).isoformat(),
        "published_case_count_unverified_in_sumstats": source["n_cases_published"],
        "published_control_count_unverified_in_sumstats": source["n_controls_published"],
        "genome_build_from_filename": "GRCh37",
        "full_source_status": "HOLD",
        "science_status": "NOT_RUN",
    }
    try:
        result.update(stream_download(result["source_url"], dest))
        if result["transfer_status"] != "COMPLETE_TRANSFER_SIZE_CHECK":
            return result
        result["prior_head_size_bytes"] = source["remote_size_at_prefix_audit_bytes"]
        result["prior_head_size_matches"] = (
            result["bytes_downloaded"] == source["remote_size_at_prefix_audit_bytes"]
        )
        result.update(scan_source_gzip(dest))
        if result.get("full_gzip_crc_integrity") == "PASS" and result.get("row_qc_status") == "STRUCTURE_QC_PASS":
            result["full_source_status"] = "TRANSFER_AND_STRUCTURE_QC_PASS"
            # Still source integrity not independently anchored and phenotype not confirmed.
            result["research_g1"] = "PARTIAL_HOLD"
        else:
            result["full_source_status"] = "QC_REVIEW"
            result["research_g1"] = "HOLD"
    except (OSError, RuntimeError, ValueError, urllib.error.URLError, EOFError) as exc:
        result["error_type"] = type(exc).__name__
        result["error_short"] = str(exc)[:220]
        result["full_source_status"] = "SOURCE_UNAVAILABLE_OR_CORRUPT"
        result["research_g1"] = "HOLD"
    finally:
        try:
            dest.unlink(missing_ok=True)
        except OSError:
            pass
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out-dir", default="artifacts/is_eas_gwas_full_qc")
    p.add_argument("--ancestries", choices=["EUR", "EAS", "both"], default="both")
    args = p.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    requested = SOURCES.items() if args.ancestries == "both" else [(args.ancestries, SOURCES[args.ancestries])]
    with tempfile.TemporaryDirectory(prefix="is_gwas_full_") as tmp:
        results = [audit_source(anc, item, Path(tmp)) for anc, item in requested]
    document = {
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "evidence_type": "FULL_TRANSFER_LOCAL_SHA256_AND_GZIP_STRUCTURE_QC_NOT_PUBLISHER_HASH_VERIFICATION",
        "research_ready": False,
        "results": results,
    }
    (out / "source_full_qc.json").write_text(
        json.dumps(document, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    lines = [
        "# GIGASTROKE EUR/EAS full original-file QC",
        "", "| Ancestry | Complete transfer & CRC | Valid rows | Local SHA-256 | G1 |",
        "|---|---|---:|---|---|",
    ]
    for result in results:
        lines.append(
            f"| {result['ancestry']} | {result.get('full_source_status', 'HOLD')} | "
            f"{result.get('valid_autosomal_snp_rows', 'unknown')} | "
            f"{result.get('sha256_calculated_local_not_publisher_verified', 'missing')[:16]}… | "
            f"{result.get('research_g1', 'HOLD')} |"
        )
    lines.extend(["", "Calculated SHA-256 is a reproducible local file fingerprint, not comparison to a publisher checksum.",
                  "No rsID or per-SNP N in sampled headers; population-specific LD mapping, cohort ascertainment,",
                  "and the sample-size model remain HOLD. No MAGMA or LDSC-SEG has been run."])
    summary = "\n".join(lines) + "\n"
    (out / "source_full_qc.md").write_text(summary, encoding="utf-8")
    print(summary)
    step_file = os.getenv("GITHUB_STEP_SUMMARY")
    if step_file:
        with open(step_file, "a", encoding="utf-8") as stream:
            stream.write(summary)
    # Fail closed in the evidence, while retaining diagnostics if a remote host fails.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
