#!/usr/bin/env python3
"""Bounded verification of E-GEAD-1054 ZIP->TAR->GZIP cis-nominal source.

Only 256KiB of original HTTP body is requested (206 mandatory). No bulk download,
no genotype access, and NO claim that rs671 or its 3 AIS CS proxies are present.
"""
import argparse
import csv
import hashlib
import io
import json
import re
import struct
import tarfile
import urllib.request
import zlib
from pathlib import Path

DATASET = "E-GEAD-1054"
BASE = "https://ddbj.nig.ac.jp/public/ddbj_database/gea/experiment/E-GEAD-1000/E-GEAD-1054/"
URL = BASE + "E-GEAD-1054.processed.zip"
EXPECTED_ZIP_SIZE = 39271327981
EXPECTED_ZIP_MD5 = "312f19beade7c98afdd19d8d3965c270"
EXPECTED_TAR_SIZE = 39259422720
EXPECTED_TAR_MD5 = "c80a2f484906597bb2da5b7966fa8cf6"
REQUIRED_COLUMNS = (
    "phenotype_id", "gene", "variant_id", "tss_distance",
    "pval_nominal", "slope", "slope_se", "af", "ma_samples", "ma_count"
)
RANGE_BYTES = 262144
MAX_ZIP_INFLATE = 1048576
MAX_GZIP_INFLATE = 300000
USER_AGENT = "IS-G0022-small-source-metadata-audit/1.0"


def safe_extract_prefix(blob):
    """Decode only the first ZIP compressed member, one TAR header and GZIP head."""
    if len(blob) > RANGE_BYTES or len(blob) < 4096:
        raise ValueError("Unexpected HTTP byte-range prefix size")
    if blob[:4] != b"PK\x03\x04":
        raise ValueError("ZIP local file header not found")
    flag, zip_method = struct.unpack_from("<HH", blob, 6)
    name_len, extra_len = struct.unpack_from("<HH", blob, 26)
    member_name = blob[30:30 + name_len].decode("utf-8")
    start = 30 + name_len + extra_len
    if zip_method != 8 or not member_name.endswith(".tar"):
        raise ValueError("Expected one DEFLATEd TAR member, no direct region access")
    if start >= len(blob) - 4096:
        raise ValueError("No member prefix data")
    zip_inflater = zlib.decompressobj(-15)
    partial_tar = zip_inflater.decompress(blob[start:], MAX_ZIP_INFLATE)
    if len(partial_tar) < 1024:
        raise ValueError("Cannot read first TAR entry")
    header = tarfile.TarInfo.frombuf(partial_tar[:512], encoding="utf-8", errors="replace")
    inner_name = header.name
    if not inner_name.endswith(".cis_nominal.txt.gz") or not re.match(r"^[A-Za-z0-9_.-]+$", inner_name):
        raise ValueError("First TAR entry is not expected cis_nominal .txt.gz")
    if header.size <= 500:
        raise ValueError("First TAR entry too small")
    inner_gzip_head = partial_tar[512:]
    if inner_gzip_head[:2] != b"\x1f\x8b":
        raise ValueError("Nested gzip header invalid")
    gzip_inflater = zlib.decompressobj(16 + zlib.MAX_WBITS)
    partial_text = gzip_inflater.decompress(inner_gzip_head, MAX_GZIP_INFLATE)
    text_data = partial_text.decode("utf-8")
    lines = text_data.splitlines()
    if not lines:
        raise ValueError("Missing nominal SNP-gene header")
    cols = lines[0].split("\t")
    if cols != list(REQUIRED_COLUMNS):
        raise ValueError(f"Unexpected nominal QTL schema: {cols}")
    complete_lines = lines[1:-1]  # trailing line might be an incomplete compressed prefix
    samples = []
    for line in complete_lines[:80]:
        fields = line.split("\t")
        if len(fields) != len(cols):
            raise ValueError("Incomplete nominal association line")
        item = dict(zip(cols, fields))
        pval = float(item["pval_nominal"])
        af = float(item["af"])
        if not (0 <= pval <= 1 and 0 <= af <= 1):
            raise ValueError("Invalid p or frequency in QTL rows")
        samples.append(item)
    if len(samples) < 5:
        raise ValueError("No enough nominal records in prefix for audit")
    n_non_sig = sum(float(x["pval_nominal"]) > .05 for x in samples)
    if n_non_sig == 0:
        raise ValueError("Cannot verify non-significant nominal tests in prefix")
    return {
        "zip_member_name": member_name,
        "zip_local_flags": flag,
        "zip_compression_method": zip_method,
        "tar_first_member_name": inner_name,
        "tar_first_member_compressed_gzip_size_bytes": header.size,
        "nominal_qtl_columns": cols,
        "complete_sample_rows_examined": len(samples),
        "sample_p_greater_than_0_05": n_non_sig,
        "sample_p_le_0_05": len(samples) - n_non_sig,
        "first_5_row_p": [float(x["pval_nominal"]) for x in samples[:5]],
        "first_5_row_chr_variants": [x["variant_id"] for x in samples[:5]],
        "sample_gene": samples[0]["gene"],
        "sample_contains_g0022_verified": False,
        "first_member_locus_completeness_verified": False,
        "all_40_celltypes_contents_verified": False,
        "scientific_state": "UNFILTERED_NOMINAL_PREFIX_CONFIRMED_ALDH2_REGION_UNCHECKED",
    }


def get_bounded_http_prefix(url=URL):
    head = urllib.request.Request(url, headers={"User-Agent": USER_AGENT}, method="HEAD")
    with urllib.request.urlopen(head, timeout=20) as response:
        size = int(response.headers.get("Content-Length", "0"))
        ranges = response.headers.get("Accept-Ranges", "")
        if response.status != 200:
            raise ValueError("Source HEAD failed")
    if size != EXPECTED_ZIP_SIZE or "bytes" not in ranges.lower():
        raise ValueError("Unexpected source size or no byte ranges: refuse download")
    ranged = urllib.request.Request(
        url, headers={"User-Agent": USER_AGENT, "Range": f"bytes=0-{RANGE_BYTES-1}"})
    with urllib.request.urlopen(ranged, timeout=30) as response:
        if response.status != 206:
            raise ValueError("No partial HTTP 206; refuse unbounded source response")
        if response.headers.get("Content-Range", "") != f"bytes 0-{RANGE_BYTES-1}/{size}":
            raise ValueError("Invalid HTTP Content-Range")
        prefix = response.read(RANGE_BYTES + 1)
        if len(prefix) != RANGE_BYTES:
            raise ValueError("HTTP response is not exactly bounded")
    return prefix, size


def build(args):
    if args.probe_remote:
        prefix, size = get_bounded_http_prefix()
    else:
        if not args.local_prefix:
            raise ValueError("Need --local-prefix fixture or --probe-remote")
        prefix = args.local_prefix.read_bytes()
        size = None
    audit = safe_extract_prefix(prefix)
    audit.update({
        "audit": "IS_G0022_EGEAD1054_BOUNDED_NOMINAL_SOURCE_AUDIT_V1",
        "dataset": DATASET,
        "source_url": URL,
        "source_archive_bytes": size if size is not None else "synthetic_fixture",
        "official_zip_md5_not_full_file_verified": EXPECTED_ZIP_MD5,
        "official_internal_tar_uncompressed_bytes": EXPECTED_TAR_SIZE,
        "official_internal_tar_md5_not_full_file_verified": EXPECTED_TAR_MD5,
        "body_bytes_downloaded": RANGE_BYTES if args.probe_remote else 0,
        "http_range_enforced": bool(args.probe_remote),
        "allele_schema_only_variant_id_ref_alt_not_explicit_effect_allele": True,
        "low_freq_exclusion_in_file_name": "MAF0.05",
        "qtl_genotype_ld_available": False,
        "eas_AIS_four_cs_variant_coverage_verified": None,
        "eas_AIS_four_cs_variant_coverage_status": "NOT_ASSESSED",
        "full_cis_coloc_ready": False,
        "scientific_warning": (
            "First B_Activated cell file starts with nominal rows with p>0.05; "
            "does not verify the remaining 39 cell-type files or the chr12 rs671 region. "
            "File name applies MAF0.05 filtering, and effect-allele conventions, "
            "study covariance and complete rs671 locus denominator require separate validation."
        ),
    })
    args.outdir.mkdir(exist_ok=True, parents=True)
    (args.outdir / "G0022_EGEAD1054_SOURCE_NOMINAL_PROBE.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + "\n")
    return audit


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", required=True, type=Path)
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--probe-remote", action="store_true")
    group.add_argument("--local-prefix", type=Path, help="Previously collected bounded bytes for offline tests")
    args = p.parse_args()
    print(json.dumps(build(args), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
