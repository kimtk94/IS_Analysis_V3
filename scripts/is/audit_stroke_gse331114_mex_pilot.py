#!/usr/bin/env python3
"""Bounded real 10x-MEX integrity pilot for human ischemic stroke GSE331114.

Select ONE AIS and ONE non-stroke sample from official GEO SOFT metadata;
download processed MEX files within strict byte caps; stream Matrix Market
header/body, barcodes and features; check shapes, nnz, coordinates and values.
No hypothesis testing, cluster labels or stroke-specific expression claims.
Never include unverified participant IDs or patient-level data as artifacts.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
from urllib.parse import urlparse
from urllib.request import Request, urlopen

from audit_stroke_versus_rosmap_sources import (
    geo_soft_url, download_metadata, parse_soft, derive_group
)

CAPS = {
    "matrix": 180 * 1024 * 1024,
    "barcodes": 35 * 1024 * 1024,
    "features": 20 * 1024 * 1024,
}
MAX_MATRIX_NONZERO = 100_000_000


def classify_file(url):
    name = urlparse(url).path.lower().rsplit("/", 1)[-1]
    if re.search(r"(matrix|count|mtx)\.mtx(\.gz)?$", name):
        return "matrix"
    if re.search(r"(barcodes|barcode)\.tsv(\.gz)?$", name):
        return "barcodes"
    if re.search(r"(features|genes)\.tsv(\.gz)?$", name):
        return "features"
    return None


def group_files(sample):
    out = {}
    for f in sample["files"]:
        if not f or f == "NONE":
            continue
        typ = classify_file(f)
        if typ is not None:
            if typ in out:
                raise ValueError("Multiple " + typ + " files for sample")
            out[typ] = f.replace("ftp://ftp.ncbi.nlm.nih.gov/", "https://ftp.ncbi.nlm.nih.gov/")
    return out


def download(url, target, limit):
    parsed = urlparse(url)
    if parsed.scheme != "https" or parsed.hostname != "ftp.ncbi.nlm.nih.gov":
        raise ValueError("Data are only fetched from NCBI GEO")
    req = Request(url, headers={"User-Agent": "IS-Analysis-V3-MEX-Structure-QC/1.0"})
    n = 0
    digest = hashlib.sha256()
    with urlopen(req, timeout=90) as response:
        final = urlparse(response.geturl())
        if final.scheme != "https" or final.hostname != "ftp.ncbi.nlm.nih.gov":
            raise ValueError("Data redirected off NCBI")
        size = response.headers.get("Content-Length", "")
        if size.isdigit() and int(size) > limit:
            raise ValueError("Source exceeds cap")
        with open(target, "wb") as stream:
            while True:
                part = response.read(1 << 20)
                if not part:
                    break
                n += len(part)
                if n > limit:
                    raise ValueError("Download cap exceeded")
                stream.write(part)
                digest.update(part)
    if size.isdigit() and n != int(size):
        raise ValueError("Transfer size mismatch")
    return {"downloaded_bytes": n, "sha256_local": digest.hexdigest()}


def iter_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, encoding="utf-8")


def count_rows(path, min_fields=1):
    n = 0
    with iter_text(path) as stream:
        for line in stream:
            if not line.strip():
                continue
            if len(line.rstrip("\n").split("\t")) < min_fields:
                raise ValueError("Malformed feature/barcode line")
            n += 1
            if n > 20_000_000:
                raise ValueError("Excessive TSV rows")
    return n


def inspect_mtx(path):
    """Sparse Matrix Market coordinate integer file, without building dense matrix."""
    with iter_text(path) as stream:
        banner = stream.readline().strip()
        if not banner.lower().startswith("%%matrixmarket matrix coordinate"):
            raise ValueError("Not MatrixMarket coordinate data")
        chunks = banner.split()
        if len(chunks) < 5 or chunks[3].lower() not in {"integer", "real"}:
            raise ValueError("Unsupported MatrixMarket value format")
        header = None
        for line in stream:
            if not line.startswith("%") and line.strip():
                header = line.strip()
                break
        if header is None:
            raise ValueError("Missing MatrixMarket dimension line")
        dims = header.split()
        if len(dims) != 3 or not all(x.isdigit() for x in dims):
            raise ValueError("Malformed MatrixMarket dimensions")
        n_features, n_barcodes, n_nnz = map(int, dims)
        if n_features <= 0 or n_barcodes <= 0 or n_nnz < 0 or n_nnz > MAX_MATRIX_NONZERO:
            raise ValueError("Unexpected dimensions")
        observed = 0
        total_nonzero = 0.0
        for line in stream:
            if not line.strip() or line.startswith("%"):
                continue
            values = line.split()
            if len(values) != 3:
                raise ValueError("Bad coordinate line")
            i, j = int(values[0]), int(values[1])
            val = float(values[2])
            if not (1 <= i <= n_features and 1 <= j <= n_barcodes and val >= 0 and val < float("inf")):
                raise ValueError("Invalid sparse coordinate or count")
            observed += 1
            total_nonzero += val
            if observed > n_nnz:
                raise ValueError("More sparse entries than claimed")
        if observed != n_nnz:
            raise ValueError("MatrixMarket truncated")
    return {"genes_rows": n_features, "barcodes_columns": n_barcodes,
            "nonzero_entries": observed, "sum_matrix_entries": total_nonzero,
            "matrix_value_type": chunks[3].lower(),
            "unit_of_analysis": "single_cell_columns_not_independent_donors"}


def inspect_sample(sample, group, basedir):
    paths = group_files(sample)
    record = {
        "group": group,
        "gsm": sample["accession"],
        "title": sample["title"],
        "status": "HOLD",
        "files_detected": sorted(paths),
        "matrix_downloaded": False,
        "cells_are_independent_patient_replicates": False,
        "differential_expression_analysis": "NOT_RUN",
    }
    if set(paths) != {"matrix", "features", "barcodes"}:
        record["status"] = "MISSING_REQUIRED_MEX_FILES"
        return record
    local = {}
    for part in ("matrix", "features", "barcodes"):
        target = basedir / (sample["accession"] + "_" + part + ".gz")
        try:
            evidence = download(paths[part], target, CAPS[part])
            record[part + "_file_bytes"] = evidence["downloaded_bytes"]
            record[part + "_sha256_local_not_publisher"] = evidence["sha256_local"]
            local[part] = target
        except (OSError, ValueError, TimeoutError) as exc:
            record["status"] = "TRANSFER_HOLD"
            record["error"] = type(exc).__name__ + ": " + str(exc)[:150]
            return record
    try:
        mx = inspect_mtx(local["matrix"])
        n_barcodes = count_rows(local["barcodes"])
        n_features = count_rows(local["features"])
        record.update(mx)
        record["barcodes_tsv_rows"] = n_barcodes
        record["features_tsv_rows"] = n_features
        record["matrix_downloaded"] = True
        record["status"] = (
            "MEX_STRUCTURE_PASS" if mx["genes_rows"] == n_features and
            mx["barcodes_columns"] == n_barcodes else "MEX_DIMENSION_MISMATCH"
        )
    except (OSError, ValueError, EOFError) as exc:
        record["status"] = "MEX_FORMAT_HOLD"
        record["error"] = type(exc).__name__ + ": " + str(exc)[:150]
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", default="artifacts/is_stroke_mex_pilot")
    args = parser.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    result = {
        "source": "GSE331114",
        "audit_time_utc": __import__("datetime").datetime.now(__import__("datetime").timezone.utc).isoformat(),
        "scope": "ONE_AIS_ONE_CONTROL_PROCESSED_MEX_STRUCTURE_QC",
        "full_cohort_analysis_completed": False,
        "true_celltype_labels": "NOT_QC_VERIFIED",
        "donor_independence": "NOT_VERIFIED",
        "DEG_analysis": "NOT_RUN",
        "cases_and_controls": [],
    }
    try:
        text, _ = download_metadata(geo_soft_url("GSE331114"))
        info, samples = parse_soft(text)
        if info.get("accession") != "GSE331114":
            raise ValueError("Wrong GEO series")
        grouped = {"AIS": [], "Control": []}
        for sample in samples.values():
            group = derive_group(sample)
            if group in grouped:
                grouped[group].append(sample)
        if len(grouped["AIS"]) != 4 or len(grouped["Control"]) != 4:
            raise ValueError("Original 4+4 labels differ from verified metadata")
        with tempfile.TemporaryDirectory(prefix="stroke_mex_") as folder:
            for group in ("AIS", "Control"):
                sample = sorted(grouped[group], key=lambda x: x["accession"])[0]
                result["cases_and_controls"].append(inspect_sample(sample, group, Path(folder)))
    except (OSError, ValueError, TimeoutError, EOFError) as exc:
        result["fatal_error"] = type(exc).__name__ + ": " + str(exc)[:180]
    (out / "mex_pilot.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    lines = ["# Human ischemic stroke GSE331114 pilot MEX integrity", "",
             "This is **one AIS / one non-stroke sample**, not a complete n=4+4 expression analysis.",
             "", "| Group | GSM | Outcome | Genes | Cells | Nonzero |", "|---|---|---|---:|---:|---:|"]
    for x in result["cases_and_controls"]:
        lines.append(f"| {x['group']} | {x['gsm']} | {x['status']} | "
                     f"{x.get('genes_rows', '-')} | {x.get('barcodes_columns', '-')} | "
                     f"{x.get('nonzero_entries', '-')} |")
    lines += ["", "**Cells are not independent donor replicates. No differential expression, eQTL, or stroke-causal inference performed.**"]
    summary = "\n".join(lines) + "\n"
    (out / "summary.md").write_text(summary, encoding="utf-8")
    print(summary)
    if os.getenv("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a", encoding="utf-8") as fd:
            fd.write(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
