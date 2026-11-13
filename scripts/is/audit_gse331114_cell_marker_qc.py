#!/usr/bin/env python3
"""Cell-quality and exploratory immune-marker co-detection for GSE331114 (8 human samples).

Only computes per-sample diagnostics from original public GEO 10x MEX. No
reference-derived cell-type annotations exist in the input. "Marker support"
is NOT a cell-type call, clustering, pseudobulk DE, or an eQTL result.

Input files are downloaded only to ephemeral GitHub Actions worker disk.
All output is aggregated over each sample (no barcode-level matrix exported).
"""
from __future__ import annotations
from array import array
import argparse
from collections import Counter
import csv
import json
import math
import os
from pathlib import Path
import statistics
import tempfile
from datetime import datetime, timezone

from audit_stroke_versus_rosmap_sources import (
    geo_soft_url, download_metadata, parse_soft, derive_group,
)
from audit_stroke_gse331114_mex_pilot import (
    group_files, download, CAPS, iter_text,
)

SOURCE_RUN_ID = 43010007530
PANELS = {
    "T_cell_marker_support": ("CD3D", "CD3E", "TRAC"),
    "B_cell_marker_support": ("CD79A", "CD79B", "MS4A1"),
    "NK_marker_support": ("NKG7", "GNLY", "KLRD1"),
    "myeloid_marker_support": ("LST1", "LYZ", "FCN1", "CD14"),
    "microglia_marker_support": ("P2RY12", "TMEM119", "SALL1"),
}
QC_MIN_GENES = 200
QC_MAX_GENES = 7500
QC_MAX_MITO_PCT = 20.0
MAX_ENTRIES = 100_000_000


def feature_map(features_path):
    mt_rows = set()
    marker_rows = {}
    features = 0
    all_symbols = set()
    for panel in PANELS:
        for symbol in PANELS[panel]:
            if symbol in all_symbols:
                raise ValueError("Marker symbol reused between panels: " + symbol)
            all_symbols.add(symbol)
    counts = Counter()
    with iter_text(features_path) as f:
        for row in f:
            if not row.strip():
                continue
            fields = row.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError("Malformed features TSV")
            features += 1
            symbol = fields[1]
            if symbol.startswith("MT-"):
                mt_rows.add(features)
            if symbol in all_symbols:
                counts[symbol] += 1
                if counts[symbol] > 1:
                    raise ValueError("Duplicated marker gene symbol: " + symbol)
                marker_rows[features] = symbol
    if features != 36601:
        raise ValueError("Unexpected MEX features: " + str(features))
    if len(mt_rows) < 10:
        raise ValueError("Mitochondrial feature annotation unexpectedly absent")
    return features, mt_rows, marker_rows, counts


def stream_cell_qc(matrix_path, nfeatures, mitochondrial_rows, marker_rows):
    with iter_text(matrix_path) as stream:
        first = stream.readline().strip().lower()
        if first != "%%matrixmarket matrix coordinate integer general":
            raise ValueError("Expected unmodified integer 10x MEX coordinate matrix")
        line = next((l for l in stream if l.strip() and not l.startswith("%")), None)
        if line is None:
            raise ValueError("Missing matrix dimensions")
        dims = line.split()
        if len(dims) != 3:
            raise ValueError("Invalid 10x matrix header")
        nrows, ncells, expected = (int(v) for v in dims)
        if nrows != nfeatures or ncells < 1 or ncells > 300000 or expected > MAX_ENTRIES:
            raise ValueError("Unexpected matrix dimensions")
        total = array("Q", [0]) * ncells
        detected = array("I", [0]) * ncells
        mito = array("Q", [0]) * ncells
        markers = {symbol: bytearray(ncells) for symbol in marker_rows.values()}
        observed = 0
        total_umis = 0
        for l in stream:
            if not l.strip() or l.startswith("%"):
                continue
            parts = l.split()
            if len(parts) != 3:
                raise ValueError("Invalid sparse entry")
            row, col, count = (int(v) for v in parts)
            if not 1 <= row <= nrows or not 1 <= col <= ncells or count <= 0:
                raise ValueError("Invalid row, barcode column or nonpositive UMI")
            observed += 1
            if observed > expected:
                raise ValueError("Sparse matrix has more entries than header")
            col -= 1
            detected[col] += 1
            total[col] += count
            total_umis += count
            if row in mitochondrial_rows:
                mito[col] += count
            if row in marker_rows:
                markers[marker_rows[row]][col] = 1
        if observed != expected:
            raise ValueError("Matrix truncated before declared nnz")
    return {
        "nfeatures": nrows,
        "cells": ncells,
        "nnz": observed,
        "total_umi": total_umis,
        "total_umi_by_cell": total,
        "detected_features_by_cell": detected,
        "mitochondrial_umi_by_cell": mito,
        "marker_detected_by_cell": markers,
    }


def classify_marker_support(cell_index, marker_detected):
    """Require >=2 distinct markers per panel, and EXACTLY one supported panel.

    Output is provisional marker evidence, NOT validated identity.
    """
    supported = []
    for panel, genes in PANELS.items():
        found = sum(marker_detected.get(gene, bytearray())[cell_index] for gene in genes
                    if gene in marker_detected)
        if found >= 2:
            supported.append(panel)
    return supported[0] if len(supported) == 1 else ("MULTIPLE_PANELS" if supported else "UNASSIGNED")


def aggregate_cell_qc(qc, expected_barcodes):
    if qc["cells"] != expected_barcodes:
        raise ValueError("Barcode count and matrix dimensions disagree")
    features, totals, mt = (qc["detected_features_by_cell"],
                            qc["total_umi_by_cell"], qc["mitochondrial_umi_by_cell"])
    markers = qc["marker_detected_by_cell"]
    n = qc["cells"]
    percent_mt = [100.0 * mt[i] / totals[i] if totals[i] > 0 else 100.0 for i in range(n)]
    pass_genes = [QC_MIN_GENES <= features[i] <= QC_MAX_GENES for i in range(n)]
    pass_mt = [percent_mt[i] < QC_MAX_MITO_PCT for i in range(n)]
    keep = [pass_genes[i] and pass_mt[i] and totals[i] > 0 for i in range(n)]
    labels = Counter(classify_marker_support(i, markers) for i in range(n) if keep[i])
    n_pass = sum(keep)
    result = {
        "total_barcodes": n,
        "cells_pass_qc": n_pass,
        "cells_fail_qc": n - n_pass,
        "pct_pass_qc": round(100 * n_pass / n, 3),
        "qc_thresholds_exploratory": {
            "min_n_genes": QC_MIN_GENES, "max_n_genes": QC_MAX_GENES,
            "mito_percent_lt": QC_MAX_MITO_PCT,
        },
        "failed_gene_count_filter": sum(not x for x in pass_genes),
        "failed_mito_filter": sum(not x for x in pass_mt),
        "median_n_genes_all_barcodes": statistics.median(features),
        "median_n_umis_all_barcodes": statistics.median(totals),
        "median_mito_percent_all_barcodes": round(statistics.median(percent_mt), 3),
        "marker_panel_support_unvalidated_counts_after_QC": dict(labels),
        "marker_panel_support_unvalidated_fraction_after_QC": {
            label: round(count / n_pass, 6) if n_pass else None
            for label, count in labels.items()
        },
        "not_a_celltype_annotation": True,
        "requires_reference_validation_before_pseudobulk": True,
        "no_p_values_or_degs": True,
    }
    if n_pass == 0:
        result["status"] = "NO_CELLS_PASS_EXPLORATORY_QC"
    else:
        result["status"] = "SAMPLE_CELL_QC_AND_MARKER_SUPPORT_PASS"
    return result


def run(group, index, original_root, out):
    raw, _ = download_metadata(geo_soft_url("GSE331114"))
    info, samples = parse_soft(raw)
    if info.get("accession") != "GSE331114":
        raise ValueError("GEO series mismatch")
    group_samples = sorted((x for x in samples.values() if derive_group(x) == group),
                           key=lambda s: s["accession"])
    if len(group_samples) != 4 or index not in range(4):
        raise ValueError("Expected 4 AIS and 4 control records")
    sample = group_samples[index]
    gsm = sample["accession"]
    folder = original_root / f"is-gse331114-{group}-{index}"
    prior = json.loads((folder / "gene_counts.json").read_text(encoding="utf-8"))
    if prior["status"] != "SAMPLE_MEX_GENE_QC_PASS" or prior["accession"] != gsm:
        raise ValueError("Previous verified source evidence unavailable or mismatched")
    paths = group_files(sample)
    if set(paths) != {"matrix", "features", "barcodes"}:
        raise ValueError("Original GEO record lacks complete 10x MEX files")
    with tempfile.TemporaryDirectory(prefix="is_cell_marker_") as tmp:
        base = Path(tmp)
        downloaded = {}
        fingerprints = {}
        for part in ("features", "barcodes", "matrix"):
            dest = base / f"{part}.gz"
            fingerprint = download(paths[part], dest, CAPS[part])
            downloaded[part] = dest
            fingerprints[part] = fingerprint["sha256_local"]
            # Extra verification against *identical original processed files* previously checked.
            if prior.get(part + "_sha256_not_source_verified") != fingerprint["sha256_local"]:
                raise ValueError("GEO file SHA256 changed from prior verified source: " + part)
        nfeatures, mt_rows, marker_rows, counts = feature_map(downloaded["features"])
        qc = stream_cell_qc(downloaded["matrix"], nfeatures, mt_rows, marker_rows)
        with iter_text(downloaded["barcodes"]) as f:
            nbarcodes = sum(1 for line in f if line.strip())
        if prior["cell_barcodes"] != nbarcodes or prior["genes"] != nfeatures:
            raise ValueError("Prior sample cell/feature dimensions don't match")
        if prior["total_umi_count"] != qc["total_umi"]:
            raise ValueError("Prior sample UMI total mismatch")
        summary = aggregate_cell_qc(qc, nbarcodes)
    summary.update({
        "source":"GSE331114",
        "GSM":gsm, "group":group, "index":index,
        "original_qc_run_id":SOURCE_RUN_ID,
        "source_file_sha256_matches_prior":True,
        "source_celltype_annotation": "NOT_AVAILABLE_IN_10X_MEX",
        "donor_independence": "NOT_VERIFIED_WITH_INDIVIDUAL_IDENTIFIERS",
        "celltype_marker_scoring": "TWO_DISTINCT_GENES_PER_LINEAGE_AND_EXCLUSIVE_SUPPORT",
        "subset": "CD45_ENRICHED_HUMAN_BRAIN",
        "celltype_pseudobulk": "NOT_PERFORMED",
        "statistical_testing": "NOT_PERFORMED",
        "marker_genes_observed": sorted(marker_rows.values()),
        "marker_genes_missing": sorted({g for panel in PANELS.values() for g in panel}
                                       - set(marker_rows.values())),
        "metadata_characteristic_field_names": sorted({
            x.split(":",1)[0].strip().lower()
            for x in sample.get("characteristics", []) if ":" in x
        }),
    })
    out.mkdir(parents=True, exist_ok=True)
    (out/"marker_qc.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    lines=[
        f"# {gsm} — {group} 10x cell-level QC & marker co-detection",
        "",
        f"**{summary['status']}** — source files SHA256 matched prior QC run, {summary['total_barcodes']} original barcodes.",
        "",
        f"QC-pass barcodes: **{summary['cells_pass_qc']}** "
        f"({summary['pct_pass_qc']}%); median mitochondrial percentage: "
        f"{summary['median_mito_percent_all_barcodes']}%.",
        "",
        "| Exclusive exploratory marker support | QC-passing barcodes |",
        "|---|---:|",
    ]
    for label in sorted(summary["marker_panel_support_unvalidated_counts_after_QC"]):
        lines.append(f"| {label} | {summary['marker_panel_support_unvalidated_counts_after_QC'][label]} |")
    lines += [
        "",
        "**Interpretation HOLD:** marker support is NOT validated cell type; single-cell barcodes",
        "are NOT patient replicates. Formal donor × cell-type pseudobulk is not done.",
    ]
    md="\n".join(lines)+"\n"
    (out/"summary.md").write_text(md, encoding="utf-8")
    print(md)
    if os.getenv("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"],"a",encoding="utf-8") as f:
            f.write(md)
    return summary


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--group", choices=("AIS","Control"), required=True)
    p.add_argument("--index", choices=range(4), type=int, required=True)
    p.add_argument("--source-artifacts", default="downloaded_gse331114")
    p.add_argument("--out-dir", default="artifacts/is_gse331114_cell_marker_qc")
    args=p.parse_args()
    run(args.group,args.index,Path(args.source_artifacts),Path(args.out_dir))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
