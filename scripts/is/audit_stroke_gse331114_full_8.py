#!/usr/bin/env python3
"""Full GSE331114 per-GSM 10x gene count QC and prespecified IS-gene summary.

One real GEO sample per GitHub Actions matrix job (8 jobs, parallel).
All comparisons are descriptive; no donor-level differential expression test,
cell-type labels, batch correction, eQTL or causal analysis. Sparse counts
are processed via streaming and discarded from runner afterwards.

The output is *aggregate public sample gene count evidence*, not individual
cell barcodes/count matrices.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import gzip
import json
import os
from pathlib import Path
import tempfile

from audit_stroke_versus_rosmap_sources import (
    download_metadata, geo_soft_url, parse_soft, derive_group
)
from audit_stroke_gse331114_mex_pilot import (
    group_files, download, CAPS, iter_text
)

GENES = ("ABO", "ALDH2", "CD40", "F11", "FGF5", "FURIN", "GRK5", "IL6R",
         "KIAA0319", "MMP12", "PROCR", "SCARA5", "SWAP70", "TFPI", "TMPRSS5")
LEGACY_ALIAS = {"IL6RA": "IL6R"}
OTHER_CONTEXT = ("PTPRC", "AIF1", "TYROBP", "CD3D", "CD14", "MS4A1",
                 "CX3CR1", "TREM2", "PECAM1", "GFAP", "RBFOX3")
TARGETS = tuple(dict.fromkeys((*GENES, *OTHER_CONTEXT)))


def parse_features(path):
    ids = {}
    row_count = 0
    duplicated_symbols = set()
    mapping = defaultdict(list)
    with iter_text(path) as inp:
        for line in inp:
            if not line.strip():
                continue
            row_count += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError("Missing features symbol")
            if fields[1] in TARGETS:
                mapping[fields[1]].append(row_count)
            if row_count > 100000:
                raise ValueError("Unexpectedly many gene rows")
    for name, positions in mapping.items():
        if len(positions) > 1:
            duplicated_symbols.add(name)
    return row_count, mapping, sorted(duplicated_symbols)


def stream_gene_counts(matrix_path, nfeatures, selected_rows):
    totals = defaultdict(float)
    detected_cells = defaultdict(int)
    parsed_records = 0
    counts_total = 0.0
    with iter_text(matrix_path) as stream:
        banner = stream.readline().strip().lower()
        if not banner.startswith("%%matrixmarket matrix coordinate"):
            raise ValueError("Unexpected sparse matrix type")
        parts = banner.split()
        if len(parts) < 5 or parts[3] not in {"integer", "real"}:
            raise ValueError("Invalid MatrixMarket value type")
        dims = None
        for line in stream:
            if not line.startswith("%") and line.strip():
                dims = line.split()
                break
        if not dims or len(dims) != 3:
            raise ValueError("No matrix dimensions")
        n_genes, cells, expected_entries = (int(x) for x in dims)
        if n_genes != nfeatures or not 0 < cells <= 300000 or expected_entries > 100000000:
            raise ValueError("Invalid dimensions or gene counts")
        # Each matrix column represents a cell. Count distinct cells for each
        # selected gene via set to avoid accidental duplicate sparse coordinates.
        seen_selected = {symbol: set() for symbol in TARGETS}
        row_to_gene = {}
        for gene, rows in selected_rows.items():
            for idx in rows:
                row_to_gene[idx] = gene
        for line in stream:
            if not line.strip() or line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) != 3:
                raise ValueError("Malformed sparse line")
            row, column = int(parts[0]), int(parts[1])
            count = float(parts[2])
            if not (1 <= row <= n_genes and 1 <= column <= cells and 0 <= count < float("inf")):
                raise ValueError("Invalid coordinate or UMI")
            parsed_records += 1
            counts_total += count
            if row in row_to_gene:
                gene = row_to_gene[row]
                totals[gene] += count
                if count > 0:
                    seen_selected[gene].add(column)
            if parsed_records > expected_entries:
                raise ValueError("Sparse matrix has too many entries")
    if parsed_records != expected_entries:
        raise ValueError("MatrixMarket EOF mismatch")
    for key, cells_set in seen_selected.items():
        detected_cells[key] = len(cells_set)
    return {
        "genes": n_genes,
        "cell_barcodes": cells,
        "nonzero_coordinate_records": parsed_records,
        "total_umi_count": counts_total,
        "total_umi_per_cell": counts_total / cells,
        "selected_gene_umis": dict(totals),
        "selected_genes_detected_cells": dict(detected_cells),
        "cell_barcodes_are_not_independent_donors": True,
    }


def inspect_one(index, group):
    url = geo_soft_url("GSE331114")
    raw, _ = download_metadata(url)
    info, samples = parse_soft(raw)
    if info.get("accession") != "GSE331114":
        raise ValueError("Wrong GEO series accession")
    cohort = sorted(
        [r for r in samples.values() if derive_group(r) == group],
        key=lambda row: row["accession"]
    )
    if len(cohort) != 4 or not 0 <= index < 4:
        raise ValueError("GSE331114 original case-control cohort size mismatch")
    sample = cohort[index]
    files = group_files(sample)
    record = {
        "accession": sample["accession"],
        "group": group, "group_position_0based": index,
        "source": "NCBI_GEO_GSE331114",
        "status": "HOLD",
        "matched_genotype": False,
        "clinical_gwas_patient_matching": False,
        "sample_description": sample["title"],
        "donor_unique_identity": "NOT_INDEPENDENTLY_VERIFIED",
    }
    if set(files) != {"matrix", "features", "barcodes"}:
        record["status"] = "MEX_FILES_MISSING"
        record["available_types"] = sorted(files)
        return record
    with tempfile.TemporaryDirectory(prefix="gse331114_") as temp:
        root = Path(temp)
        paths = {}
        for part in ("features", "barcodes", "matrix"):
            target = root / (part + ".gz")
            try:
                evidence = download(files[part], target, CAPS[part])
                record[part + "_sha256_not_source_verified"] = evidence["sha256_local"]
                record[part + "_compressed_bytes"] = evidence["downloaded_bytes"]
                paths[part] = target
            except (OSError, ValueError, TimeoutError) as exc:
                record["status"] = "TRANSFER_HOLD"
                record["error"] = str(exc)[:160]
                return record
        try:
            nfeatures, symbols, duplicates = parse_features(paths["features"])
            if duplicates:
                record["status"] = "DUPLICATE_TARGET_SYMBOL_REVIEW"
                record["duplicated_symbols"] = duplicates
                return record
            matrix = stream_gene_counts(paths["matrix"], nfeatures, symbols)
            with iter_text(paths["barcodes"]) as inp:
                barcodes = sum(1 for x in inp if x.strip())
            if barcodes != matrix["cell_barcodes"]:
                record["status"] = "BARCODE_DIMENSION_MISMATCH"
                return record
            record.update({key: value for key, value in matrix.items()
                           if key not in ("selected_gene_umis", "selected_genes_detected_cells")})
            record["status"] = "SAMPLE_MEX_GENE_QC_PASS"
            record["target_gene_rows"] = []
            for symbol in TARGETS:
                count = matrix["selected_gene_umis"].get(symbol, 0.0)
                cells = matrix["selected_genes_detected_cells"].get(symbol, 0)
                record["target_gene_rows"].append({
                    "gene_symbol": symbol, "gene_present_in_features": symbol in symbols,
                    "umis_sum": int(count) if count.is_integer() else count,
                    "umis_per_cell": round(count / barcodes, 6),
                    "fraction_cells_nonzero": round(cells / barcodes, 6),
                    "sample_group": group,
                    "sample_GSM": sample["accession"],
                })
        except (OSError, ValueError, EOFError) as exc:
            record["status"] = "MATRIX_QC_HOLD"
            record["error"] = str(exc)[:160]
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--group", choices=["AIS", "Control"], required=True)
    parser.add_argument("--index", type=int, choices=list(range(4)), required=True)
    parser.add_argument("--out-dir", default="artifacts/is_stroke_single_sample")
    args = parser.parse_args()
    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    try:
        record = inspect_one(args.index, args.group)
    except (OSError, ValueError, TimeoutError, EOFError) as exc:
        record = {"group": args.group, "index": args.index, "status": "SOURCE_UNAVAILABLE",
                  "error": type(exc).__name__ + ": " + str(exc)[:180]}
    (output / "gene_counts.json").write_text(json.dumps(
        record, indent=2, ensure_ascii=False
    ) + "\n", encoding="utf-8")
    if record.get("target_gene_rows"):
        columns = ["sample_GSM", "sample_group", "gene_symbol", "gene_present_in_features",
                   "umis_sum", "umis_per_cell", "fraction_cells_nonzero"]
        with open(output / "gene_counts.tsv", "w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            writer.writerows(record["target_gene_rows"])
    summary = (
        f"## GSE331114 {args.group} index {args.index}\n\n"
        f"**{record['status']}** | GSM {record.get('accession','unknown')} | "
        f"barcodes {record.get('cell_barcodes','-')} | genes {record.get('genes','-')} | "
        f"total UMIs {record.get('total_umi_count','-')}\n\n"
        "This is one **public GEO sample**; no DE/association test, no clinical individual matching."
    )
    (output / "summary.md").write_text(summary + "\n", encoding="utf-8")
    print(summary)
    if os.getenv("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a", encoding="utf-8") as stream:
            stream.write(summary + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
