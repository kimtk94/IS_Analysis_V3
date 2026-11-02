#!/usr/bin/env python3
"""P4 pilot: streaming raw-count pseudobulk for paired DEX/vehicle BeamA in pool C.

GSE313736 public in-vitro primary human TM cells, NOT clinical steroid IOP.
Raw matrices and donor IDs are never written to artifacts or the repository.
"""
import argparse
import csv
import gzip
import json
import statistics
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from ocular_p1_integrity import BASE, download, file_url, bar, lines

SAMPLES = {
    "vehicle": ("GSM9373278", "pool_C", "no-dex", "HTMC_pool-C_GRCh38"),
    "dex": ("GSM9373287", "pool_dexC", "dex", "HTMC_pool-dexC_GRCh38"),
}
ANNOTATION = "celltypist_chen2025_fine"
SCORE = "celltypist_chen2025_fine_score"
GENES = 36601


def read_metadata(path, min_confidence):
    by_arm = defaultdict(dict)
    seen = Counter()
    with gzip.open(path, "rt", encoding="utf8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        need = {"barcode", "library_id", "donor_id_demuxlet",
                "DEX Treatment", "droplet_type", ANNOTATION, SCORE}
        if not need.issubset(reader.fieldnames or []):
            raise ValueError("MISSING_METADATA_COLUMNS " +
                             str(sorted(need - set(reader.fieldnames or []))))
        for row in reader:
            cond = row["DEX Treatment"].strip().lower()
            lib = row["library_id"].strip()
            for arm, (_, _, expected_cond, expected_lib) in SAMPLES.items():
                if cond != expected_cond or lib != expected_lib:
                    continue
                seen[arm] += 1
                if row["droplet_type"] != "SNG" or row[ANNOTATION] != "BeamA":
                    continue
                try:
                    confidence = float(row[SCORE])
                except (TypeError, ValueError):
                    continue
                if confidence < min_confidence:
                    continue
                donor = row["donor_id_demuxlet"].strip()
                bc = bar(row["barcode"])
                if not donor or donor.lower() in ("na", "doublet", "unassigned"):
                    continue
                if bc in by_arm[arm]:
                    raise ValueError("DUPLICATE_METADATA_BARCODE")
                by_arm[arm][bc] = donor
    if not all(seen[arm] for arm in SAMPLES):
        raise ValueError("EXPECTED_METADATA_LIBRARY_NOT_FOUND " + str(dict(seen)))
    return by_arm, dict(seen)


def eligible_donors(by_arm, minimum_cells):
    a = Counter(by_arm["vehicle"].values())
    b = Counter(by_arm["dex"].values())
    selected = sorted(d for d in a.keys() & b.keys()
                      if a[d] >= minimum_cells and b[d] >= minimum_cells)
    return selected, {"vehicle": a, "dex": b}


def matrix_stream(path, sample_cols, group_for_column, gene_count,
                  n_groups, bytes_per_chunk=4 * 1024 * 1024):
    """Process entire coordinate file with bounded RAM and verify gzip CRC/NNZ."""
    pooled = np.zeros((gene_count, n_groups), dtype=np.int64)
    total_lines = selected_entries = 0
    selected_umis = 0
    with gzip.open(path, "rb") as stream:
        header = stream.readline().strip()
        if header != b"%%MatrixMarket matrix coordinate integer general":
            raise ValueError("UNSUPPORTED_MATRIX_FORMAT")
        while True:
            line = stream.readline()
            if not line:
                raise ValueError("MISSING_MATRIX_HEADER")
            if not line.startswith(b"%"):
                break
        rows, cols, nnz = map(int, line.split())
        if rows != gene_count or cols != sample_cols or nnz < 1:
            raise ValueError(f"MATRIX_DIMS_BAD {(rows, cols, nnz)}")
        carry = b""
        while True:
            chunk = stream.read(bytes_per_chunk)
            if not chunk:
                if carry:
                    chunk = carry + b"\n"
                    carry = b""
                else:
                    break
            else:
                chunk = carry + chunk
            end = chunk.rfind(b"\n")
            if end < 0:
                if len(chunk) > 1024:
                    raise ValueError("UNEXPECTED_LONG_MATRIX_LINE")
                carry = chunk
                continue
            ready, carry = chunk[:end + 1], chunk[end + 1:]
            numbers = np.fromstring(ready.decode("ascii"), sep=" ", dtype=np.int64)
            if len(numbers) == 0:
                continue
            if len(numbers) % 3:
                raise ValueError("INVALID_SPARSE_TRIPLES")
            triples = numbers.reshape(-1, 3)
            row_idx = triples[:, 0] - 1
            col_idx = triples[:, 1] - 1
            val = triples[:, 2]
            if (np.any(row_idx < 0) or np.any(row_idx >= rows)
                or np.any(col_idx < 0) or np.any(col_idx >= cols)
                or np.any(val <= 0)):
                raise ValueError("SPARSE_ENTRY_OUT_OF_RANGE_OR_NONPOSITIVE")
            total_lines += triples.shape[0]
            groups = group_for_column[col_idx]
            use = groups >= 0
            if np.any(use):
                flat_idx = row_idx[use] * n_groups + groups[use]
                sums = np.bincount(flat_idx, weights=val[use].astype(np.float64),
                                   minlength=gene_count * n_groups)
                # Counts remain exact as long as chunk-level group sums < 2^53.
                if sums.max() >= 2**53:
                    raise ValueError("FLOAT_ACCUMULATOR_PRECISION_RISK")
                pooled += np.rint(sums).astype(np.int64).reshape(
                    gene_count, n_groups)
                selected_entries += int(np.count_nonzero(use))
                selected_umis += int(val[use].sum())
        if total_lines != nnz:
            raise ValueError(f"NNZ_MISMATCH expected={nnz} actual={total_lines}")
    if int(pooled.sum()) != selected_umis:
        raise ValueError("UMI_CONSERVATION_FAILED")
    return pooled, {"nnz": int(total_lines),
                    "selected_nnz": int(selected_entries),
                    "selected_umis": int(selected_umis),
                    "gzip_crc": "PASS"}


def summary_numbers(vals):
    v = sorted(int(x) for x in vals)
    return {"min": min(v), "median": statistics.median(v), "max": max(v)}


def pilot(outdir, min_confidence, min_cells):
    report = {
        "dataset": "GSE313736",
        "stage": "P4_pool_C_paired_raw_count_pseudobulk",
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "status": "HOLD",
        "cell_annotation": ANNOTATION,
        "confidence_cutoff": min_confidence,
        "minimum_cells_per_donor_per_arm": min_cells,
        "study_scope": "In-vitro human TM cellular DEX response only; no clinical IOP phenotype",
        "arms": {},
        "no_donor_ids_saved": True,
        "raw_counts_not_uploaded": True,
    }
    with tempfile.TemporaryDirectory(prefix="ocular-p4-") as tmp:
        root = Path(tmp)
        meta = root / "metadata.gz"
        report["metadata_transfer"] = download(
            BASE + "series/GSE313nnn/GSE313736/suppl/"
            "GSE313736_HTMC_scRNA_cell_metadata.tsv.gz",
            meta, 12 * 1024 * 1024)
        maps, all_cell_counts = read_metadata(meta, min_confidence)
        donors, eligible_counts = eligible_donors(maps, min_cells)
        report["metadata_cell_counts_in_selected_libraries"] = all_cell_counts
        report["eligible_paired_donors"] = len(donors)
        if len(donors) < 3:
            raise ValueError("INSUFFICIENT_PAIRED_DONORS_FOR_PILOT")
        report["eligible_cells_per_arm"] = {
            arm: int(sum(eligible_counts[arm][d] for d in donors))
            for arm in SAMPLES
        }
        sums_by_arm = {}
        for arm, (gsm, pool, cond, lib) in SAMPLES.items():
            bar_path = root / (arm + "_barcodes.gz")
            feat_path = root / (arm + "_features.gz")
            matrix_path = root / (arm + "_matrix.mtx.gz")
            bar_transfer = download(file_url(gsm, pool, "barcodes.tsv.gz"),
                                    bar_path, 12 * 1024 * 1024)
            feat_transfer = download(file_url(gsm, pool, "features.tsv.gz"),
                                     feat_path, 12 * 1024 * 1024)
            barcodes = lines(bar_path)
            features = lines(feat_path)
            if len(features) != GENES or len({x.split("\t")[0] for x in features}) != GENES:
                raise ValueError("GENE_FEATURE_COUNT_OR_DUPLICATE")
            report.setdefault("gene_order_digests", {})[arm] = __import__("hashlib").sha256(
                ("\n".join(x.split("\t")[0] for x in features) + "\n").encode()
            ).hexdigest()
            if len(set(map(bar, barcodes))) != len(barcodes):
                raise ValueError("DUPLICATE_10X_BARCODES")
            donor_indexes = {d: i for i, d in enumerate(donors)}
            lookup = maps[arm]
            ix = np.asarray([donor_indexes.get(lookup.get(bar(bc)), -1)
                             for bc in barcodes], dtype=np.int32)
            n_map = int(np.count_nonzero(ix >= 0))
            expected = report["eligible_cells_per_arm"][arm]
            # Metadata singlets must be present in the source matrix. No silent loss.
            report["arms"][arm] = {
                "gsm": gsm,
                "n_source_barcodes": len(barcodes),
                "selected_barcodes_in_matrix": n_map,
                "selected_cells_in_metadata": expected,
                "selected_barcode_recall": round(n_map / max(1, expected), 5),
                "barcode_sha256": bar_transfer["sha256"],
                "features_sha256": feat_transfer["sha256"],
            }
            if n_map < .95 * expected:
                raise ValueError(f"{arm}_SELECTED_BARCODE_RECALL_LT_95PCT")
            report["arms"][arm]["matrix_transfer"] = download(
                file_url(gsm, pool, "matrix.mtx.gz"), matrix_path,
                620 * 1024 * 1024)
            totals, integrity = matrix_stream(
                matrix_path, len(barcodes), ix, GENES, len(donors))
            report["arms"][arm]["matrix_integrity"] = integrity
            report["arms"][arm]["cell_depth_selected"] = n_map
            report["arms"][arm]["pseudobulk_umi_per_donor"] = summary_numbers(
                totals.sum(axis=0))
            report["arms"][arm]["pseudobulk_detected_genes_per_donor"] = summary_numbers(
                np.count_nonzero(totals, axis=0))
            if np.any(totals.sum(axis=0) <= 0):
                raise ValueError("EMPTY_PSEUDOBULK_DONOR")
            sums_by_arm[arm] = totals
            matrix_path.unlink(missing_ok=True)
        if report["gene_order_digests"]["vehicle"] != report["gene_order_digests"]["dex"]:
            raise ValueError("GENE_ORDER_DIFFERENT_BETWEEN_ARMS")
        report["paired_donor_count"] = len(donors)
        report["total_umi_in_selected_cells"] = int(
            sum(x.sum() for x in sums_by_arm.values()))
        report["genes_expressed_in_either_arm"] = int(
            np.count_nonzero((sums_by_arm["vehicle"] +
                              sums_by_arm["dex"]).sum(axis=1)))
        report["genes_detected_in_at_least_three_donors_each_arm"] = int(
            np.count_nonzero(
                (np.count_nonzero(sums_by_arm["vehicle"], axis=1) >= 3)
                & (np.count_nonzero(sums_by_arm["dex"], axis=1) >= 3)))
        report["status"] = "GO_POOL_C_RAW_PSEUDOBULK_QC"
        # No per-donor counts, genotype or identifiers are written to results.
    return report


def self_test():
    # 3 genes x 4 cells x 5 nonzeros, two biological groups.
    # Demonstrates true raw integer UMI summation, not mean expression/CPM.
    with tempfile.TemporaryDirectory() as td:
        fp = Path(td) / "fixture.mtx.gz"
        with gzip.open(fp, "wb") as w:
            w.write(b"%%MatrixMarket matrix coordinate integer general\n"
                    b"% fixture\n3 4 5\n1 1 2\n2 1 1\n1 2 3\n3 3 4\n2 4 5\n")
        selected = np.array([0, 0, 1, -1], dtype=np.int32)
        sums, check = matrix_stream(fp, 4, selected, 3, 2, bytes_per_chunk=16)
        expected = np.array([[5, 0], [1, 0], [0, 4]])
        if not np.array_equal(sums, expected) or check["selected_umis"] != 10:
            raise ValueError("SYNTHETIC_PSEUDOBULK_TEST_FAILED")
    return True


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", default="ocular-p4-pilot.json")
    p.add_argument("--confidence", type=float, default=0.75)
    p.add_argument("--min-cells", type=int, default=30)
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    if args.self_test:
        self_test()
        print("P4_SELF_TEST_PASS")
        return
    report = {"status": "HOLD", "run_utc": datetime.now(timezone.utc).isoformat()}
    try:
        report = pilot(Path(args.output).parent, args.confidence, args.min_cells)
    except Exception as exc:
        report["status"] = "HOLD_P4_PILOT"
        report["error"] = type(exc).__name__ + ": " + str(exc)[:300]
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print("P4_RESULT", report["status"],
          "paired_donors", report.get("paired_donor_count"),
          "error", report.get("error"))
    if report["status"] != "GO_POOL_C_RAW_PSEUDOBULK_QC":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
