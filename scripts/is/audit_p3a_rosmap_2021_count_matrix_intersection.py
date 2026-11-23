#!/usr/bin/env python3
"""Audit public ROSMAP 2021 count-matrix sample intersection, not Release7 sdY.

All pseudonymous IDs remain in process memory only. The 2021 processed count
matrix is not the exact Release7 normalized qtlmap phenotype input.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib
import io
import json
import re
import zlib
from pathlib import Path

P = importlib.import_module("scripts.is.audit_p3a_36_original_zenodo_metadata_gwas_beta_model")
COUNT_NAME = "ROSMAP_merged_counts_pseudo.tsv.gz"
COUNT_BYTES = 25484860
COUNT_MD5 = "d215dc801f687e15d5b938356d80d183"
R7_URL = ("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/"
          "master/data_tables/dataset_metadata_r7.tsv")
GENE = re.compile(r"^ENSG\d{11}(?:\.\d+)?$")


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def sample_sets(raw):
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), delimiter="\t")
    require(reader.fieldnames is not None and
            {"pseudo_id", "qtl_group", "rna_qc_passed", "genotype_qc_passed"}.issubset(reader.fieldnames),
            "Original ROSMAP 2021 metadata schema absent")
    all_ids, passed = set(), set()
    for row in reader:
        if row["qtl_group"] != "brain_naive":
            continue
        sid = row["pseudo_id"].strip()
        require(sid and sid not in all_ids, "Duplicate or empty original pseudonymous ID")
        r = P.source_qc_flag(row["rna_qc_passed"])
        g = P.source_qc_flag(row["genotype_qc_passed"])
        require(r is not None and g is not None, "Original QC flags not recognized")
        all_ids.add(sid)
        if r and g:
            passed.add(sid)
    require(len(all_ids) == 581 and len(passed) == 569, "Original 2021 ROSMAP QC counts drift")
    return all_ids, passed


def release7(raw):
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), delimiter="\t")
    require(reader.fieldnames is not None and
            {"dataset_id", "study_label", "sample_group", "sample_size"}.issubset(reader.fieldnames),
            "Release7 source schema changed")
    rows = [r for r in reader if r["dataset_id"] == "QTD000434"]
    require(len(rows) == 1 and rows[0]["study_label"] == "ROSMAP" and
            rows[0]["sample_group"] == "brain_naive" and
            int(rows[0]["sample_size"]) == 560, "Release7 official ROSMAP QTL identity/N mismatch")
    return {"dataset_id": "QTD000434", "sample_group": "brain_naive",
            "study": "ROSMAP", "Release7_declared_tissue_N": 560}


def manifest_groups(m):
    g = m.get("groups", [])
    require(m.get("schema") == "is-p3a-36-real-group-keys-from-legacy-393-v1"
            and len(g) == 36 and len({(x["QTL_dataset"], x["gene_id"]) for x in g}) == 36
            and len({x["gene_id"] for x in g}) == 20,
            "Original 36 groups / 20 genes not authentic")
    ros = [x for x in g if x["QTL_dataset"] == "QTD000434"]
    require(len(ros) == 15, "Original ROSMAP 15-gene set changed")
    return ros


def count_matrix_aggregates(compressed, all_ids, passed, ros):
    require(isinstance(compressed, bytes) and len(compressed) > 20, "Original counts unavailable")
    require(len(all_ids) == 581 and len(passed) == 569 and len(ros) == 15,
            "Original sample and gene source mismatch")
    targets = {r["gene_id"] for r in ros}
    with gzip.GzipFile(fileobj=io.BytesIO(compressed)) as gz:
        with io.TextIOWrapper(gz, encoding="utf-8-sig", newline="") as stream:
            reader = csv.reader(stream, delimiter="\t")
            header = next(reader, None)
            require(header is not None and 2 < len(header) <= 2000,
                    "Original count matrix header invalid")
            mids = [x for x in header if x in all_ids]
            require(len(mids) > 0 and len(mids) == len(set(mids)),
                    "Original matrix pseudonymous sample column overlap missing or duplicated")
            matrix_ids = set(mids)
            genes_found = set()
            rows = 0
            for fields in reader:
                if not fields:
                    continue
                rows += 1
                require(len(fields) == len(header), "Count matrix row shape changed")
                for label in fields[:min(4, len(fields))]:
                    if GENE.fullmatch(label):
                        key = label.split(".")[0]
                        if key in targets:
                            genes_found.add(key)
                            break
            require(rows > 0, "Count matrix has no phenotype rows")
    # All donor IDs were needed only for set operations above and are never output.
    intersect = matrix_ids & passed
    return {
        "2021_original_matrix_role": "PROCESSED_COUNTS_NOT_RELEASE7_NORMALIZED_QTL_PHENOTYPE",
        "2021_matrix_metadata_pseudo_ID_column_count": len(matrix_ids),
        "2021_matrix_and_RNA_genotype_QC_passed_intersection_count": len(intersect),
        "2021_matrix_and_QC_failed_metadata_column_count": len(matrix_ids - passed),
        "2021_QC_passed_metadata_samples_absent_from_matrix_count": len(passed - matrix_ids),
        "2021_count_matrix_feature_rows": rows,
        "original_15_ROSMAP_gene_ids_in_2021_count_matrix": len(genes_found),
        "per_ROSMAP_gene_2021_counts_presence_only": [
            {"gene_id": x["gene_id"], "gene_symbol": x["gene_symbol"],
             "present": x["gene_id"] in genes_found} for x in ros],
        "2021_intersection_count_equals_Release7_560_DESCRIPTIVE": len(intersect) == 560,
        "release7_qtlmap_sample_intersection_directly_proven": False,
        "any_9_missing_donors_attributed_to_Release7_exclusion": False,
        "gene_specific_original_normalized_sdY_measured": False,
        "gene_specific_Release7_regression_N_verified": False,
        "sample_identifiers_exported": False,
    }


def assemble(manifest, meta, count, r7, hashes):
    ros = manifest_groups(manifest)
    all_ids, passed = sample_sets(meta)
    official = release7(r7)
    matrix = count_matrix_aggregates(count, all_ids, passed, ros)
    return {
        "schema": "is-p3a-36-original-2021-ROSMAP-count-matrix-intersection-v1",
        "technical_status": "PASS_ORIGINAL_2021_ROSMAP_COUNT_MATRIX_AGGREGATION",
        "scientific_status": "HOLD_EXACT_RELEASE7_QTLMAP_SAMPLE_N_SDY_AND_GWAS_UNITS",
        "original_source_integrity": hashes,
        "2021_metadata_total_unique_samples": len(all_ids),
        "2021_metadata_RNA_and_genotype_QC_passed": len(passed),
        "Release7_reference": official,
        "2021_matrix_aggregates": matrix,
        "n_gene_tissue_groups": 36, "n_distinct_genes": 20, "n_ROSMAP_assays": 15,
        "n_groups_direct_qtlmap_sdY_verified": 0,
        "n_groups_original_GIGASTROKE_beta_log_OR_verified": 0,
        "n_formal_coloc_abf_eligible": 0, "n_formal_coloc_abf_executed": 0,
        "n_formal_coloc_susie_eligible": 0, "n_formal_coloc_susie_executed": 0,
        "groups": [
            {"gene_id": x["gene_id"], "gene_symbol": x["gene_symbol"],
             "QTL_dataset": x["QTL_dataset"], "status": "HOLD",
             "coloc_eligible": False, "normalized_sdY_proven": False}
            for x in manifest["groups"]],
        "any_source_sample_ids_exported": False,
    }


def run_original_sources(manifest):
    record_raw, _ = P.get(P.RECORD, P.ALLOWED_RECORD_HOSTS, P.MAX_RECORDBYTES)
    entries, provenance = P.require_original_record(record_raw)
    obj = entries.get(COUNT_NAME)
    require(obj is not None and obj.get("size") == COUNT_BYTES and
            obj.get("checksum") == "md5:" + COUNT_MD5,
            "Original Zenodo record counts MD5 or size changed")
    meta, _ = P.get_official_file("ROSMAP.tsv", entries["ROSMAP.tsv"])
    count_url = P.official_file_url(obj, COUNT_NAME)
    count, cm = P.get(count_url, P.ALLOWED_RECORD_HOSTS, 27000000)
    require(len(count) == COUNT_BYTES and
            hashlib.md5(count, usedforsecurity=False).hexdigest() == COUNT_MD5,
            "Downloaded original count matrix publisher MD5 mismatch")
    r7, r7meta = P.get(R7_URL, {"raw.githubusercontent.com"}, 500000)
    return assemble(manifest, meta, count, r7, {
        "2021_Zenodo_record_SHA256": provenance["official_record_original_SHA256"],
        "2021_ROSMAP_metadata_MD5": P.FILES["ROSMAP.tsv"]["md5"],
        "2021_ROSMAP_count_matrix_MD5": COUNT_MD5,
        "2021_ROSMAP_count_matrix_compressed_bytes": COUNT_BYTES,
        "2021_ROSMAP_count_matrix_actual_SHA256": cm["SHA256"],
        "Release7_dataset_metadata_SHA256": r7meta["SHA256"],
        "Release7_dataset_metadata_url": R7_URL,
        "downloads_verified_against_publisher_hashes": True,
    })


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = {"schema": "is-p3a-36-original-2021-ROSMAP-count-matrix-intersection-v1",
           "technical_status": "HOLD_PUBLIC_ORIGINAL_SOURCE",
           "scientific_status": "HOLD", "n_formal_coloc_abf_eligible": 0,
           "n_formal_coloc_abf_executed": 0}
    try:
        manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
        out = run_original_sources(manifest)
    except (OSError, ValueError, KeyError, TypeError, UnicodeError,
            EOFError, gzip.BadGzipFile, zlib.error, csv.Error) as error:
        # Exception messages could embed raw pseudonymous source data: never print them.
        out["error_class"] = type(error).__name__
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "technical_status": out["technical_status"],
        "scientific_status": out["scientific_status"],
        "2021_metadata_qc_passed": out.get("2021_metadata_RNA_and_genotype_QC_passed"),
        "2021_matrix_aggregates": out.get("2021_matrix_aggregates"),
        "n_coloc_abf_eligible": 0}, ensure_ascii=False))
    return 0 if out["technical_status"].startswith("PASS_") else 1


if __name__ == "__main__":
    raise SystemExit(main())
