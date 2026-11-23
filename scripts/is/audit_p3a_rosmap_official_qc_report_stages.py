#!/usr/bin/env python3
"""Audit independent official ROSMAP QC-report counts against 2021 metadata and R7.

The 2021 metadata, a historical QC report, and release-7 summaries are
DIFFERENT analysis-stage publications. Numerical differences do not identify
individuals or an exclusion reason; never promote into QTL N or sdY.
"""
from __future__ import annotations

import argparse
import hashlib
from html.parser import HTMLParser
import importlib
import json
import re
from pathlib import Path

SMALL = importlib.import_module("scripts.is.audit_p3a_36_original_zenodo_metadata_gwas_beta_model")
MATRIX = importlib.import_module("scripts.is.audit_p3a_rosmap_2021_count_matrix_intersection")

QC_URL = "https://www.ebi.ac.uk/eqtl/static/QC_reports/ROSMAP_QC_report.html"
MAX_QC_BYTES = 4_500_000
REPORT_FIELDS = {
    "qc_before": ("Number of samples in study before quality control", 581),
    "qc_after": ("Number of samples in study after quality control", 576),
    "genotype_sample_count": ("Number of genotyped individuals in the current study ROSMAP", 581),
    "matched_genotypes": ("Number of matching genotypes", 576),
    "nonmatching_genotypes": ("Number of non-matching genotypes", 0),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


class TextExtractor(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.parts = []
        self.suppressed = 0

    def handle_starttag(self, tag, attrs):
        if tag in {"script", "style"}:
            self.suppressed += 1

    def handle_endtag(self, tag):
        if tag in {"script", "style"} and self.suppressed:
            self.suppressed -= 1

    def handle_data(self, data):
        if not self.suppressed:
            self.parts.append(data)


def parse_original_qc_report(raw: bytes) -> dict:
    require(isinstance(raw, bytes) and len(raw) > 30, "Missing official ROSMAP QC source")
    parser = TextExtractor()
    parser.feed(raw.decode("utf-8-sig"))
    text = " ".join(" ".join(parser.parts).split())
    require("ROSMAP QC report" in text and
            "Concordance between expression and genotype data" in text,
            "Publisher report identity or schema unrecognized")
    counts = {}
    for key, (field, expected) in REPORT_FIELDS.items():
        pattern = r"\b" + re.escape(field) + r"\s*:\s*(\d+)"
        hits = re.findall(pattern, text, flags=re.I)
        require(len(hits) == 1, "Publisher QC report field missing or ambiguous: " + key)
        value = int(hits[0])
        require(value == expected, "Publisher QC report sample count changed: " + key)
        counts[key] = value
    require(counts["qc_after"] <= counts["qc_before"] and
            counts["matched_genotypes"] + counts["nonmatching_genotypes"] ==
            counts["qc_after"], "Publisher QC counts inconsistent")
    return {
        "source_status": "PASS_PUBLISHER_HISTORIC_QC_REPORT_COUNTS",
        "publisher": "EMBL-EBI eQTL Catalogue",
        "report_title": "ROSMAP QC report",
        "source_url": QC_URL,
        "counts_only": counts,
        "report_stage_source_release_certified_same_as_Release7": False,
        "report_stage_source_release_certified_same_as_2021_Zenodo_metadata": False,
        "source_individual_identifiers_exported": False,
    }


def get_original_2021_counts():
    raw, record_meta = SMALL.get(SMALL.RECORD, SMALL.ALLOWED_RECORD_HOSTS,
                                 SMALL.MAX_RECORDBYTES)
    entries, provenance = SMALL.require_original_record(raw)
    meta, file_meta = SMALL.get_official_file("ROSMAP.tsv", entries["ROSMAP.tsv"])
    _, passed = MATRIX.sample_sets(meta)
    require(len(passed) == 569, "Original ROSMAP 2021 both QC PASS drift")
    return {"2021_original_Zenodo_RNA_genotype_QC_pass_N": len(passed),
            "metadata_publisher_MD5": SMALL.FILES["ROSMAP.tsv"]["md5"],
            "metadata_publisher_bytes": SMALL.FILES["ROSMAP.tsv"]["bytes"],
            "metadata_actual_SHA256": file_meta["SHA256"],
            "Zenodo_record_response_SHA256_unpinned": record_meta["SHA256"],
            "Zenodo_manifest_verified": provenance["official_record_n_public_files"] == 32}


def assemble(groups_manifest: dict, report: dict, original2021: dict,
             release7: dict, report_sha: str, release7_sha: str) -> dict:
    groups = groups_manifest.get("groups", [])
    selected = MATRIX.manifest_groups(groups_manifest)
    require(len(selected) == 15 and report["source_status"] ==
            "PASS_PUBLISHER_HISTORIC_QC_REPORT_COUNTS" and
            original2021["2021_original_Zenodo_RNA_genotype_QC_pass_N"] == 569 and
            release7["Release7_declared_tissue_N"] == 560,
            "Original study/QTL manifest or verified stages mismatched")
    return {
        "schema": "is-p3a-36-ROSMAP-original-official-QC-report-multistage-v1",
        "technical_status": "PASS_OFFICIAL_ROSMAP_INDEPENDENT_QC_REPORT_STAGE_COUNTS",
        "scientific_status": "HOLD_RELEASE7_EXACT_SAMPLE_INTERSECTION_SDY_AND_GWAS_EFFECT_UNIT",
        "n_original_QTL_gene_tissue_assays": 36,
        "n_original_distinct_genes": 20,
        "n_ROSMAP_QTD000434_groups": 15,
        "independent_official_historical_report": report,
        "observed_2021_Zenodo_original_metadata_aggregate": original2021,
        "release7_public_reference": release7,
        "source_files_SHA256": {
            "official_EBI_QC_report_raw_SHA256": report_sha,
            "Release7_dataset_metadata_raw_SHA256": release7_sha,
        },
        "cross_source_descriptive_only": {
            "official_report_after_QC_minus_2021_metadata_both_QC": 576 - 569,
            "2021_metadata_both_QC_minus_Release7_nominal_N": 569 - 560,
            "official_report_after_QC_minus_Release7_nominal_N": 576 - 560,
            "these_differences_are_not_verified_sequential_sample_exclusions": True,
            "historical_QC_report_and_2021_metadata_same_analytical_release": False,
            "2021_metadata_and_Release7_exact_sample_identity_verified": False,
            "Release7_qtlmap_normalized_expression_observed": False,
            "Release7_qtlmap_sample_meta_genotype_covariates_intersection_observed": False,
            "ROSMAP_15_gene_expression_sdY_measured": False,
            "ROSMAP_15_gene_regression_N_verified": False,
            "GIGASTROKE_exact_beta_log_OR_verified": False,
        },
        "n_coloc_abf_eligible": 0, "n_coloc_abf_executed": 0,
        "n_coloc_susie_eligible": 0, "n_coloc_susie_executed": 0,
        "original_sample_IDs_in_results": False,
        "groups": [
            {"QTL_dataset": g["QTL_dataset"], "gene_id": g["gene_id"],
             "gene_symbol": g["gene_symbol"],
             "status": "HOLD_EXACT_RELEASE7_N_SDY_AND_GWAS_UNITS",
             "scientifically_eligible_for_coloc": False}
            for g in groups],
    }


def execute(manifest: dict):
    print("SOURCE_STAGE_1_EBI_OFFICIAL_REPORT", flush=True)
    raw, metadata = SMALL.get(QC_URL, {"www.ebi.ac.uk"}, MAX_QC_BYTES)
    report = parse_original_qc_report(raw)
    print("SOURCE_STAGE_2_2021_ORIGINAL_METADATA_MD5", flush=True)
    original = get_original_2021_counts()
    print("SOURCE_STAGE_3_RELEASE7_REFERENCE", flush=True)
    r7, r7metadata = SMALL.get(MATRIX.R7_URL, {"raw.githubusercontent.com"}, 500000)
    release7 = MATRIX.release7(r7)
    print("SOURCE_STAGE_4_36_GROUP_SCIENTIFIC_HOLD", flush=True)
    return assemble(manifest, report, original, release7,
                    metadata["SHA256"], r7metadata["SHA256"])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    result = {
        "schema": "is-p3a-36-ROSMAP-original-official-QC-report-multistage-v1",
        "technical_status": "HOLD_OFFICIAL_PUBLISHER_SOURCE",
        "scientific_status": "HOLD",
        "n_coloc_abf_eligible": 0,
        "n_coloc_abf_executed": 0,
    }
    try:
        manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
        result = execute(manifest)
    except (OSError, ValueError, TypeError, KeyError, UnicodeError, TimeoutError) as err:
        # Never log exception text: a source error could include sample metadata.
        result["error_class_only"] = type(err).__name__
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n",
                        encoding="utf-8")
    print(json.dumps({
        "technical_status": result["technical_status"],
        "scientific_status": result["scientific_status"],
        "report_counts": result.get("independent_official_historical_report", {}).get("counts_only"),
        "n_original_QTL_gene_tissue_assays": result.get("n_original_QTL_gene_tissue_assays"),
        "n_coloc_abf_executed": 0,
        "error_class_only": result.get("error_class_only"),
    }, ensure_ascii=False))
    return 0 if result["technical_status"].startswith("PASS_") else 1


if __name__ == "__main__":
    raise SystemExit(main())
