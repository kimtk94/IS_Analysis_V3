#!/usr/bin/env python3
"""IS P3A: audit ROSMAP QTL sample N across R6/R7/R8 and source-input inventory.

All sources are official small public metadata. Never infer the availability of
qtlmap's *private* exact sample_meta, normalized matrix, or covariates from the
published output manifest. Keep scientific coloc 36/36 HOLD.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import io
import json
import re
from pathlib import Path

SMALL = importlib.import_module(
    "scripts.is.audit_p3a_36_original_zenodo_metadata_gwas_beta_model")
MATRIX = importlib.import_module(
    "scripts.is.audit_p3a_rosmap_2021_count_matrix_intersection")
SOURCE_BASE = ("https://raw.githubusercontent.com/eQTL-Catalogue/"
               "eQTL-Catalogue-resources/master/")
QTLMAP_DOC = ("https://raw.githubusercontent.com/eQTL-Catalogue/"
              "qtlmap/master/docs/usage.md")
RELEASE_PATHS = {
    "R6": "data_tables/dataset_metadata_r6.tsv",
    "R7": "data_tables/dataset_metadata_r7.tsv",
    "R8": "data_tables/dataset_metadata_r8.tsv",
}
FTP_MANIFEST = "tabix/tabix_ftp_paths.tsv"
TARGET_STUDY = "QTS000025"
TARGET_DATASET = "QTD000434"
TARGET_TISSUE = "brain_naive"
TARGET_N = 560
QUANTIFICATION = {"ge", "exon", "tx", "txrev", "leafcutter"}
OUTPUT_FIELDS = {"ftp_path", "ftp_cs_path", "ftp_lbf_path"}
INPUT_FIELDS = {"count_matrix", "pheno_meta", "sample_meta", "vcf", "covariates"}


def require(value, reason):
    if not value:
        raise ValueError(reason)


def parse_tsv(raw, required_columns):
    text = raw.decode("utf-8-sig")
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    require(reader.fieldnames and len(reader.fieldnames) == len(set(reader.fieldnames))
            and set(required_columns).issubset(reader.fieldnames),
            "Official source metadata schema missing required fields")
    rows = list(reader)
    require(bool(rows) and all(None not in row and None not in row.values()
                               for row in rows), "Malformed source metadata rows")
    return set(reader.fieldnames), rows


def release_audit(version, raw):
    require(version in RELEASE_PATHS, "Unsupported release version")
    fields, rows = parse_tsv(raw, ("study_id", "dataset_id", "study_label",
                                   "sample_group", "sample_size", "quant_method"))
    study = [x for x in rows if x["study_id"] == TARGET_STUDY and
             x["study_label"] == "ROSMAP" and
             x["sample_group"] == TARGET_TISSUE]
    require(study and len({x["dataset_id"] for x in study}) == len(study),
            "ROSMAP per-release study row identities missing or duplicate")
    require(all(int(x["sample_size"]) == TARGET_N for x in study),
            "ROSMAP study sample N differs from 560 within official release")
    genes = [x for x in study if x["dataset_id"] == TARGET_DATASET]
    require(len(genes) == 1 and genes[0]["quant_method"] == "ge",
            "QTD000434 not one and the same ROSMAP gene expression dataset")
    found_methods = {x["quant_method"] for x in study}
    require(QUANTIFICATION.issubset(found_methods), "R6-R8 core ROSMAP quantification missing")
    require(len(study) <= 15, "Unbounded ROSMAP release sample group")
    return {
        "source_release": version, "study_id": TARGET_STUDY,
        "gene_expression_dataset_id": TARGET_DATASET,
        "gene_expression_nominal_tissue_N": TARGET_N,
        "qtl_group": TARGET_TISSUE,
        "ROSMAP_study_quantification_count": len(study),
        "ROSMAP_study_quantification_methods": sorted(found_methods),
        "official_metadata_declares_individual_sample_IDs": False,
        "qtlmap_input_sample_intersection_measured": False,
        "direct_quantitative_trait_sdY_measured": False,
    }


def ftp_audit(raw):
    fields, rows = parse_tsv(raw, ("study_id", "dataset_id", "study_label",
                                   "sample_group", "sample_size", "quant_method",
                                   *sorted(OUTPUT_FIELDS)))
    ros = [x for x in rows if x["dataset_id"] == TARGET_DATASET]
    require(len(ros) == 1 and ros[0]["study_id"] == TARGET_STUDY and
            ros[0]["study_label"] == "ROSMAP" and
            ros[0]["sample_group"] == TARGET_TISSUE and
            int(ros[0]["sample_size"]) == TARGET_N and
            ros[0]["quant_method"] == "ge",
            "Wrong original ROSMAP gene expression FTP dataset")
    row = ros[0]
    for col in OUTPUT_FIELDS:
        path = row[col]
        require(path.startswith(
            "ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/"),
            "Official FTP source path is not trusted publisher")
        require("/QTS000025/QTD000434/" in path,
                "Official FTP source path does not resolve exact ROSMAP dataset")
    require(all(col not in fields for col in INPUT_FIELDS),
            "Unexpected input type columns exist: review inventory classification")
    return {
        "manifest_dataset_id": TARGET_DATASET,
        "manifest_study_id": TARGET_STUDY,
        "manifest_sample_N": TARGET_N,
        "manifest_column_names": sorted(fields),
        "official_publisher_manifest_output_columns": sorted(OUTPUT_FIELDS),
        "has_direct_qtlmap_count_matrix_column": "count_matrix" in fields,
        "has_direct_qtlmap_sample_meta_column": "sample_meta" in fields,
        "has_direct_qtlmap_genotype_vcf_column": "vcf" in fields,
        "has_direct_qtlmap_covariate_column": "covariates" in fields,
        "official_manifest_contains_full_cis_sumstats_path": True,
        "official_manifest_contains_susie_credible_set_path": True,
        "official_manifest_contains_susie_lbf_path": True,
        "output_manifest_not_inventory_of_all_possible_qtlmap_inputs": True,
        "Release7_normalized_expression_public_availability_determined": False,
        "Release7_sample_meta_public_availability_determined": False,
        "file_absence_outside_manifest_claimed": False,
    }


def qtlmap_model_audit(raw):
    prior = SMALL.PRIOR.parse_qtlmap(raw)
    require(prior["source_status"] ==
            "PASS_OFFICIAL_QTLMAP_NORMALIZED_MATRIX_AND_SAMPLE_INTERSECTION",
            "Official qtlmap input model no longer understood")
    return {
        "input_requirements_proven_by_official_workflow": {
            "normalized_count_matrix": True, "phenotype_metadata": True,
            "sample_meta": True, "genotype_vcf": True,
            "intersection_by_sample_ids": True,
        },
        "original_QTD000434_studyFile_recovered": False,
        "original_QTD000434_normalized_matrix_recovered": False,
        "original_QTD000434_sample_meta_recovered": False,
        "original_QTD000434_genotype_model_recovered": False,
        "original_QTD000434_covariate_matrix_recovered": False,
    }


def assemble(manifest, releases, ftp, qtl, provenance):
    groups = manifest.get("groups", [])
    ros = MATRIX.manifest_groups(manifest)
    require(len(groups) == 36 and len(ros) == 15 and
            set(releases) == set(RELEASE_PATHS) and
            all(releases[version]["gene_expression_nominal_tissue_N"] == 560
                for version in releases) and
            ftp["manifest_dataset_id"] == TARGET_DATASET and
            qtl["input_requirements_proven_by_official_workflow"]["sample_meta"],
            "Required 36-group original metadata crosswalk missing")
    require(all(releases[k]["study_id"] == TARGET_STUDY for k in releases),
            "ROSMAP official study ID differs by release")
    return {
        "schema": "is-p3a-36-ROSMAP-release6-7-8-input-public-inventory-v1",
        "technical_status": "PASS_OFFICIAL_R6_R7_R8_N_AND_OUTPUT_MANIFEST_AUDIT",
        "scientific_status": "HOLD_ORIGINAL_RELEASE7_QTLMAP_INPUT_AND_GWAS_BETA_UNIT",
        "official_publisher_sources": provenance,
        "release_comparison": releases,
        "r6_r7_r8_all_declare_560_for_QTD000434": True,
        "560_sample_N_already_in_R6_not_newly_R7": True,
        "R6_R7_R8_sample_overlap_independently_verified": False,
        "2021_Zenodo_QC_pass_569_to_560_exclusion_reason_proven": False,
        "ROSMAP_2021_QC_pass_count_reference_only": 569,
        "historical_2021_QC_count_source_reverified_in_this_run": False,
        "public_file_path_inventory": ftp,
        "official_qtlmap_input_model": qtl,
        "n_original_gene_tissue_groups": 36,
        "n_original_ROSMAP_groups": 15,
        "n_original_distinct_genes": 20,
        "n_gene_specific_original_Release7_sdY_measured": 0,
        "n_gene_specific_original_Release7_effective_N_verified": 0,
        "n_original_GIGASTROKE_beta_scale_verified": 0,
        "n_coloc_abf_eligible": 0, "n_coloc_abf_executed": 0,
        "n_coloc_susie_eligible": 0, "n_coloc_susie_executed": 0,
        "donor_or_sample_identifiers_exported": False,
        "groups": [
            {"gene_id": row["gene_id"], "gene_symbol": row["gene_symbol"],
             "QTL_dataset": row["QTL_dataset"],
             "scientific_status": "HOLD_R7_QTLMAP_N_SDY_GWAS_BETA",
             "coloc_eligible": False}
            for row in groups
        ],
    }


def fetch_official(url, cap=400000):
    raw, meta = SMALL.get(url, {"raw.githubusercontent.com"}, cap)
    return raw, {"source_url": url, "SHA256": meta["SHA256"],
                 "original_bytes": meta["read_bytes"]}


def execute(manifest):
    results, provenance = {}, {}
    for version, path in RELEASE_PATHS.items():
        print("SOURCE_STAGE_OFFICIAL_" + version, flush=True)
        raw, info = fetch_official(SOURCE_BASE + path)
        results[version] = release_audit(version, raw)
        provenance[version] = info
    print("SOURCE_STAGE_OFFICIAL_FTP_MANIFEST", flush=True)
    raw, info = fetch_official(SOURCE_BASE + FTP_MANIFEST)
    ftp = ftp_audit(raw)
    provenance["FTP_manifest"] = info
    print("SOURCE_STAGE_OFFICIAL_QTLMAP_INPUT_SCHEMA", flush=True)
    raw, info = fetch_official(QTLMAP_DOC)
    qtl = qtlmap_model_audit(raw)
    provenance["qtlmap_specification"] = info
    print("SOURCE_STAGE_36_GROUP_HOLD", flush=True)
    return assemble(manifest, results, ftp, qtl, provenance)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    out = {
        "schema": "is-p3a-36-ROSMAP-release6-7-8-input-public-inventory-v1",
        "technical_status": "HOLD_OFFICIAL_RELEASE_METADATA",
        "scientific_status": "HOLD",
        "n_coloc_abf_eligible": 0, "n_coloc_abf_executed": 0,
    }
    try:
        manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
        out = execute(manifest)
    except (OSError, TypeError, ValueError, KeyError, UnicodeError, TimeoutError) as err:
        out["error_class_only"] = type(err).__name__
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, ensure_ascii=False, indent=2) + "\n",
                        encoding="utf-8")
    print(json.dumps({
        "technical_status": out["technical_status"],
        "scientific_status": out["scientific_status"],
        "release_QTD000434_N": {k:v["gene_expression_nominal_tissue_N"]
            for k,v in out.get("release_comparison", {}).items()},
        "n_coloc_abf_executed": 0,
        "error_class_only": out.get("error_class_only"),
    }, ensure_ascii=False))
    return 0 if out["technical_status"].startswith("PASS_") else 1


if __name__ == "__main__":
    raise SystemExit(main())
