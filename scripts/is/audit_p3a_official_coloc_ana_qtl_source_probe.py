#!/usr/bin/env python3
"""IS P3A: official eQTL Catalogue AN/MAF quantitative coloc path and live QTL source probe.

Official coloc demonstrates beta, se^2, MAF, AN/2, type=quant. This is an
example, NOT approval to substitute AN/2 for a validated gene-wise regression N,
nor to model an ischemic-stroke case-control outcome as type=quant.
Only a bounded *partial* gzip probe of official QTD000434 source is performed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import io
import json
import math
import re
import urllib.parse
import urllib.request
import zlib
from pathlib import Path

PRIOR = importlib.import_module("scripts.is.audit_p3a_rosmap_release_input_inventory")
SMALL = importlib.import_module("scripts.is.audit_p3a_36_original_zenodo_metadata_gwas_beta_model")
REPO = "https://raw.githubusercontent.com/eQTL-Catalogue/"
SOURCES = {
    "official_columns": (REPO+"eQTL-Catalogue-resources/master/tabix/Columns.md"),
    "official_coloc_tutorial": (REPO+"eQTL-Catalogue-resources/master/tutorials/tabix_use_case.md"),
    "official_coloc_implementation": (REPO+"colocalisation/master/bin/run_coloc_batch.R"),
    "original_public_paths": (REPO+"eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv"),
}
RELEASE_N = 560
DATASET_ID = "QTD000434"
STUDY_ID = "QTS000025"
MAX_GZIP_RANGE = 131072
REQUIRED_QTL_COLUMNS = {
    "molecular_trait_id", "gene_id", "variant", "chromosome", "position", "ref", "alt",
    "maf", "an", "beta", "se", "pvalue", "rsid", "type"
}


def require(ok, msg):
    if not ok:
        raise ValueError(msg)


def official_method_checks(cols, tutorial, implementation):
    # Exact method implemented and the official published column definitions.
    for symbol in ("**an**", "**maf**", "**beta**", "**se**", "**alt**"):
        require(symbol in cols, "Publisher column definitions changed: " + symbol)
    require("sample size is equal to an/2" in cols.lower(),
            "Official AN/2 interpretation missing from columns specification")
    require("also the effect allele" in cols.lower(),
            "Official ALT effect allele specification missing")
    for name, source in (("tutorial", tutorial), ("implementation", implementation)):
        compact = re.sub(r"\s+", "", source)
        require("N=(eqtl_sumstats$an)[1]/2" in compact,
                "Original "+name+" no longer explicitly derives QTL N as first AN/2")
        require("varbeta=eqtl_sumstats$se^2" in compact and
                "MAF=eqtl_sumstats$maf" in compact and
                "beta=eqtl_sumstats$beta" in compact and
                'type="quant"' in compact,
                "Original "+name+" not the documented quantitative coloc model")
    return {
        "source_proof": "OFFICIAL_EQTL_CATALOGUE_COLOC_EXAMPLE_BETA_SE_MAF_AN_DIV2_QUANT",
        "QTL_beta_SE_and_MAF_proven_official_fields": True,
        "QTL_effect_allele_is_GRCh38_ALT": True,
        "QTL_AN_defined_total_alleles": True,
        "AN_div2_common_autosomal_approximation_official_documentation": True,
        "original_tutorial_uses_first_SNP_AN_div2_as_dataset_N": True,
        "official_tutorial_estimated_sdY_without_normalized_phenotype": True,
        "AN_div2_validated_exact_regression_sample_N_for_QTD000434": False,
        "official_tutorial_GWAS_type_quant_not_validated_for_stroke_case_control": True,
        "binary_stroke_GWAS_requires_type_cc_and_verified_case_fraction_s": True,
        "formal_coloc_H4_proven": False,
    }


def ftp_url_from_manifest(raw):
    output = PRIOR.ftp_audit(raw)
    _, rows = PRIOR.parse_tsv(raw, ("study_id", "dataset_id", "ftp_path"))
    selected = [r for r in rows if r["dataset_id"] == DATASET_ID]
    require(len(selected) == 1 and selected[0]["study_id"] == STUDY_ID,
            "Original public ROSMAP source path mismatch")
    original = selected[0]["ftp_path"]
    parsed = urllib.parse.urlsplit(original)
    require(parsed.scheme == "ftp" and parsed.hostname == "ftp.ebi.ac.uk" and
            parsed.path.endswith("/QTD000434.all.tsv.gz") and
            "/QTS000025/QTD000434/" in parsed.path,
            "Publisher path not QTD000434 original all cis stats")
    https = urllib.parse.urlunsplit(("https", "ftp.ebi.ac.uk", parsed.path, "", ""))
    return https, output


def parse_qtl_header_rows(raw_partial_gzip):
    require(isinstance(raw_partial_gzip, bytes) and len(raw_partial_gzip) > 16 and
            raw_partial_gzip[:2] == b"\x1f\x8b",
            "Original publisher source is not gzip")
    inflater = zlib.decompressobj(wbits=31)
    text = inflater.decompress(raw_partial_gzip, max_length=5_000_000).decode(
        "utf-8", "strict")
    lines = text.splitlines()
    if not inflater.eof:  # ignore partial trailing row in range response
        lines = lines[:-1]
    require(len(lines) >= 3, "Insufficient original source TSV rows in bounded range")
    header = lines[0].split("\t")
    require(len(header) == len(set(header)) and REQUIRED_QTL_COLUMNS.issubset(header),
            "Original QTD000434 summary statistic columns missing or duplicated")
    reader = csv.DictReader(io.StringIO("\n".join(lines)), delimiter="\t")
    items = [r for _, r in zip(range(200), reader)]
    require(items and all(None not in row and None not in row.values()
                          for row in items), "Malformed publisher QTL preview rows")
    summaries = set()
    maf_samples, an_samples = [], []
    for row in items:
        # These rows may fall OUTSIDE chr10 and are not evidence for 15 candidate genes.
        try:
            an = int(row["an"])
            maf = float(row["maf"])
            beta, se = float(row["beta"]), float(row["se"])
        except (ValueError, TypeError) as exc:
            raise ValueError("Source preview has unparseable AN/MAF/beta/SE") from exc
        require(math.isfinite(beta) and math.isfinite(se) and se > 0 and
                math.isfinite(maf) and 0 < maf <= 0.5 and
                an > 0 and an % 2 == 0,
                "QTL preview has nonfinite data, zero SE or incompatible AN/MAF")
        maf_samples.append(maf)
        an_samples.append(an)
        summaries.add((row["chromosome"], row["gene_id"]))
    return {
        "headers_only": header,
        "preview_rows_checked": len(items),
        "preview_sample_count_min_AN_over_2": min(an_samples) // 2,
        "preview_sample_count_max_AN_over_2": max(an_samples) // 2,
        "preview_rows_variable_AN": len(set(an_samples)) > 1,
        "preview_maf_min": min(maf_samples), "preview_maf_max": max(maf_samples),
        "preview_different_chr_gene_pairs": len(summaries),
        "preview_is_NOT_chr10_15_gene_per_locus_verification": True,
        "actual_15_ROSMAP_gene_QTL_AN_MAF_beta_SE_verified": 0,
        "actual_15_ROSMAP_gene_QTL_sdY_estimated": 0,
        "raw_variant_ids_or_genotypes_exported": False,
    }


def bounded_source_get(url, opener=urllib.request.urlopen):
    uri = urllib.parse.urlsplit(url)
    require(uri.scheme == "https" and uri.hostname == "ftp.ebi.ac.uk" and
            uri.path == ("/pub/databases/spot/eQTL/sumstats/"
                         "QTS000025/QTD000434/QTD000434.all.tsv.gz"),
            "Unexpected original eQTL dataset URL")
    req = urllib.request.Request(url, headers={
        "Range": "bytes=0-"+str(MAX_GZIP_RANGE-1),
        "Accept-Encoding": "identity",
        "User-Agent": "Research-Ischemic-Stroke-Original-QTL-Source-Audit/1.0",
    })
    with opener(req, timeout=40) as res:
        require(res.status == 206, "Server did not support bounded HTTP byte range")
        loc = urllib.parse.urlsplit(res.geturl())
        require(loc.scheme == "https" and loc.hostname == "ftp.ebi.ac.uk",
                "Publisher source redirected to untrusted host")
        r = res.headers.get("Content-Range", "")
        require(re.match(r"^bytes 0-\d+/\d+$",r) is not None,
                "Publisher range response missing unambiguous content-range")
        raw = res.read(MAX_GZIP_RANGE + 1)
        require(len(raw) == MAX_GZIP_RANGE, "Original range read underflow/overflow")
        return raw, {
            "official_HTTPS_url": url,
            "compressed_prefix_SHA256": hashlib.sha256(raw).hexdigest(),
            "range_read_bytes": len(raw), "reported_content_range": r,
            "whole_source_file_verified_or_downloaded": False,
            "source_ETag_observed_not_immutable_pin": res.headers.get("ETag"),
            "source_Last_Modified": res.headers.get("Last-Modified"),
        }


def assemble(manifest, proof, ftp, preview, provenance):
    groups = manifest.get("groups", [])
    ros = [g for g in groups if g["QTL_dataset"] == DATASET_ID]
    require(len(groups) == 36 and len(ros) == 15 and
            proof["QTL_beta_SE_and_MAF_proven_official_fields"] and
            ftp["manifest_dataset_id"] == DATASET_ID and
            preview["preview_rows_checked"] > 0,
            "Official source audit or 36-group manifest mismatch")
    return {
        "schema": "is-p3a-original-ROSMAP-quantitative-coloc-AN-MAF-method-and-QTL-preview-v1",
        "technical_status": "PASS_OFFICIAL_COLOC_EXAMPLE_AND_BOUNDED_REAL_QTL_SOURCE_PREVIEW",
        "scientific_status": "HOLD_GENE_SPECIFIC_FULLCIS_QTL_N_AND_STROKE_GWAS_CASE_CONTROL_MODEL",
        "official_source_provenance": provenance,
        "official_coloc_QTL_AN_MAF_method": proof,
        "exact_QTL_source_dataset": DATASET_ID,
        "exact_QTL_source_study": STUDY_ID,
        "official_full_cis_output_path_identified": ftp["official_manifest_contains_full_cis_sumstats_path"],
        "actual_live_original_QTD000434_first_rows": preview,
        "n_original_gene_tissue_groups": 36, "n_original_ROSMAP_assays": 15,
        "n_gene_specific_QTL_beta_SE_maf_AN_fullcis_evidence_pass": 0,
        "n_approved_QTL_N_via_AN_div2": 0, "n_approved_QTL_sdY_direct": 0,
        "n_original_GIGASTROKE_GWAS_beta_logOR_and_case_fraction_verified": 0,
        "n_formal_coloc_abf_eligible": 0, "n_formal_coloc_abf_executed": 0,
        "n_formal_coloc_susie_eligible": 0, "n_formal_coloc_susie_executed": 0,
        "donor_IDS_or_SNP_effect_rows_saved": False,
        "groups": [
            {"QTL_dataset": g["QTL_dataset"], "gene_id": g["gene_id"],
             "gene_symbol": g["gene_symbol"],
             "scientific_status": "HOLD_FULLCIS_QTL_AND_BINARY_GWAS_PROVENANCE",
             "coloc_abf_eligible": False}
            for g in groups
        ],
    }


def execute(manifest):
    raw, provenance = {}, {}
    for name,url in SOURCES.items():
        print("SOURCE_"+name.upper(),flush=True)
        body, info = SMALL.get(url, {"raw.githubusercontent.com"}, 500000)
        raw[name] = body
        provenance[name] = {"source_url":url,"SHA256":info["SHA256"],
                            "read_bytes":info["read_bytes"]}
    proof = official_method_checks(raw["official_columns"].decode(),
                                   raw["official_coloc_tutorial"].decode(),
                                   raw["official_coloc_implementation"].decode())
    url, ftp = ftp_url_from_manifest(raw["original_public_paths"])
    print("SOURCE_REAL_QTD000434_BOUNDED_RANGE",flush=True)
    prefix, source = bounded_source_get(url)
    preview = parse_qtl_header_rows(prefix)
    provenance["live_original_QTD000434_prefix"] = source
    return assemble(manifest, proof, ftp, preview, provenance)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--manifest",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    result={"schema":"is-p3a-original-ROSMAP-quantitative-coloc-AN-MAF-method-and-QTL-preview-v1",
            "technical_status":"HOLD_BOUNDED_REAL_QTL_SOURCE_OR_METHOD",
            "scientific_status":"HOLD",
            "n_formal_coloc_abf_eligible":0,
            "n_formal_coloc_abf_executed":0}
    try:
        m=json.loads(args.manifest.read_text(encoding="utf-8"))
        result=execute(m)
    except (OSError, ValueError, TypeError, KeyError, UnicodeError,
            TimeoutError, zlib.error) as exc:
        result["error_class_only"]=type(exc).__name__
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "technical_status":result["technical_status"],
        "scientific_status":result["scientific_status"],
        "rows_preview":result.get("actual_live_original_QTD000434_first_rows",{}).get("preview_rows_checked"),
        "n_formal_coloc_abf_executed":0,
        "error_class_only":result.get("error_class_only"),
    },ensure_ascii=False),flush=True)
    return 0 if result["technical_status"].startswith("PASS_") else 1


if __name__=="__main__":
    raise SystemExit(main())
