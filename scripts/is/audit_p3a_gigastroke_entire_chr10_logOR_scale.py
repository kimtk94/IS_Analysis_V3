#!/usr/bin/env python3
"""Read entire original GIGASTROKE AIS files, audit ALL chr10 beta/log(OR)/CI.

Ancestry-specific GRCh37 files are downloaded only to disposable local disk.
Original compressed size/MD5/SHA256 must match publisher-pinned evidence.
Read all rows to confirm original total, but only aggregate chr10 SNPs. No
individual row, rsID, allele or effect values saved; no replication/Coloc.
"""
from __future__ import annotations

import argparse
import collections
import csv
import gzip
import importlib
import json
import math
import tempfile
from pathlib import Path

SOURCE=importlib.import_module("scripts.is.audit_p3a_chr10_original_publisher_gwas")
PINS=importlib.import_module("scripts.is.p3a_publisher_md5_gate")
ORIGINAL=importlib.import_module("scripts.is.audit_p3a_gigastroke_original_beta_se_p_scale")

REFERENCE_CHR10_CONTEXT={"EAS":3913,"EUR":4397}
FULL_SOURCE_ROWS={"EAS":6_725_216,"EUR":7_482_032}
MAX_PREVIEW_COLUMNS=50
MAX_FILE_ROWS=9_000_000
GWAS_COLS=("chromosome","base_pair_location","beta","standard_error",
            "p_value","effect_allele","other_allele","effect_allele_frequency",
            "odds_ratio","ci_lower","ci_upper")


def require(ok, reason):
    if not ok:raise ValueError(reason)


def empty_counts():
    return collections.Counter({
       "entire_chr10_rows":0,"chr10_rows_with_parsed_beta_se":0,
       "chr10_rows_OR_and_CI_populated":0,
       "chr10_rows_OR_missing_or_invalid":0,
       "chr10_rows_CI_missing_or_invalid":0,
       "chr10_rows_OR_not_exp_beta":0,
       "chr10_rows_CI_not_beta_plusminus_1_96_SE":0,
       "chr10_rows_p_zero_underflow":0,
       "chr10_rows_z_p_comparable":0,
       "chr10_rows_z_p_discrepant_log10_gt_0_1":0,
       "chr10_rows_invalid_beta_SE_p_or_alleles":0,
       "chr10_source_context_valid_2Mb":0,
       "chr10_source_context_invalid_2Mb":0,
       "chr10_target_legacy_interval_valid":0,
       "chr10_source_context_OR_consistent":0,
       "chr10_source_context_CI_consistent":0,
    })


def audit_all_chr10_rows(path, ancestry, min_chr10=1000, expected_total=True,
                         min_total_rows=1_000_000):
    """Streaming full gzip; fail-scientifically on bad OR/CI, not fake observed N."""
    require(ancestry in PINS.SOURCES, "Unknown exact source ancestry")
    counts=empty_counts()
    total=0
    max_or_error=0.0
    max_ci_error=0.0
    max_p_error=0.0
    se_rows_not_finite=0
    with gzip.open(path,"rt",encoding="utf-8-sig",newline="") as stream:
        reader=csv.DictReader(stream,delimiter="\t")
        fields=reader.fieldnames or []
        require(len(fields)==len(set(fields)) and set(GWAS_COLS).issubset(fields)
                and len(fields)<=MAX_PREVIEW_COLUMNS,
                "Original GWAS file missing populated OR/CI/beta/SE/allele columns")
        for row in reader:
            total+=1
            require(total<=MAX_FILE_ROWS,"Original GWAS source total rows over safety limit")
            chrom=(row.get("chromosome") or "").strip().upper().removeprefix("CHR")
            if chrom!="10":
                continue
            counts["entire_chr10_rows"]+=1
            require(None not in row and None not in row.values(),
                    "Original publisher chr10 row has incorrect tabular structure")
            try:
                pos=int(row["base_pair_location"])
                beta=float(row["beta"])
                se=float(row["standard_error"])
                p=float(row["p_value"])
                ea=row["effect_allele"].strip().upper()
                oa=row["other_allele"].strip().upper()
                eaf=float(row["effect_allele_frequency"])
            except (TypeError,ValueError,KeyError):
                counts["chr10_rows_invalid_beta_SE_p_or_alleles"]+=1
                continue
            if not (pos>0 and math.isfinite(beta) and math.isfinite(se)
                    and se>0 and math.isfinite(p) and 0<=p<=1
                    and math.isfinite(eaf) and 0<=eaf<=1
                    and ea in "ACGT" and oa in "ACGT" and ea!=oa
                    and len(ea)==len(oa)==1):
                counts["chr10_rows_invalid_beta_SE_p_or_alleles"]+=1
                continue
            counts["chr10_rows_with_parsed_beta_se"]+=1
            in_context=(SOURCE.REGION["source_context_start"]<=pos<=
                        SOURCE.REGION["source_context_end"])
            if in_context:
                counts["chr10_source_context_valid_2Mb"]+=1
                if (SOURCE.REGION["original_region_start"]<=pos<=
                    SOURCE.REGION["original_region_end"]):
                    counts["chr10_target_legacy_interval_valid"]+=1
            try:
                orval=ORIGINAL.optional_positive_float(row["odds_ratio"])
                low=ORIGINAL.optional_positive_float(row["ci_lower"])
                high=ORIGINAL.optional_positive_float(row["ci_upper"])
            except (TypeError,ValueError):
                counts["chr10_rows_OR_missing_or_invalid"]+=1
                counts["chr10_rows_CI_missing_or_invalid"]+=1
                continue
            or_consistent=False
            ci_consistent=False
            if orval is None:
                counts["chr10_rows_OR_missing_or_invalid"]+=1
            else:
                diff=abs(math.log(orval)-beta)
                max_or_error=max(max_or_error,diff)
                if diff>0.015:
                    counts["chr10_rows_OR_not_exp_beta"]+=1
                else:or_consistent=True
            if low is None or high is None or low>=high:
                counts["chr10_rows_CI_missing_or_invalid"]+=1
            else:
                diff=max(
                    abs(math.log(low)-(beta-1.96*se)),
                    abs(math.log(high)-(beta+1.96*se)))
                max_ci_error=max(max_ci_error,diff)
                if diff>0.04:
                    counts["chr10_rows_CI_not_beta_plusminus_1_96_SE"]+=1
                else:ci_consistent=True
            if or_consistent and ci_consistent:
                counts["chr10_rows_OR_and_CI_populated"]+=1
            if in_context:
                counts["chr10_source_context_OR_consistent"]+=or_consistent
                counts["chr10_source_context_CI_consistent"]+=ci_consistent
            if p==0:
                counts["chr10_rows_p_zero_underflow"]+=1
            else:
                pz=math.erfc(abs(beta/se)/math.sqrt(2))
                if pz>0:
                    err=abs(math.log10(pz)-math.log10(p))
                    max_p_error=max(max_p_error,err)
                    counts["chr10_rows_z_p_comparable"]+=1
                    if err>0.1:counts["chr10_rows_z_p_discrepant_log10_gt_0_1"]+=1
    require(total>=min_total_rows and counts["entire_chr10_rows"]>=min_chr10,
            "Original source file or chr10 incomplete")
    require(not expected_total or total==FULL_SOURCE_ROWS[ancestry],
            "Original file variant count different from verified published source")
    require(not expected_total or counts["chr10_source_context_valid_2Mb"]==
            REFERENCE_CHR10_CONTEXT[ancestry],
            "Original chr10 context count no longer matches immutable prior original scan")
    count=dict(sorted(counts.items()))
    all_or=(counts["chr10_rows_with_parsed_beta_se"]==
            counts["chr10_rows_OR_and_CI_populated"]==
            counts["entire_chr10_rows"]
            and counts["chr10_rows_OR_not_exp_beta"]==0
            and counts["chr10_rows_CI_not_beta_plusminus_1_96_SE"]==0)
    return {
        "original_publisher_accession":PINS.SOURCES[ancestry]["accession"],
        "source_build":"GRCh37",
        "original_complete_GWAS_rows_scanned":total,
        "original_total_chr10_rows":counts["entire_chr10_rows"],
        "original_legacy_chr10_2MB_context_rows_matches_prior":(
            counts["chr10_source_context_valid_2Mb"]==REFERENCE_CHR10_CONTEXT[ancestry]),
        "source_beta_lnOR_and_CI_consistency_chr10":(
            "PASS_ALL_OBSERVED_CHR10" if all_or else "HOLD_SOME_CHR10_SOURCE_ROWS"),
        "all_source_chr10_OR_and_CI_logscale_consistent":all_or,
        "max_abs_beta_minus_ln_original_OR_all_chr10":max_or_error,
        "max_abs_log_OR_CI_vs_beta_pm_1_96_SE_chr10":max_ci_error,
        "max_abs_log10_source_p_vs_beta_SE_chr10":max_p_error,
        "aggregates":count,
        "effect_allele_and_other_allele_checked":True,
        "original_allele_orientation_preserved_not_rewritten":True,
        "all_original_SNP_effect_rows_or_IDs_emitted":False,
        "exact_source_file_cohort_case_fraction_verified":False,
        "GWAS_per_SNP_effective_N_verified":False,
        "BBJ_independent_replication_approved":False,
        "source_beta_units_approved_for_automated_formal_coloc":False,
        "scientific_status":"HOLD_CASE_FRACTION_AND_QTL_MODEL_N",
    }


def collect_original(ancestry,groups):
    require(len(groups)==36 and len({(x["QTL_dataset"],x["gene_id"]) for x in groups})==36,
            "Original 36 assay evidence contract broken")
    with tempfile.TemporaryDirectory(prefix="is_p3a_original_full_chr10_") as tmp:
        dest=Path(tmp)/"official_original_ais.tsv.gz"
        print("PUBLISHER_ORIGINAL_"+ancestry+"_DOWNLOAD_AND_VERIFY_MD5_SHA256",flush=True)
        source=SOURCE.stream_source(ancestry,dest)
        require(source["publisher_md5_and_prior_sha256_verified"] is True,
                "Original source full digest failed")
        print("PUBLISHER_ORIGINAL_"+ancestry+"_WHOLE_FILE_SCAN_AND_ALL_CHR10_BETA_OR_CI",flush=True)
        result=audit_all_chr10_rows(dest,ancestry)
    return {
        "schema":"is-p3a-original-GIGASTROKE-entire-chr10-EAS-EUR-logOR-scale-v1",
        "technical_status":"PASS_PUBLISHER_WHOLE_GWAS_MD5_SHA256_AND_ALL_CHR10_ROWS",
        "scientific_status":"HOLD_CASE_FRACTION_QTL_REGRESSION_N_AND_EXACT_REGIONAL_HARMONIZATION",
        "ancestry":ancestry,
        "source_publisher_integrity":source,
        "entire_chr10_GWAS_scale_audit":result,
        "n_original_assays":36,
        "n_original_distinct_genes":20,
        "n_original_independent_stroke_loci":1,
        "n_formal_coloc_abf_eligible":0,
        "n_formal_coloc_abf_executed":0,
        "n_formal_coloc_susie_eligible":0,
        "n_formal_coloc_susie_executed":0,
        "original_2022_case_fraction_s_not_approved":True,
        "sample_or_SNP_effect_records_exported":False,
        "target_assay_keys":[
            {"QTL_dataset":x["QTL_dataset"],"gene_id":x["gene_id"],
             "gene_symbol":x["gene_symbol"],"coloc_eligible":False}
            for x in groups
        ]
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--ancestry",required=True,choices=("EAS","EUR"))
    p.add_argument("--manifest",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    result={"schema":"is-p3a-original-GIGASTROKE-entire-chr10-EAS-EUR-logOR-scale-v1",
            "technical_status":"HOLD_ORIGINAL_SOURCE_OR_CHR10_QC",
            "scientific_status":"HOLD",
            "ancestry":args.ancestry,
            "n_formal_coloc_abf_eligible":0,"n_formal_coloc_abf_executed":0}
    try:
        groups=SOURCE.checked_manifest(args.manifest)
        result=collect_original(args.ancestry,groups)
    except (OSError,ValueError,TypeError,KeyError,EOFError,OverflowError) as exc:
        result["error_class_only"]=type(exc).__name__
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",
                        encoding="utf-8")
    print(json.dumps({
        "ancestry":args.ancestry,
        "technical_status":result["technical_status"],
        "scientific_status":result["scientific_status"],
        "total_original_chr10_SNPs":result.get("entire_chr10_GWAS_scale_audit",{}).get(
            "original_total_chr10_rows"),
        "all_chr10_original_OR_and_CI_valid":result.get("entire_chr10_GWAS_scale_audit",{}).get(
            "all_source_chr10_OR_and_CI_logscale_consistent"),
        "error_class_only":result.get("error_class_only"),
        "coloc_executed":0},ensure_ascii=False),flush=True)
    return 0 if result["technical_status"].startswith("PASS_") else 1


if __name__=="__main__":
    raise SystemExit(main())
