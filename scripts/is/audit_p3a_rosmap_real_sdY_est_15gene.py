#!/usr/bin/env python3
"""Original 15 ROSMAP gene QTL: exact coloc::sdY.est estimator and diagnostics.

Inputs: original QTD000434 all-association source, 2023 Ensembl105 phenotype
metadata and 36 original assay manifest. The 2023 locus +/-1Mb is an audit
window, not certified original Release 7 runtime parameters.

Mathematics exactly reproduces coloc/R/claudia.R:
  oneover=1/vbeta; nvx=2*N*maf*(1-maf)
  cf=coef(lm(nvx ~ oneover - 1))[["oneover"]]
  estimated_sdY=sqrt(cf)
N=AN/2=560 is a *source allele-count-based scenario*, NOT a final regression N.
All formal coloc and scientific approvals stay HOLD. Only per-gene aggregates
are exported, never SNP-level β/SE/MAF data or donor identifiers.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import importlib
import json
import math
import statistics
from pathlib import Path

QTL = importlib.import_module("scripts.is.audit_p3a_rosmap_15_gene_single_tabix_qtl")
RANGE = importlib.import_module("scripts.is.qtl_original_indexed_http_range")
SOURCE_URL = ("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
              "QTS000025/QTD000434/QTD000434.all.tsv.gz")
PINNED_SOURCE_INDEX_SHA256 = (
    "6cbf12826e19a61594662972a1205ba38cdd3071984e51c24261bff8ace9806a")
PINNED_SOURCE_RANGE_SHA256 = (
    "5a54c1faadb0530ecdca37497f82f2e28f1003c000d7a24bb9fa8422de3fa2f1")
EXPECTED_CHR10_UNION_REGION = "10:101854259-105353833"
EXPECTED_SOURCE_ROWS_IN_UNION = 367_753
EXPECTED_TARGET_SOURCE_ROWS = 100_822
EXPECTED_15_UNIQUE_GENE_SNP_RECORDS = 89_844
SOURCE_AN = 1120
AN_DIV2_PROVISIONAL_N = SOURCE_AN // 2
MIN_GENE_INPUT_ROWS = 3000
MIN_STRATUM_ROWS = 30
MAX_QTL_ROWS = 3_500_000
EXACT_15_RECORD_COUNTS = {
    "ENSG00000065613":5854, "ENSG00000065621":6093,
    "ENSG00000107954":5961, "ENSG00000120051":6213,
    "ENSG00000138172":6018, "ENSG00000148798":6031,
    "ENSG00000148834":6103, "ENSG00000148841":6202,
    "ENSG00000148842":5809, "ENSG00000166275":5728,
    "ENSG00000173915":6041, "ENSG00000197748":6099,
    "ENSG00000214435":5723, "ENSG00000228261":6163,
    "ENSG00000272912":5806,
}


def require(value, message):
    if not value:
        raise ValueError(message)


def estimate_coloc_sdY(allele_maf_and_se, n):
    """Reproduce coloc R lm(nvx ~ oneover - 1) exactly by OLS coefficient."""
    require(isinstance(n,(int,float)) and math.isfinite(n) and 1<n<100_000,
            "N hypothesis must be an explicitly declared finite value")
    require(len(allele_maf_and_se)>=2,
            "Need at least two valid original source variants")
    x=[]
    y=[]
    per_variant_sd=[]
    for maf,se in allele_maf_and_se:
        require(math.isfinite(maf) and math.isfinite(se) and
                0<maf<=0.5 and se>0, "Source MAF/SE invalid")
        vbeta=se*se
        oneover=1/vbeta
        nvx=2*n*maf*(1-maf)
        require(all(map(math.isfinite,(vbeta,oneover,nvx))) and
                oneover>0 and nvx>0, "Nonfinite coloc sdY estimator inputs")
        x.append(oneover)
        y.append(nvx)
        per_variant_sd.append(math.sqrt(vbeta*nvx))
    denominator=math.fsum(v*v for v in x)
    numerator=math.fsum(a*b for a,b in zip(x,y))
    require(math.isfinite(denominator) and math.isfinite(numerator)
            and denominator>0,"Coloc zero-intercept regression unstable")
    slope=numerator/denominator
    require(math.isfinite(slope) and slope>0,
            "Estimated coloc sdY variance not positive finite")
    sdY=math.sqrt(slope)
    residuals=[(u-slope*v) for u,v in zip(y,x)]
    syy=math.fsum(u*u for u in y)
    require(syy>0,"Invalid QTL variance regression scale")
    zero_intercept_r2=1-math.fsum(v*v for v in residuals)/syy
    return {
        "sdY_est_coloc_zero_intercept_OLS":sdY,
        "sdY_estimate_squared_OLS_coefficient":slope,
        "source_rows_used_for_estimator":len(x),
        "zero_intercept_OLS_R_squared_diagnostic":zero_intercept_r2,
        "median_single_SNP_implied_sdY_diagnostic":statistics.median(per_variant_sd),
        "p10_single_SNP_implied_sdY_diagnostic":quantile(per_variant_sd,.10),
        "p90_single_SNP_implied_sdY_diagnostic":quantile(per_variant_sd,.90),
        "this_is_original_normalized_expression_sdY_direct_measurement":False,
        "assumed_effective_N_approved":False,
    }


def quantile(values, p):
    a=sorted(values)
    require(a and 0<=p<=1,"Invalid diagnostic quantile")
    k=(len(a)-1)*p
    lo=math.floor(k)
    hi=math.ceil(k)
    return a[lo]+(a[hi]-a[lo])*(k-lo)


def source_record_qc(lines, manifest, region, gene_metadata):
    """Strict original source per-gene dedup, then estimate candidate sdY."""
    assay_rows=[g for g in manifest["groups"] if g["QTL_dataset"]==QTL.DATASET]
    require(len(manifest["groups"])==36 and len(assay_rows)==15,
            "Expected exactly 36 true assays and 15 ROSMAP groups")
    targets={g["gene_id"]:g for g in assay_rows}
    require(set(targets)==set(EXACT_15_RECORD_COUNTS)==set(gene_metadata),
            "ROSMAP 15 original ENSG identity list differs")
    raw_gene_records=collections.Counter()
    dedup={gid:{} for gid in targets}
    scanned=0
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        parts=line.rstrip("\r\n").split("\t")
        if tuple(parts)==QTL.COLUMNS:
            continue
        require(len(parts)==len(QTL.COLUMNS),
                "Original QTL field count differs from approved source header")
        scanned+=1
        require(scanned<=MAX_QTL_ROWS,
                "Original source stream exceeded bounded source row limit")
        row=dict(zip(QTL.COLUMNS,parts))
        gid=QTL.stable_gene_id(row["gene_id"])
        if gid not in targets:
            continue
        require(QTL.stable_gene_id(row["molecular_trait_id"])==gid,
                "QTL gene phenotype and association gene identifier differ")
        variant,pos,alt=QTL.source_variant_identity(row,region)
        require(abs(pos-gene_metadata[gid]["phenotype_pos"])<=QTL.AUDIT_CIS_BP,
                "Original phenotype window changed")
        an,maf,beta,se,p=QTL.numeric_fields(row)
        require(an==SOURCE_AN, "One or more original gene QTL AN !=1120")
        raw_gene_records[gid]+=1
        key=(an,maf,beta,se,p,row["ref"],row["alt"],row["type"])
        if variant in dedup[gid]:
            require(dedup[gid][variant][0]==key,
                    "Conflicting duplicate QTL variant source association")
            continue
        dedup[gid][variant]=(key,(maf,se),pos)
        require(len(dedup[gid])<=QTL.MAX_GENE_UNIQUE_VARIANTS,
                "Original gene variant count exceeds source bound")
    require(scanned==EXPECTED_SOURCE_ROWS_IN_UNION and
            sum(raw_gene_records.values())==EXPECTED_TARGET_SOURCE_ROWS,
            "Original QTL publisher source union rows changed from main #323")
    report=[]
    for gid,g in targets.items():
        records=dedup[gid]
        require(len(records)==EXACT_15_RECORD_COUNTS[gid] and
                len(records)>=MIN_GENE_INPUT_ROWS,
                "Original per-gene QTL variant completeness differs from verified main")
        ordered=sorted(records.items(),key=lambda x:x[0])
        vals=[v[1][1] for v in ordered] # each = (maf, SE)
        full=estimate_coloc_sdY(vals,AN_DIV2_PROVISIONAL_N)
        low=[z for z in vals if z[0]<.1]
        high=[z for z in vals if z[0]>=.1]
        low_result=(estimate_coloc_sdY(low,AN_DIV2_PROVISIONAL_N)
                    if len(low)>=MIN_STRATUM_ROWS else None)
        high_result=(estimate_coloc_sdY(high,AN_DIV2_PROVISIONAL_N)
                     if len(high)>=MIN_STRATUM_ROWS else None)
        # Odd-even deterministic by sorted exact GRCh38 REF/ALT variant ID.
        # NOT independent replication: neighboring variants are in LD.
        even=estimate_coloc_sdY(vals[::2],AN_DIV2_PROVISIONAL_N)
        odd=estimate_coloc_sdY(vals[1::2],AN_DIV2_PROVISIONAL_N)
        sensitivity={
            str(n): estimate_coloc_sdY(vals,n)["sdY_est_coloc_zero_intercept_OLS"]
            for n in (540,560,580)
        }
        baseline=full["sdY_est_coloc_zero_intercept_OLS"]
        require(math.isclose(sensitivity["560"],baseline,rel_tol=1e-13),
                "Coloc estimator changed across nominal N scenario")
        # N scaling is inherent to method, NOT newly observed sample counts.
        for n in (540,580):
            require(math.isclose(sensitivity[str(n)]/baseline,
                                 math.sqrt(n/560),rel_tol=1e-12),
                    "Zero-intercept coloc sdY scaling with assumed N incorrect")
        ratio=(low_result["sdY_est_coloc_zero_intercept_OLS"]/
               high_result["sdY_est_coloc_zero_intercept_OLS"]
               if low_result and high_result else None)
        split_relative=abs(even["sdY_est_coloc_zero_intercept_OLS"]-
                           odd["sdY_est_coloc_zero_intercept_OLS"])/baseline
        report.append({
            "gene_id":gid,"gene_symbol":g["gene_symbol"],
            "QTL_dataset":QTL.DATASET,
            "original_source_raw_QTL_rows":raw_gene_records[gid],
            "original_source_unique_QTL_variants":len(records),
            "source_exact_AN_across_all_variants":SOURCE_AN,
            "source_AN_div2_nominal_N_hypothesis":AN_DIV2_PROVISIONAL_N,
            "original_qtl_run_gene_N_independently_verified":False,
            "source_MAF_min":min(z[0] for z in vals),
            "source_MAF_max":max(z[0] for z in vals),
            "source_SE_min":min(z[1] for z in vals),
            "source_SE_max":max(z[1] for z in vals),
            "source_QTL_beta_SE_MAF_AN_all_numeric_and_allelic":True,
            "official_coloc_sdY_estimator_N560_hypothesis":full,
            "MAF_below_0_1_rows":len(low),
            "MAF_at_least_0_1_rows":len(high),
            "MAF_below_0_1_sdY_est":(low_result or {}).get(
                "sdY_est_coloc_zero_intercept_OLS"),
            "MAF_at_least_0_1_sdY_est":(high_result or {}).get(
                "sdY_est_coloc_zero_intercept_OLS"),
            "MAF_below_vs_above_0_1_estimate_ratio":ratio,
            "odd_even_sorted_variant_split_sdY_relative_difference":
                split_relative,
            "odd_even_split_is_NOT_independent_due_to_LD":True,
            "hypothetical_N_540_560_580_sdY_sensitivity_UNAPPROVED":sensitivity,
            "original_expression_matrix_sdY_directly_measured":False,
            "original_R7_qtlmap_run_cis_config_verified":False,
            "sdY_estimate_approved_for_formal_coloc":False,
            "scientific_status":"HOLD_GENEWISE_N_AND_FULLCIS_ORIGINAL_CONFIG_AND_BINARY_GWAS_CC",
            "formal_coloc_eligible":False,
        })
    return {
        "real_original_all_QTL_rows_in_union_region":scanned,
        "real_original_target_gene_QTL_source_rows":sum(raw_gene_records.values()),
        "real_original_gene_wise_unique_SNP_sum":sum(len(x) for x in dedup.values()),
        "n_genes_real_sdY_OLS_estimated_N560_scenario":len(report),
        "n_genes_original_normalized_expression_sdY_directly_measured":0,
        "n_genes_N560_regression_sample_size_approved":0,
        "n_genes_estimator_approved_for_coloc":0,
        "n_genes_with_both_MAF_strata_estimable":sum(
            x["MAF_below_vs_above_0_1_estimate_ratio"] is not None for x in report),
        "no_variant_effect_rows_or_sample_IDs_exported":True,
        "groups":report,
    }


def audit(manifest):
    source_metadata,meta=QTL.source_annotations()
    target_ids={g["gene_id"] for g in manifest["groups"]
                if g["QTL_dataset"]==QTL.DATASET}
    require(target_ids==set(EXACT_15_RECORD_COUNTS),
            "Unexpected genuine 15 ROSMAP gene IDs")
    phenotypes=QTL.source_gene_metadata(source_metadata,target_ids)
    region=QTL.one_combined_region(phenotypes,target_ids)
    require(region["region"]==EXPECTED_CHR10_UNION_REGION and
            region["R7_actual_qtlmap_run_cis_parameter_independently_verified"] is False,
            "2023 QTL source union locus drift or false R7 run approval")
    print("SOURCE_ONLY_ONE_ORIGINAL_2023_PHENOTYPE_MD5_AND_TBI_RANGE",flush=True)
    lines,range_prov=RANGE.source_indexed_range(SOURCE_URL,region)
    require(range_prov["tabix_index_SHA256"]==PINNED_SOURCE_INDEX_SHA256
            and range_prov["raw_source_metadata"][
                "publisher_compressed_range_SHA256"]==PINNED_SOURCE_RANGE_SHA256,
            "Original publisher QTL source TBI/byte-range SHA256 changed from #323")
    print("SOURCE_15_REAL_ROSMAP_QTL_GENEWISE_COLOC_SDY_EST_AND_MAF_SENSITIVITY",
          flush=True)
    result=source_record_qc(lines,manifest,region,phenotypes)
    require(result["n_genes_real_sdY_OLS_estimated_N560_scenario"]==15 and
            result["real_original_gene_wise_unique_SNP_sum"]==
                EXPECTED_15_UNIQUE_GENE_SNP_RECORDS,
            "15-gene original QTL records differ from previous source evidence")
    return {
        "schema":"is-p3a-real-ROSMAP-15gene-coloc-sdY-est-N560-MAF-QC-20261123-v1",
        "technical_status":"PASS_ORIGINAL_COLLOC_SDY_FORMULA_AND_15GENE_SOURCE_QC",
        "scientific_status":"HOLD_ORIGINAL_N_AND_SDY_DIRECT_VALIDATION_AND_GWAS_CASE_FRACTION",
        "original_coloc_source_reference":{
            "repo":"chr1swallace/coloc",
            "R_file":"R/claudia.R",
            "github_original_blob_sha1":"043783ac5ddec3477ff3b91e3e1faf271c51ce09",
            "function":"sdY.est(vbeta,maf,n)",
            "method":"lm(nvx ~ oneover - 1), oneover=1/vbeta; nvx=2*N*maf*(1-maf)",
            "sdY_is_sqrt_zero_intercept_OLS_slope":True,
            "variance_estimate_not_observed_original_expression_standard_deviation":True,
        },
        "original_source_provenance":{
            "original_publisher_metadata_Zenodo_MD5":QTL.Z_MD5,
            "original_publisher_QTL_TBI_SHA256":range_prov["tabix_index_SHA256"],
            "original_publisher_QTL_range_SHA256":
                range_prov["raw_source_metadata"]["publisher_compressed_range_SHA256"],
            "original_source_region":region,
            "original_source_number_EBI_HTTP_requests":range_prov[
                "requests_to_original_EBI_source"],
            "original_full_4GB_QTL_file_downloaded":False,
        },
        "original_15_gene_sdY_diagnostics":result,
        "original_ROSMAP_QTL_AN_div2_regression_N_approval":False,
        "original_ROSMAP_normalized_expression_sample_sdY_measured":False,
        "original_GIGASTROKE_binary_case_fraction_approved":False,
        "n_original_genuine_assays":36,
        "n_original_ROSMAP_gene_assays":15,
        "n_gene_tissue_assays_formal_coloc_eligible":0,
        "n_coloc_abf_executed":0,
        "n_coloc_susie_executed":0,
        "no_raw_source_SNP_effects_MAF_SE_or_donor_IDs_exported":True,
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--manifest",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    out={
        "schema":"is-p3a-real-ROSMAP-15gene-coloc-sdY-est-N560-MAF-QC-20261123-v1",
        "technical_status":"HOLD_ORIGINAL_QTL_SOURCE_OR_SD_EST",
        "scientific_status":"HOLD",
        "n_gene_tissue_assays_formal_coloc_eligible":0,
        "n_coloc_abf_executed":0,
    }
    try:
        manifest=json.loads(args.manifest.read_text(encoding="utf-8"))
        out=audit(manifest)
    except (OSError,ValueError,KeyError,TypeError,UnicodeError,TimeoutError,
            OverflowError,ArithmeticError) as error:
        out["error_class_only"]=type(error).__name__
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n",
                        encoding="utf-8")
    print(json.dumps({
        "technical_status":out["technical_status"],
        "scientific_status":out["scientific_status"],
        "real_genes_estimated_N560_hypothesis":out.get(
            "original_15_gene_sdY_diagnostics",{}).get(
            "n_genes_real_sdY_OLS_estimated_N560_scenario"),
        "n_formal_coloc_executed":0,
        "error_class_only":out.get("error_class_only"),
    },ensure_ascii=False),flush=True)
    return 0 if out["technical_status"].startswith("PASS_") else 1


if __name__=="__main__":
    raise SystemExit(main())
