#!/usr/bin/env python3
"""Compute original coloc sdY.est on all 21 genuine GTEx Aorta/Cortex gene-QTLs.

Original eQTL Catalogue datasets QTD000131 (GTEx artery aorta) and QTD000171
(GTEx brain cortex) are read with ONE original TBI+ONE BGZF Range per dataset.
Validate all source region rows against earlier original indexed source SHA256
and the preverified per-assay strict nonpal variants and AN/2 source census.

R coloc sdY.est: oneover=1/SE^2; nvx=2*N*MAF*(1-MAF);
sdY=sqrt(sum(oneover*nvx)/sum(oneover^2)).
N=387 or 205 is a CATALOGUE-N SCENARIO, not certified regression model N.
Variant AN/2 is heterogeneous across all 21 gene assays.
No original per-SNP data exported. No formal coloc run or posterior fabricated.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import importlib
import json
import math
from pathlib import Path

ORIGINAL=importlib.import_module("scripts.is.qtl_original_indexed_http_range")
QTL=importlib.import_module("scripts.is.audit_p3a_indexed_qtl_real_chr10")
MATCH=importlib.import_module("scripts.is.audit_p3a_36_real_gwas_qtl_reciprocal_snp")
CIS=importlib.import_module("scripts.is.audit_p3a_official_cis_union_36")
OLS=importlib.import_module("scripts.is.audit_p3a_rosmap_real_sdY_est_15gene")

DATASETS={
  "QTD000131":{
    "url":"https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000131/QTD000131.all.tsv.gz",
    "catalogue_n_hypothesis":387,
    "region_source_rows":241786,
    "region_stream_SHA256":"451ba8ba15493986f431b1751aff9ea512ab3dbc04ef41e38f95217d49cf3f53",
    "expected_assays":10,
    "expected_nonpal_source_assay_records":44693,
    "gene_variant_counts":{
       "ENSG00000065613":4537, "ENSG00000065621":4851,
       "ENSG00000076685":4287, "ENSG00000107954":4383,
       "ENSG00000138172":4419, "ENSG00000148834":4872,
       "ENSG00000148843":4386, "ENSG00000156384":4755,
       "ENSG00000166275":3748, "ENSG00000273108":4455
    }
  },
  "QTD000171":{
    "url":"https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000171/QTD000171.all.tsv.gz",
    "catalogue_n_hypothesis":205,
    "region_source_rows":264947,
    "region_stream_SHA256":"5964ba14a274f59247eebf88602da21706fdaea6c6e5df57990a3157f7854953",
    "expected_assays":11,
    "expected_nonpal_source_assay_records":45570,
    "gene_variant_counts":{
       "ENSG00000065621":4477, "ENSG00000076685":4045,
       "ENSG00000138172":4115, "ENSG00000148795":3527,
       "ENSG00000148834":4501, "ENSG00000148841":4546,
       "ENSG00000166275":3565, "ENSG00000197748":4501,
       "ENSG00000228261":4521, "ENSG00000272912":3637,
       "ENSG00000273108":4135
    }
  },
}
MIN_GENE_VARIANTS=3000
MIN_FULL_AN_VARIANTS=3000
MIN_MAF_STRATUM=30
MAX_SOURCE_REGION_ROWS=500000


def require(cond, explanation):
    if not cond:
        raise ValueError(explanation)


def original_region():
    x=CIS.EXPECTED_WINDOW
    require(x=={"chrom":"10","start":101837472,"end":105353833},
            "Pinned original 20-gene source reference union changed")
    return {"start":x["start"],"end":x["end"],"region":
            f"10:{x['start']}-{x['end']}"}


def manifest_assays(manifest,dataset):
    groups=manifest.get("groups",[])
    require(len(groups)==36 and
            len({(x["QTL_dataset"],x["gene_id"]) for x in groups})==36 and
            len({x["gene_id"] for x in groups})==20 and
            sum(x["QTL_dataset"]=="QTD000434" for x in groups)==15,
            "Original 36 source assay manifest identity changed")
    d=DATASETS[dataset]
    selected={x["gene_id"]:x for x in groups if x["QTL_dataset"]==dataset}
    require(len(selected)==d["expected_assays"] and
            set(selected)==set(d["gene_variant_counts"]),
            "Original GTEx 21 assay source identity changed")
    return selected


def scan_source(lines,dataset,manifest,*,expected_source_counts=True):
    require(dataset in DATASETS,"Only two pinned original GTEx QTL datasets allowed")
    d=DATASETS[dataset]
    target=manifest_assays(manifest,dataset)
    region=original_region()
    bins={g:{} for g in target}
    scanned=0
    # The ORIGINAL 2026-11-22 source digest protocol includes the 19-column
    # physical header BEFORE tabix source rows (even though tabix excludes it).
    # Match the prior independently verified stream hash byte-for-byte.
    digest=hashlib.sha256(("\t".join(QTL.OFFICIAL_COLUMN_ORDER)+"\n").encode("utf-8"))
    common_counts=collections.Counter()
    filtered_candidate_rows=0
    source_an_missing=0
    for original_line in lines:
        if isinstance(original_line,bytes):
            original_line=original_line.decode("utf-8")
        if not original_line.strip() or original_line.startswith("#"):
            continue
        parts=original_line.rstrip("\r\n").split("\t")
        if tuple(parts)==QTL.OFFICIAL_COLUMN_ORDER:
            continue
        require(len(parts)==len(QTL.OFFICIAL_COLUMN_ORDER),
                "Original eQTL publisher TSV column order or count drift")
        row=dict(zip(QTL.OFFICIAL_COLUMN_ORDER,parts))
        try:
            pos=int(row["position"])
        except (TypeError,ValueError) as exc:
            raise ValueError("Original variant position malformed") from exc
        chrom=row["chromosome"].removeprefix("chr")
        if chrom!="10" or not region["start"]<=pos<=region["end"]:
            continue
        scanned+=1
        require(scanned<=MAX_SOURCE_REGION_ROWS,"Original source GTEx interval row cap")
        digest.update((original_line.rstrip("\r\n")+"\n").encode("utf-8"))
        gid=row["gene_id"].split(".")[0]
        if gid not in target:
            continue
        item=MATCH.qtl_record(row)
        if item is None:
            continue
        require(item["molecular_trait_id"].split(".")[0]==gid,
                "Source QTL molecular phenotype ID drift")
        filtered_candidate_rows+=1
        try:
            an=int(row["an"]);ac=float(row["ac"])
        except (TypeError,ValueError,OverflowError) as exc:
            source_an_missing+=1
            raise ValueError("Original GTEx QTL allele number missing/invalid") from exc
        require(an>0 and an%2==0 and math.isfinite(ac) and
                0<=ac<=an and an//2<=d["catalogue_n_hypothesis"],
                "Source GTEx genotype AN/2 unphysical or exceeds original study N")
        n=an//2
        key=(item["pos38"],item["ref"],item["alt"],row["molecular_trait_id"])
        datum=(item["maf"],item["se"],n,item["beta"],row["pvalue"])
        prior=bins[gid].get(key)
        if prior is not None:
            require(prior==datum,
                    "Original GTEx duplicate variant differs in QTL beta SE MAF AN or p")
            continue
        bins[gid][key]=datum
        common_counts[gid]+=1
    # Only safe whole-stream aggregate diagnostics, not SNP IDs or values.
    print("ORIGINAL_GTEX_"+dataset+"_REGION_AGGREGATE_ROWS="+str(scanned)+
          " BASELINE_ROWS="+str(d["region_source_rows"])+
          " STREAM_SHA256="+digest.hexdigest()+
          " BASELINE_SHA256="+d["region_stream_SHA256"],flush=True)
    if expected_source_counts:
        require(scanned==d["region_source_rows"] and
                digest.hexdigest()==d["region_stream_SHA256"],
                "Original GTEx complete 20-gene region stream SHA256/count changed")
    require(sum(len(b) for b in bins.values())==
            d["expected_nonpal_source_assay_records"],
            "Original 21 assay genuine nonpal source SNP total changed")
    groups=[]
    for gid,g in target.items():
        records=bins[gid]
        require(len(records)==d["gene_variant_counts"][gid] and
                len(records)>=MIN_GENE_VARIANTS,"GTEx gene source SNP count mismatches original baseline")
        sorted_records=sorted(records.items())
        maf_se=[(v[0],v[1]) for _,v in sorted_records]
        N_cat=d["catalogue_n_hypothesis"]
        estimated=OLS.estimate_coloc_sdY(maf_se,N_cat)
        hist=collections.Counter(v[2] for _,v in sorted_records)
        require(min(hist)>0 and max(hist)==N_cat and len(hist)>1,
                "Expected real GTEx SNP genotype N variability absent")
        modal_only=[(v[0],v[1]) for _,v in sorted_records if v[2]==N_cat]
        require(len(modal_only)>=MIN_FULL_AN_VARIANTS,"Too few source SNPs with nominal AN/2")
        fullN_only=OLS.estimate_coloc_sdY(modal_only,N_cat)
        small=[v for v in maf_se if v[0]<.1]
        common=[v for v in maf_se if v[0]>=.1]
        low=(OLS.estimate_coloc_sdY(small,N_cat) if len(small)>=MIN_MAF_STRATUM else None)
        high=(OLS.estimate_coloc_sdY(common,N_cat) if len(common)>=MIN_MAF_STRATUM else None)
        odd=OLS.estimate_coloc_sdY(maf_se[::2],N_cat)
        even=OLS.estimate_coloc_sdY(maf_se[1::2],N_cat)
        base=estimated["sdY_est_coloc_zero_intercept_OLS"]
        minN=min(hist)
        minN_est=OLS.estimate_coloc_sdY(maf_se,minN)[
            "sdY_est_coloc_zero_intercept_OLS"]
        require(math.isclose(minN_est/base,math.sqrt(minN/N_cat),rel_tol=1e-12),
                "R coloc zero intercept OLS N scaling wrong")
        groups.append({
            "QTL_dataset":dataset,
            "gene_id":gid,"gene_symbol":g["gene_symbol"],
            "source_variant_class":"canonical_nonpalindromic_SNP_only",
            "source_original_cis_gene_assay_SNPs":len(records),
            "source_an_over_2_min":minN,
            "source_an_over_2_max":max(hist),
            "source_an_over_2_modal_study_N_records":hist[N_cat],
            "source_an_over_2_distinct_values":len(hist),
            "source_AN_per_SNP_heterogeneous":len(hist)>1,
            "source_an_over_2_hypothesis_only":N_cat,
            "published_QTL_model_gene_N_directly_verified":False,
            "source_direct_normalized_expression_sdY_observed":False,
            "original_coloc_sdY_est_at_nominal_study_N_UNAPPROVED":estimated,
            "sdY_est_at_min_observed_variant_AN_over_2_scenario_UNAPPROVED":minN_est,
            "sdY_est_nominal_N_on_only_full_AN_SNPs_UNAPPROVED":
                fullN_only["sdY_est_coloc_zero_intercept_OLS"],
            "all_vs_full_AN_variant_subset_relative_difference":
                abs(fullN_only["sdY_est_coloc_zero_intercept_OLS"]-base)/base,
            "source_MAF_below_0_1_rows":len(small),
            "source_MAF_at_least_0_1_rows":len(common),
            "MAF_low_sdY_est":low["sdY_est_coloc_zero_intercept_OLS"] if low else None,
            "MAF_high_sdY_est":high["sdY_est_coloc_zero_intercept_OLS"] if high else None,
            "MAF_low_high_ratio":(
                low["sdY_est_coloc_zero_intercept_OLS"]/
                high["sdY_est_coloc_zero_intercept_OLS"] if low and high else None),
            "sorted_odd_even_estimate_relative_difference":abs(
                odd["sdY_est_coloc_zero_intercept_OLS"]-
                even["sdY_est_coloc_zero_intercept_OLS"])/base,
            "sorted_odd_even_NOT_independent_due_to_LD":True,
            "source_expression_phenotype_sdY_approved_for_coloc":False,
            "source_variant_heterogeneous_N_not_accepted_as_model_N":True,
            "scientific_status":"HOLD_ORIGINAL_MODEL_N_SDY_AND_QTL_GWAS_CIS_HARMONIZATION",
            "formal_coloc_abf_eligible":False,
        })
    return {
        "source_region":"GRCh38_"+region["region"],
        "source_genuine_candidate_group_count":len(groups),
        "source_all_region_original_rows":scanned,
        "source_all_region_original_stream_SHA256":digest.hexdigest(),
        "source_gene_assay_nonpal_SNP_total":sum(len(v) for v in bins.values()),
        "source_original_allele_AN_missing_in_candidate_rows":source_an_missing,
        "source_variant_AN_heterogeneous_gene_assays":sum(
            x["source_AN_per_SNP_heterogeneous"] for x in groups),
        "source_dataset_nominal_sample_N_is_not_regression_N":N_cat,
        "n_QTL_estimates_only_UNAPPROVED":len(groups),
        "n_QTL_direct_phenotype_sdY_measured":0,
        "n_QTL_N_regression_approved":0,
        "n_formal_coloc_abf_approved":0,
        "raw_source_QTL_effects_or_variant_IDs_exported":False,
        "groups":groups,
    }


def audit(ancestry,manifest):
    require(ancestry in DATASETS,"Must be source-allowed GTEx QTL dataset")
    d=DATASETS[ancestry]
    source_region=original_region()
    print("SOURCE_ONE_"+ancestry+"_ORIGINAL_TBI_AND_BGZF_HTTP_RANGE",flush=True)
    lines,source_prov=ORIGINAL.source_indexed_range(d["url"],source_region)
    result=scan_source(lines,ancestry,manifest)
    return {
        "schema":"is-p3a-GTEx-21gene-original-QTL-coloc-sdY-est-AN-sensitivity-v1",
        "technical_status":"PASS_ORIGINAL_GTEx_QTL_SOURCE_AND_SDY_ESTIMATOR",
        "scientific_status":"HOLD_ORIGINAL_GENE_MODEL_N_SDY_CIS_MATCH",
        "QTL_dataset":ancestry,
        "original_QTL_accession":ancestry,
        "original_index_sha256":source_prov["tabix_index_SHA256"],
        "original_range_sha256":source_prov["raw_source_metadata"][
            "publisher_compressed_range_SHA256"],
        "original_source_read_once_only":True,
        "original_bulk_GWAS_or_11GB_QTL_file_downloaded":False,
        "original_source_source_record_AUDIT":result,
        "n_original_gene_assays":d["expected_assays"],
        "n_original_assay_groups_total":36,
        "n_formal_coloc_abf_executed":0,
        "n_formal_coloc_susie_executed":0,
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--dataset",choices=tuple(DATASETS),required=True)
    parser.add_argument("--manifest",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    result={
        "schema":"is-p3a-GTEx-21gene-original-QTL-coloc-sdY-est-AN-sensitivity-v1",
        "technical_status":"HOLD_GTEx_ORIGINAL_SOURCE_OR_ESTIMATOR",
        "scientific_status":"HOLD",
        "QTL_dataset":args.dataset,
        "n_formal_coloc_abf_executed":0,
    }
    try:
        manifest=json.loads(args.manifest.read_text(encoding="utf-8"))
        result=audit(args.dataset,manifest)
    except (OSError,ValueError,TypeError,KeyError,UnicodeError,TimeoutError,
            OverflowError,ArithmeticError) as exc:
        result["error_class_only"]=type(exc).__name__
        result["validation_error_reason"]=str(exc)[:220]  # static guarded QC exception, no SNP effect row
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",
                        encoding="utf-8")
    print(json.dumps({
        "dataset":args.dataset,
        "technical_status":result["technical_status"],
        "original_source_gene_assays":result.get("n_original_gene_assays"),
        "original_source_SNP_assay_records":result.get(
            "original_source_source_record_AUDIT",{}).get("source_gene_assay_nonpal_SNP_total"),
        "hypothetical_coloc_sdY_estimator_N_approved":False,
        "n_formal_coloc":0,
        "error_class_only":result.get("error_class_only"),
        "validation_error_reason":result.get("validation_error_reason")
    }),flush=True)
    return 0 if result["technical_status"].startswith("PASS_") else 1


if __name__=="__main__":
    raise SystemExit(main())
