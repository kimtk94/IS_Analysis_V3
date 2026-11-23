#!/usr/bin/env python3
"""Aggregate audit of real QTD000434 cis-QTLs for all 15 original ROSMAP genes.

Never upload individual genotypes, sample IDs, SNP effect rows or matrices.
One combined remote tabix query only; no retries (EMBL-EBI DoS rate guidance).
Historical Ensembl 105 phenotype_pos +/-1Mb is the *audit extraction range*,
NOT proof of the exact release7 qtlmap per-phenotype cis setting.
"""
from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import importlib
import io
import json
import math
import re
import subprocess
import urllib.parse
from pathlib import Path

PRIOR = importlib.import_module(
    "scripts.is.audit_p3a_official_coloc_ana_qtl_source_probe")
SMALL = importlib.import_module(
    "scripts.is.audit_p3a_36_original_zenodo_metadata_gwas_beta_model")
DATASET = "QTD000434"
STUDY = "QTS000025"
Z_FILE = "gene_counts_Ensembl_105_phenotype_metadata.tsv.gz"
Z_URL = "https://zenodo.org/records/7808390/files/" + Z_FILE + "?download=1"
Z_MD5 = "b999b21dcffeecbfe0bf886a5a91997b"
Z_MAX = 3_500_000
AUDIT_CIS_BP = 1_000_000
MAX_REGION_BP = 20_000_000
MAX_SOURCE_LINES = 3_500_000
MAX_GENE_UNIQUE_VARIANTS = 220_000
# Pinned exact column order was observed in original 2026-11-23 CI
# from an HTTP 206 gzip source prefix (SHA256 a138f18b65b3a8b1...)
COLUMNS = (
    "molecular_trait_id", "chromosome", "position", "ref", "alt", "variant",
    "ma_samples", "maf", "pvalue", "beta", "se", "type", "ac", "an", "r2",
    "molecular_trait_object_id", "gene_id", "median_tpm", "rsid",
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def stable_gene_id(gene_id):
    require(isinstance(gene_id, str), "Missing ENSG gene ID")
    stable = re.sub(r"\.\d+$", "", gene_id)
    require(re.fullmatch(r"ENSG\d{11}", stable) is not None,
            "Malformed original source ENSG stable ID")
    return stable


def source_gene_metadata(raw, gene_ids):
    require(isinstance(raw, bytes) and raw[:2] == b"\x1f\x8b",
            "Official Ensembl 105 phenotype metadata is not gzip")
    require(hashlib.md5(raw, usedforsecurity=False).hexdigest() == Z_MD5,
            "Publisher phenotype metadata MD5 did not match original Zenodo")
    reader = csv.DictReader(io.StringIO(gzip.decompress(raw).decode("utf-8-sig")),
                            delimiter="\t")
    needed={"phenotype_id", "gene_id", "phenotype_pos", "chromosome",
            "gene_start", "gene_end", "strand"}
    require(reader.fieldnames and needed.issubset(set(reader.fieldnames)),
            "Original phenotype metadata schema absent")
    result = {}
    for row in reader:
        value = stable_gene_id(row["gene_id"])
        if value not in gene_ids:
            continue
        phenotype_id = stable_gene_id(row["phenotype_id"])
        require(value == phenotype_id, "Gene expression phenotype ID differs from gene ID")
        require(value not in result, "Ambiguous 2023 gene phenotype metadata row")
        chrom = row["chromosome"].replace("chr", "")
        require(chrom == "10", "Original candidate phenotype not on chr10")
        pos = int(row["phenotype_pos"])
        start, end = int(row["gene_start"]), int(row["gene_end"])
        require(pos > 1_000_000 and 0 < min(start,end) <= pos <= max(start,end),
                "Original phenotype_pos outside annotated gene extent")
        result[value] = {"chromosome":"10", "phenotype_pos":pos,
                         "gene_interval_start":min(start,end),
                         "gene_interval_end":max(start,end),
                         "strand":row["strand"]}
    return result


def one_combined_region(meta, gene_ids):
    require(set(meta) == set(gene_ids), "Not all 15 gene phenotype positions resolved")
    positions=[x["phenotype_pos"] for x in meta.values()]
    lo=max(1,min(positions)-AUDIT_CIS_BP)
    hi=max(positions)+AUDIT_CIS_BP
    require(hi-lo+1 <= MAX_REGION_BP,
            "Candidate gene union region unexpectedly wide: manual review required")
    for meta_row in meta.values():
        require(lo <= meta_row["phenotype_pos"]-AUDIT_CIS_BP and
                hi >= meta_row["phenotype_pos"]+AUDIT_CIS_BP,
                "Union region does not fully contain a proposed phenotype cis window")
    return {"chromosome":"10", "start":lo,"end":hi,
            "region":f"10:{lo}-{hi}", "span_bp":hi-lo+1,
            "cis_half_window_assumed_from_official_qtlmap_default_bp":AUDIT_CIS_BP,
            "R7_actual_qtlmap_run_cis_parameter_independently_verified":False,
            "all_15_historical_phenotype_windows_contained":True}


def numeric_fields(row):
    try:
        an=int(row["an"])
        maf=float(row["maf"])
        beta=float(row["beta"])
        se=float(row["se"])
        p=float(row["pvalue"])
    except (KeyError, TypeError, ValueError, OverflowError) as exc:
        raise ValueError("Source variant contains invalid AN, MAF, beta, SE or p") from exc
    require(an>0 and an%2==0 and an<=2*100000 and
            math.isfinite(maf) and 0 < maf <= 0.5 and
            math.isfinite(beta) and math.isfinite(se) and se>0 and
            math.isfinite(p) and 0<=p<=1,
            "Source variant AN, MAF, effect or p out of bounds")
    return an,maf,beta,se,p


def source_variant_identity(row, allowed_positions):
    require(row["chromosome"] in {"10","chr10"},
            "Tabix returned non-chr10 variant in chr10 regional query")
    pos=int(row["position"])
    require(allowed_positions["start"]<=pos<=allowed_positions["end"],
            "Tabix output outside declared source extraction window")
    ref,alt=row["ref"].upper(),row["alt"].upper()
    require(ref and alt and re.fullmatch("[ACGT]+",ref) is not None and
            re.fullmatch("[ACGT]+",alt) is not None and ref!=alt,
            "Source REF/ALT invalid or ambiguous")
    expected=f"chr10_{pos}_{ref}_{alt}"
    require(row["variant"]==expected,
            "QTL variant ID does not equal GRCh38 position and REF/ALT")
    return expected,pos,alt


def aggregate_tabix_stream(lines, manifest, region, meta):
    gene_rows=[g for g in manifest["groups"] if g["QTL_dataset"]==DATASET]
    require(len(manifest["groups"])==36 and len(gene_rows)==15,
            "Wrong original ROSMAP manifest groups")
    targets={g["gene_id"]:g for g in gene_rows}
    require(set(meta)==set(targets), "Phenotype metadata != 15 ROSMAP group IDs")
    entries={g:{"raw_rows":0, "duplicate_rsid_rows":0,
                "per_variant":{}, "sites":collections.defaultdict(set)}
             for g in targets}
    scanned=0
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        parts=line.rstrip("\r\n").split("\t")
        if tuple(parts)==COLUMNS:  # acceptable if tabix -h includes header
            continue
        require(len(parts)==len(COLUMNS), "Tabix QTL source row has wrong column count")
        scanned+=1
        require(scanned<=MAX_SOURCE_LINES, "Original QTL streamed region exceeded hard row budget")
        row=dict(zip(COLUMNS,parts))
        gid=stable_gene_id(row["gene_id"])
        if gid not in targets:
            continue
        details=entries[gid]
        details["raw_rows"]+=1
        require(stable_gene_id(row["molecular_trait_id"])==gid,
                "Gene expression molecular_trait_id not identical to gene_id")
        variant,pos,alt=source_variant_identity(row,region)
        p0=meta[gid]["phenotype_pos"]
        require(abs(pos-p0)<=AUDIT_CIS_BP,
                "Original target variant outside candidate phenotype +/-1Mb window")
        an,maf,beta,se,p=numeric_fields(row)
        key=(an,maf,beta,se,p,row["ref"],row["alt"],row["type"])
        if variant in details["per_variant"]:
            require(details["per_variant"][variant]==key,
                    "Duplicated rsID source effect values inconsistent")
            details["duplicate_rsid_rows"]+=1
            continue
        details["per_variant"][variant]=key
        details["sites"][pos].add(alt)
        require(len(details["per_variant"]) <= MAX_GENE_UNIQUE_VARIANTS,
                "Unexpectedly large gene locus query or possible corrupted source")
    out=[]
    for gid,g in targets.items():
        info=entries[gid]
        data=list(info["per_variant"].values())
        an_values=[x[0] for x in data]
        maf_values=[x[1] for x in data]
        ctr=collections.Counter(an_values)
        n_unique=len(data)
        multi=sum(len(a)>1 for a in info["sites"].values())
        expected=g.get("official_full_cis_tested_variant_count")
        out.append({
            "gene_id":gid, "gene_symbol":g["gene_symbol"],
            "QTL_dataset":DATASET,
            "source_gene_2023_phenotype_pos":meta[gid]["phenotype_pos"],
            "source_rows_in_single_tabix_query":info["raw_rows"],
            "unique_full_allele_identity_variants_in_query":n_unique,
            "duplicated_rsid_variant_rows_consistent":info["duplicate_rsid_rows"],
            "multi_alt_positions_present":multi,
            "AN_min":min(an_values) if data else None,
            "AN_max":max(an_values) if data else None,
            "AN_mode":ctr.most_common(1)[0][0] if data else None,
            "AN_mode_fraction":(ctr.most_common(1)[0][1]/n_unique) if data else None,
            "AN_distinct_values":len(ctr),
            "MAF_min":min(maf_values) if data else None,
            "MAF_max":max(maf_values) if data else None,
            "n_rows_passing_numeric_effect_QC":n_unique,
            "archived_original_eGene_tested_variant_count_REFERENCE_ONLY":expected,
            "reference_count_matches_unique_variants":(expected==n_unique)
                if isinstance(expected,int) else None,
            "source_effect_rows_recovered":n_unique>0,
            "source_R7_run_specific_cis_window_certified":False,
            "source_genewise_regression_N_certified":False,
            "source_normalized_expression_sdY_measured":False,
            "scientific_status":"HOLD_R7_QTLMAP_REGRESSION_N_SDY_GWAS_CASE_CONTROL",
            "coloc_abf_eligible":False,
        })
    require(len(out)==15, "ROSMAP original target gene count drift")
    return {
        "query_scanned_all_study_QTL_rows_in_union_range":scanned,
        "n_ROSMAP_genes_with_original_any_SNP_QTL":sum(x["source_effect_rows_recovered"] for x in out),
        "n_ROSMAP_genes_with_missing_source_SNP_QTL":sum(not x["source_effect_rows_recovered"] for x in out),
        "sum_gene_raw_rows":sum(x["source_rows_in_single_tabix_query"] for x in out),
        "sum_gene_unique_full_variant_ids":sum(x["unique_full_allele_identity_variants_in_query"] for x in out),
        "sum_rsid_duplicates_removed_for_audit":sum(x["duplicated_rsid_variant_rows_consistent"] for x in out),
        "n_genes_all_AN_div2_560":sum(x["AN_distinct_values"]==1 and x["AN_mode"]==1120 for x in out),
        "n_genes_matching_archived_variant_count":sum(x["reference_count_matches_unique_variants"] is True for x in out),
        "official_R7_qtlmap_original_window_config_certified":False,
        "n_original_genes_approved_fullcis_for_formal_Coloc":0,
        "raw_SNP_ids_effects_or_donor_ids_exported":False,
        "groups":out,
    }


def run_tabix_single_request(url, region, manifest, meta):
    parsed=urllib.parse.urlsplit(url)
    require(parsed.scheme=="https" and parsed.hostname=="ftp.ebi.ac.uk" and
            parsed.path==("/pub/databases/spot/eQTL/sumstats/"
                          "QTS000025/QTD000434/QTD000434.all.tsv.gz"),
            "Original QTD000434 publisher URL changed")
    cmd=["tabix",url,region["region"]]
    # Explicit one query, no retries. Use bounded wall clock in CI timeout.
    # If server refuses, preserve HOLD; do not hammer EBI endpoints.
    print("TABIX_STAGE_COMMAND_START_ONE_REQUEST",flush=True)
    with subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                          text=True,encoding="utf-8",errors="strict",
                          bufsize=65536) as proc:
        try:
            output=aggregate_tabix_stream(proc.stdout,manifest,region,meta)
        except ValueError:
            print("TABIX_STAGE_STREAM_AGGREGATE_SCHEMA_OR_QC_ERROR",flush=True)
            proc.kill()
            proc.wait(timeout=15)
            raise
        except Exception:
            print("TABIX_STAGE_STREAM_SYSTEM_EXCEPTION",flush=True)
            proc.kill()
            proc.wait(timeout=15)
            raise
        status=proc.wait(timeout=60)
        stderr=proc.stderr.read(16384)
        if status!=0:
            lower=stderr.lower()
            flags={
                "missing_index":("index" in lower and ("fail" in lower or "not found" in lower)),
                "http_transport":("http" in lower or "curl" in lower or "remote" in lower),
                "file_open":("open" in lower or "failed" in lower),
                "tabix_no_results":output["query_scanned_all_study_QTL_rows_in_union_range"]==0,
            }
            print("TABIX_STAGE_NONZERO_EXIT_CLASSIFICATION="+json.dumps(
                {"nonzero_exit":True,"flags":flags},sort_keys=True),flush=True)
            raise ValueError("Original EBI tabix failed: diagnostic classification in CI log")
        print("TABIX_STAGE_SUCCESSFUL_EXIT_CODE_0",flush=True)
        require(bool(output["query_scanned_all_study_QTL_rows_in_union_range"]),
                "No original EBI tabix rows in union interval")
    output["tabix_result_returncode"]=status
    output["single_tabix_region_request_only"]=True
    output["automatic_server_retry_count"]=0
    output["tabix_source_original_full_file_downloaded"]=False
    output["tabix_stdout_saved_as_data_file"]=False
    output["source_stderr_recorded_in_artifact"]=False
    return output


def source_annotations():
    print("SOURCE_ORIGINAL_2023_ENSEMBL_105_PHENOTYPE_METADATA",flush=True)
    raw, meta=SMALL.get(Z_URL,{"zenodo.org","www.zenodo.org"},Z_MAX)
    require(hashlib.md5(raw,usedforsecurity=False).hexdigest()==Z_MD5,
            "Publisher original phenotype metadata deposited MD5 changed")
    return raw,{
        "original_publisher_record":"https://zenodo.org/records/7808390",
        "original_phenotype_annotation_filename":Z_FILE,
        "original_publisher_deposited_MD5":Z_MD5,
        "original_metadata_bytes":len(raw),
        "original_metadata_SHA256":meta["SHA256"],
        "2023_en_105_phenotype_positions_are_NOT_release7_run_config":True,
    }


def execute(manifest):
    raw,source_prov=source_annotations()
    genes={g["gene_id"] for g in manifest["groups"] if g["QTL_dataset"]==DATASET}
    require(len(genes)==15, "ROSMAP target 15 genes manifest mismatch")
    phenotype=source_gene_metadata(raw,genes)
    require(len(phenotype)==15, "Missing one or more official 2023 Ensembl105 genes")
    region=one_combined_region(phenotype,genes)
    print("SOURCE_ONE_TABIX_REQUEST_FOR_ALL_15_ROSMAP_GENES_"+region["region"],
          flush=True)
    url=("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"
         "QTS000025/QTD000434/QTD000434.all.tsv.gz")
    stats=run_tabix_single_request(url,region,manifest,phenotype)
    return {
        "schema":"is-p3a-36-original-ROSMAP-15-gene-single-tabix-regional-QTL-aggregates-v1",
        "technical_status":"PASS_ORIGINAL_SINGLE_REGION_TABIX_15_GENE_AGGREGATE_SOURCE",
        "scientific_status":"HOLD_RELEASE7_CIS_PROVENANCE_QTL_N_SDY_AND_GWAS_CC",
        "source_provenance":source_prov,
        "original_summary_source_url":url,
        "original_summary_source_snp_rows_saved":False,
        "audit_query_region":region,
        "QTL_15_gene_aggregate":stats,
        "n_all_assays":36,"n_ROSMAP_assays":15,"n_distinct_genes":20,
        "n_formal_coloc_abf_eligible":0,"n_formal_coloc_abf_executed":0,
        "n_formal_coloc_susie_eligible":0,"n_formal_coloc_susie_executed":0,
        "original_21_GTEx_assays_status":"HOLD_NOT_REAUDITED_IN_THIS_WORKFLOW",
        "original_GIGASTROKE_GWAS_status":"HOLD_BINARY_LOG_OR_CASE_FRACTION_NOT_APPROVED",
        "all_source_private_donor_ids_absent_from_artifact":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    data={
        "schema":"is-p3a-36-original-ROSMAP-15-gene-single-tabix-regional-QTL-aggregates-v1",
        "technical_status":"HOLD_ORIGINAL_ENSEMBL105_OR_SINGLE_EBI_TABIX",
        "scientific_status":"HOLD",
        "n_formal_coloc_abf_eligible":0,"n_formal_coloc_abf_executed":0,
    }
    try:
        manifest=json.loads(args.manifest.read_text(encoding="utf-8"))
        data=execute(manifest)
    except (ValueError,OSError,TypeError,KeyError,UnicodeError,TimeoutError,
            subprocess.TimeoutExpired,gzip.BadGzipFile,EOFError) as err:
        data["error_class_only"]=type(err).__name__
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(data,ensure_ascii=False,indent=2)+"\n",
                        encoding="utf-8")
    print(json.dumps({
        "technical_status":data["technical_status"],
        "scientific_status":data["scientific_status"],
        "n_groups_with_any_original_QTL":data.get("QTL_15_gene_aggregate",{}).get(
            "n_ROSMAP_genes_with_original_any_SNP_QTL"),
        "query_region":data.get("audit_query_region",{}).get("region"),
        "n_formal_coloc_abf_executed":0,
        "error_class_only":data.get("error_class_only"),
    },ensure_ascii=False),flush=True)
    return 0 if data["technical_status"].startswith("PASS_") else 1


if __name__=="__main__":
    raise SystemExit(main())
