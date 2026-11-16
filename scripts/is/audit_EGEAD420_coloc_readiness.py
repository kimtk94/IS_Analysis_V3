#!/usr/bin/env python3
"""Independent, FAIL-CLOSED E-GEAD-420 IS eQTL coloc readiness audit.

Input: genuine GitHub Actions source artifacts from separate original
ImmuNexUT 13-column association QC and independently MD5-verified 416-record
SDRF. No original 43 GB downloaded. NO SNP sharing, liftover or coloc run.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import os

QTL_ID="E-GEAD-420"
SOURCE_RUN=43847910968
COHORT_RUN=43433142595
SOURCE_STATUS="OFFICIAL_PUBLIC_43GB_ZIP64_CENTRAL_AND_NOMINAL_TAR_PREFIX_VERIFIED"
COHORT_STATUS="SOURCE_SDRF_416_COHORT_METADATA_AND_PUBLISHER_MD5_PASS"


def audit_inputs(raw,metadata):
    if raw.get("source_status")!=SOURCE_STATUS:
        raise ValueError("Actual EAS source archive structural QC is not verified")
    if raw.get("published_whole_ZIP_MD5_verified_by_full_download") is not False:
        raise ValueError("43GB source full-checksum status unexpectedly changed")
    if raw.get("nominal_source_build")!="GRCh38" or raw.get("source_GWAS_build")!="GIGASTROKE GRCh37":
        raise ValueError("Source genome build identity mismatch")
    first=raw.get("first_actual_nominal_association_records",{})
    audit=first.get("actual_ImmuNexUT_13_column_source_data_QC")
    if not isinstance(audit,dict) or audit.get("source_13col_header_verified") is not True:
        raise ValueError("Missing genuine original ImmuNexUT 13-column source evidence")
    if first.get("source_nominal_file")!="CD16p_Mono_nominal.txt":
        raise ValueError("Original immune-cell context not proven")
    if first.get("whitespace_delimited_columns_per_line")!=13:
        raise ValueError("Unexpected original source col count")
    n=audit.get("n_actual_complete_source_association_rows")
    nonsig=audit.get("source_nominal_non_significant_p_gt_0_05")
    sig=audit.get("source_nominal_significant_p_le_0_05")
    if not isinstance(n,int) or n<5 or not all(isinstance(x,int) for x in (nonsig,sig)):
        raise ValueError("Actual P-value summary unverified")
    if nonsig+sig!=n or nonsig==0:
        raise ValueError("Original nominal sample doesn't have verified nonsignificant pairs")
    if (audit.get("slope_beta_present") is not True or
            audit.get("beta_effect_allele")!="ALT" or
            audit.get("ref_alt_alleles_present") is not True):
        raise ValueError("Allele effect fields absent")
    if (audit.get("slope_se_column_present") is not False or
            audit.get("allele_frequency_MAF_EAF_present") is not False or
            audit.get("per_celltype_sample_n_column_present") is not False):
        raise ValueError("Original source-field status changed; re-evaluate readiness")
    if audit.get("full_cis_association_completeness")!="NOT_VERIFIED":
        raise ValueError("Full cis coverage unexpectedly assumed")
    if metadata.get("status")!=COHORT_STATUS or metadata.get("original_SDRF_publisher_MD5_verified") is not True:
        raise ValueError("Genuine cohort public SDRF original checksum not verified")
    cohort=metadata.get("aggregate_evidence",{})
    if (cohort.get("original_public_SDRF_distinct_source_names")!=416 or
            cohort.get("healthy_control_WGS_records")!=79 or
            cohort.get("immune_disease_WGS_records")!=337 or
            cohort.get("per_celltype_effective_eQTL_N_observed") is not False):
        raise ValueError("SDRF cohort counts incorrect or per-cell n overstated")
    return {
        "source":"ImmuNexUT E-GEAD-420",
        "verified_source_run":SOURCE_RUN,
        "verified_cohort_run":COHORT_RUN,
        "original_nominal_association_rows_sampled":n,
        "original_non_significant_p_gt_0_05":nonsig,
        "original_significant_p_le_0_05":sig,
        "first_cell_type":"CD16p_Mono",
        "first_gene":audit["source_first_gene_id"],
        "source_reported_first_gene_cis_variant_count":audit["source_reported_cis_variants_for_first_gene"],
        "original_beta_slope_ALT":"PASS",
        "original_per_variant_standard_error":"MISSING_FROM_HEADER",
        "original_allele_frequency":"MISSING_FROM_HEADER",
        "actual_per_celltype_N":"NOT_ESTABLISHED_BY_416_COHORT_ROWS",
        "independent_source_IDF_SDRF_checksums":"PASS_SMALL_FILES_ONLY",
        "whole_source_43GB_checksum":"NOT_VERIFIED",
        "source_nominal_build":"GRCh38",
        "GIGASTROKE_is_GWAS_build":"GRCh37",
        "build_liftover_and_REF_ALT_normalization":"NOT_RUN",
        "complete_first_gene_cis_variants":"PENDING_INDEPENDENT_BLOCK_TEST",
        "all_stroke_locus_variant_coverage":"NOT_TESTED",
        "coloc_ABF_beta_varbeta_route":"BLOCKED_MISSING_VARBETA",
        "coloc_ABF_pvalue_MAF_N_route":"BLOCKED_MISSING_MAF_AND_CELLTYPE_N",
        "SuSiE_coloc_and_finemap":"BLOCKED_MISSING_HARMONIZED_LD_N",
        "clinical_stroke_genotype_expression_measured":False,
        "scientific_decision":"SOURCE_ROW_QC_PASS_COLOC_HOLD",
    }


def run(root,out):
    source=root/"source"/"E_GEAD_420_source_zip64_evidence.json"
    cohort=root/"cohort"/"E_GEAD_420_cohort_metadata.json"
    if not source.is_file() or not cohort.is_file():
        raise ValueError("Source or cohort GitHub evidence ZIP not downloaded")
    status=audit_inputs(json.loads(source.read_text(encoding="utf-8")),
                        json.loads(cohort.read_text(encoding="utf-8")))
    status["assembled_at_utc"]=datetime.now(timezone.utc).isoformat()
    out.mkdir(parents=True,exist_ok=True)
    (out/"EGEAD420_reproducible_go_hold.json").write_text(
        json.dumps(status,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    summary=f"""# E-GEAD-420 actual verified source evidence — coloc feasibility

**{status['scientific_decision']}** — actual source rows independently cross-checked against official SDRF MD5.

- True CD16+ monocyte eQTL source rows inspected: {status['original_nominal_association_rows_sampled']} (non-significant {status['original_non_significant_p_gt_0_05']}).
- beta/slope(ALT) + variant GRCh38 REF/ALT exist; **SE, MAF, per-cell type N absent**.
- Total study 416 (healthy 79, immune-disease 337) does **not** prove CD16+ monocyte effective eQTL sample N.
- Coloc ABF beta/varbeta route **BLOCKED**; p/MAF/N route **BLOCKED**.
- GRCh38 to GRCh37 liftOver/REF-ALT reconciliation, per-gene whole-cis coverage, ancestry-matched LD and power sensitivity not done.
- No ischemic stroke mechanism, joint association or colocalization calculated.
"""
    (out/"summary.md").write_text(summary,encoding="utf-8")
    print(summary)
    if name:=os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(name,"a",encoding="utf-8") as f:
            f.write(summary)
    return status


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",default="original_evidence")
    p.add_argument("--out-dir",default="artifacts/is_EGEAD420_coloc_readiness")
    a=p.parse_args()
    run(Path(a.source),Path(a.out_dir))
