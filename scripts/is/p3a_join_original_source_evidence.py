#!/usr/bin/env python3
"""Review published GWAS byte verification + original GWAS/QTL CI artifacts
against old 36-group P3A aggregate status.

This is a scientific acceptance GATE, not a coloc executor, and not a source
of invented 36 group identities. It fails closed on any missing old evidence,
and reports missing group-level JSON explicitly.
"""
from __future__ import annotations
import argparse
import importlib
import json
from pathlib import Path

class EvidenceError(ValueError): pass

def require(ok,reason):
    if not ok:raise EvidenceError(reason)

def parse(path):
    raw=json.loads(Path(path).read_text(encoding="utf-8"))
    require(isinstance(raw,dict),"source artifact not a JSON object")
    return raw

def inspect(ancestry, verified, old):
    md5gate=importlib.import_module("scripts.is.p3a_publisher_md5_gate")
    md5gate.validate_result(verified,ancestry)
    require(verified.get("technical_status")=="PASS_PUBLISHER_MD5_SHA256_SIZE",
            "publisher digest not passed")
    hits=old.get("results")
    require(isinstance(hits,list) and len(hits)==1,"prior full-source report missing")
    row=hits[0]
    ref=md5gate.SOURCES[ancestry]
    require(row.get("ancestry")==ancestry
            and row.get("accession")==ref["accession"]
            and row.get("full_source_status")=="TRANSFER_AND_STRUCTURE_QC_PASS"
            and row.get("total_variant_rows")==ref["original_rows_in_prior_source_qc"]
            and row.get("bytes_downloaded")==ref["bytes"]
            and row.get("sha256_calculated_local_not_publisher_verified")==ref["full_sha256_pinned"]
            and row.get("per_variant_sample_size_column_present") is False
            and row.get("source_independent_sha256_verified") is False
            and old.get("research_ready") is False,
            "prior full-GWAS source identity/coverage drift")
    return {
        "accession":ref["accession"],
        "publisher_original_full_archive_digest_verified":True,
        "original_rows_from_prior_full_QC":ref["original_rows_in_prior_source_qc"],
        "source_allele_effect_compared_to_archived_BBJ_EAS5":False,
        "source_outcome_case_control_N_validated_for_archived_stage":False,
        "science_status":"SOURCE_FILE_INTEGRITY_PASS_COLOC_HOLD"
    }

def summarize(baseline, eur_now, eas_now, eur_prior, eas_prior, eqtl, supplement,
              optional_groups=None):
    require(baseline.get("schema")=="is-candidate-wide-coloc-baseline-v1"
            and baseline.get("scientific_status")=="HOLD", "bad 2225/393 source baseline")
    g=baseline.get("assay_group_universe",{})
    require(g.get("readiness_region_gene_tissue_groups")==393
            and g.get("abf_followup_source_flagged_groups")==36
            and g.get("abf_followup_groups_with_some_valid_source_stats")==25
            and g.get("abf_followup_groups_without_valid_source_stats")==11
            and g.get("coloc_abf_eligible_groups")==0
            and g.get("coloc_abf_executed_new_stage_groups")==0,
            "wrong legacy P3A counts or false scientific approval")

    euro=inspect("EUR",eur_now,eur_prior)
    asia=inspect("EAS",eas_now,eas_prior)
    require(eqtl.get("scientific_decision")=="SOURCE_ROW_QC_PASS_COLOC_HOLD"
            and eqtl.get("coloc_ABF_beta_varbeta_route")=="BLOCKED_MISSING_VARBETA"
            and eqtl.get("coloc_ABF_pvalue_MAF_N_route")=="BLOCKED_MISSING_MAF_AND_CELLTYPE_N"
            and eqtl.get("actual_per_celltype_N")=="NOT_ESTABLISHED_BY_416_COHORT_ROWS"
            and eqtl.get("source_nominal_build")=="GRCh38",
            "original Japanese cis eQTL source effects/N approved incorrectly")
    require(supplement.get("schema")=="is-ci-public-source-gate-v1"
            and supplement.get("gigastroke_status")=="PASS"
            and supplement.get("bbj_2021",{}).get("status")=="DISCOVERED"
            and supplement.get("bbj_2021",{}).get("file_downloaded") is False
            and supplement.get("bbj_2021",{}).get("cohort_independent") is False
            and supplement.get("independent_replication")=="NOT_ESTABLISHED",
            "original BBJ discovery/supplement not proof of independent replication")

    group_keys=None
    if optional_groups is not None:
        module=importlib.import_module("scripts.is.candidate_wide_coloc_source_triage")
        by_group=module.validate_per_group(optional_groups,g)
        group_keys=[{"region_id":x["region_id"],"QTL_dataset":x["QTL_dataset"],
                     "gene_id":x["gene_id"],"source_stat_screen":x["source_stat_screen"],
                     "blocked_on":x["blockers"],"eligible":False} for x in by_group]
        require(len(group_keys)==36,"actual 36 group report not complete")
    return {
        "schema":"is-p3a-all-source-joined-readiness-v1",
        "technical_status":"PASS_ORIGINAL_GWAS_PUBLISHER_DIGEST_AND_SOURCE_EVIDENCE",
        "scientific_status":"HOLD",
        "publisher_GWAS":{
            "EUR":euro,"EAS":asia
        },
        "Japanese_eQTL":"SOURCE_QC_PASS_ABF_BLOCKED_MISSING_EFFECT_SE_N",
        "candidate_gene_universe":2225,
        "source_readiness_assay_groups":393,
        "original_source_ABF_flagged_groups":36,
        "original_source_ABF_groups_basic_stats":25,
        "original_source_ABF_groups_no_stats":11,
        "individual_original_group_IDs_loaded":group_keys is not None,
        "per_group_audit_rows":len(group_keys or []),
        "unresolved_original_group_ID_count":0 if group_keys is not None else 36,
        "genuinely_coloc_eligible_groups":0,
        "new_coloc_abf_computations":0,
        "independent_BBJ_excluded_stroke_replications":0,
        "source_effect_allele_harmonization_approved":False,
        "source_case_fraction_effective_N_approved":False,
        "publisher_md5_alone_not_sufficient_for_coloc":True,
        "per_group_blockers_if_available":group_keys,
        "required_next_inputs":[
          "36 ORIGINAL assay group identities from server archived_fasta_abf_audit.json",
          "Archived EAS5 full cis QTL–GWAS row dataset and original official GWAS effects",
          "Signed harmonized allele+build, exact phenotype case fraction/N, study QTL N and SNP overlap",
          "Predeclared priors plus limited stable coloc and multi-signal SuSiE decisions"
        ],
    }

def main():
    p=argparse.ArgumentParser()
    for field in ["baseline","eur_now","eas_now","eur_prior","eas_prior",
                  "eqtl","supplement","out"]:
        p.add_argument("--"+field.replace("_","-"),type=Path,required=True)
    p.add_argument("--original-groups",type=Path)
    a=p.parse_args()
    try:
        inputs=[parse(getattr(a,k)) for k in
                ("baseline","eur_now","eas_now","eur_prior","eas_prior","eqtl","supplement")]
        group=parse(a.original_groups) if a.original_groups else None
        status=summarize(*inputs,optional_groups=group)
    except (OSError,ValueError,TypeError,KeyError) as e:
        status={
            "schema":"is-p3a-all-source-joined-readiness-v1",
            "technical_status":"HOLD_SOURCE_VALIDATION_FAILURE",
            "scientific_status":"HOLD","genuinely_coloc_eligible_groups":0,
            "new_coloc_abf_computations":0,"error":str(e)[:300]
        }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(status,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in status.items() if k in (
       "technical_status","scientific_status","unresolved_original_group_ID_count",
       "per_group_audit_rows","genuinely_coloc_eligible_groups",
       "new_coloc_abf_computations","error")},ensure_ascii=False))
    return 0 if status["technical_status"].startswith("PASS_") else 1

if __name__=="__main__":
    raise SystemExit(main())
