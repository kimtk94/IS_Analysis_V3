#!/usr/bin/env python3
"""Merge two publisher-verifiable EUR/EAS original chr10 source QC reports.

36 candidate P3A QTL groups all belong to ONE chr10 region: no claim of
independent 36-locus GWAS replication, cohort sample N, SNP level harmonization,
signed matched LD, MR, or coloc.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

def check(p,s):
    if not p:raise ValueError(s)

def audit(eur,eas):
    pair={"EUR":eur,"EAS":eas}
    names=set()
    for a,x in pair.items():
        check(x.get("schema")=="is-p3a-chr10-original-publisher-gwas-all-36-v1"
              and x.get("technical_status")=="PASS_PUBLISHER_WHOLE_FILE_AND_CHR10_ORIGINAL_ROWS"
              and x.get("scientific_status")=="HOLD"
              and x.get("ancestry")==a
              and x.get("source",{}).get("publisher_md5_and_prior_sha256_verified") is True
              and x.get("context",{}).get("source_GWAS_effect_and_SE_original_rows_checked") is True
              and x.get("context",{}).get("source_GWAS_effects_reconciled_to_archived_fullcis_36_groups") is False
              and x.get("actual_36_assay_groups_considered")==36
              and x.get("unique_candidate_gene_ids")==20
              and x.get("new_coloc_performed")==0
              and x.get("independent_GWAS_genomic_loci_replicated")==0,
              "original EUR/EAS publisher source or scientific STOP not accepted: "+a)
        d=x["context"]
        check(d.get("source_context_GRCh37",{}).get("id")=="IS_XDATA_G0015"
              and d["source_context_GRCh37"].get("original_region_start")==104663676
              and d["source_context_GRCh37"].get("original_region_end")==105720351
              and d.get("source_context_valid_variant_rows",0)>100
              and d.get("original_region_valid_variant_rows",0)>0,
              "unexpected chr10 region source span")
        g=x.get("all_group_input_statuses",[])
        check(isinstance(g,list) and len(g)==36
              and len({(r.get("region_id"),r.get("QTL_dataset"),r.get("gene_id"))
                       for r in g})==36
              and len({r.get("gene_id") for r in g})==20
              and all(r.get("GWAS_source_ancestry")==a
                      and r.get("variant_level_QTL_GWAS_signed_effect_match_proven") is False
                      and r.get("GWAS_exact_case_fraction_and_group_N_approved") is False
                      and r.get("QTL_effect_SE_N_alleles_independently_approved") is False
                      and r.get("source_cohort_matched_LD_approved") is False
                      and r.get("new_coloc_or_replication_executed") is False for r in g),
              "GWAS technical PASS incorrectly promoted to group scientific eligibility")
        current={(r["region_id"],r["QTL_dataset"],r["gene_id"]) for r in g}
        if names:
            check(names==current,"EUR/EAS source reports do not cover the same 36 real groups")
        names=current
    lead=eas["context"].get("legacy_GWS_EAS_AIS_lead_comparison")
    return {
        "schema":"is-p3a-chr10-original-gwas-combined-36-qtl-groups-v1",
        "technical_status":"PASS_BOTH_PUBLISHER_ORIGINAL_GWAS_CHR10_CONTEXT_SCAN",
        "scientific_status":"HOLD",
        "publisher_accessions":{"EUR":eur["source"]["accession"],
                                "EAS":eas["source"]["accession"]},
        "EUR_chr10_context_valid_original_variant_rows":
            eur["context"]["source_context_valid_variant_rows"],
        "EAS_chr10_context_valid_original_variant_rows":
            eas["context"]["source_context_valid_variant_rows"],
        "EAS_legacy_regional_lead_comparison":lead,
        "real_group_keys":36,
        "real_gene_ids":20,
        "distinct_analysed_GWAS_region_count":1,
        "fully_approved_per_group_gwas_qtl_source_matches":0,
        "case_control_N_approved":False,
        "coloc_abf_executed":0,
        "coloc_susie_executed":0,
        "independent_stroke_replication_claimed":False,
        "all_groups_QTL_source_N_and_allele_still_HOLD":True,
        "required_next":"Original full cis QTL/GWAS matched signed effect+SE, cohort case fraction, donor N, official original source and matched LD",
        "per_group_status":[
            {"region_id":t[0],"QTL_dataset":t[1],"gene_id":t[2],
             "GWAS_technical_context_scanned_EUR_EAS":True,
             "QTL_GWAS_per_variant_same_source_effect_verified":False,
             "scientific_status":"HOLD_NO_PER_GROUP_ALLELE_N_AND_FULLCIS_STAT"}
            for t in sorted(names)
        ]
    }

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--eur",type=Path,required=True)
    a.add_argument("--eas",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    try:
        result=audit(json.loads(args.eur.read_text()),json.loads(args.eas.read_text()))
    except (ValueError,KeyError,TypeError,OSError) as e:
        result={"schema":"is-p3a-chr10-original-gwas-combined-36-qtl-groups-v1",
                "technical_status":"HOLD_QC_MISMATCH_OR_SOURCE_NOT_AVAILABLE",
                "scientific_status":"HOLD","coloc_abf_executed":0,
                "error":str(e)[:240]}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k in (
        "technical_status","scientific_status","real_group_keys","real_gene_ids",
        "EUR_chr10_context_valid_original_variant_rows",
        "EAS_chr10_context_valid_original_variant_rows","error")},ensure_ascii=False))
    return 0 if result["technical_status"].startswith("PASS_") else 1

if __name__=="__main__":
    raise SystemExit(main())
