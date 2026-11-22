#!/usr/bin/env python3
"""Fast retrospective reconciliation: actual EAS GIGASTROKE chr10 original source
to the preserved 2026-10-09 30-region GWAS descriptive table.

Checks source-origin SHA1 and counts / significance / minimum-p / a single
C/G PALINDROMIC regional lead, without confusing count matches with
independent 20-SNP source identity or approved coloc / GWAS cohort N.
No raw per-variant results are uploaded by this script.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib
import io
import json
import math
from pathlib import Path

JOIN=importlib.import_module("scripts.is.join_p3a_chr10_publisher_source")
SOURCE_FILE="IS_30_REGIONS_SIX_OUTCOMES_GWAS_CONTEXT_500KB.tsv"
SOURCE_SHA="b34b4ae84ea4e5b0c72e3982d530a048cc9be5ca"

def require(x,reason):
    if not x:raise ValueError(reason)

def source_table(path):
    raw=path.read_bytes()
    digest=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    require(digest==SOURCE_SHA,"pinned historical 30-region source Git blob drift")
    with io.StringIO(raw.decode("utf-8-sig")) as f:
        rdr=csv.DictReader(f,delimiter="\t")
        rows=list(rdr)
    require(len(rows)==180
            and len({(x["region_id"],x["source_dataset"]) for x in rows})==180,
            "expected all 30 regions and six original outcome sources")
    matches=[r for r in rows if r["region_id"]=="IS_XDATA_G0015"
             and r["source_dataset"]=="GWS_EAS_AIS"]
    require(len(matches)==1,"missing exact archived EAS AIS source row for chr10")
    x=matches[0]
    require(x["chromosome_GRCh37"]=="10"
            and int(x["original_region_start"])==104663676
            and int(x["original_region_end"])==105720351
            and int(x["context_window_start"])==104163676
            and int(x["context_window_end"])==106220351
            and x["group_source_lead_variant"]=="10:105596641:G:C",
            "historical GRCh37 scope or C/G regional lead drift")
    return x

def reconcile(eur,eas,old):
    new=JOIN.audit(eur,eas)
    c=eas["context"]
    source_rows=int(old["n_rows_in_context"])
    source_sig=int(old["n_genomewide_p_lt_5e8"])
    source_minp=float(old["best_p"])
    rcount=c["source_context_valid_variant_rows"]
    sig=c["source_context_p_lt_5e8_rows"]
    minp=c["source_context_min_p_descriptive"]
    count_ok=rcount==source_rows
    sig_ok=sig==source_sig
    minp_ok=minp>0 and source_minp>0 and abs(math.log10(minp)-math.log10(source_minp))<=.05
    lead=c.get("legacy_GWS_EAS_AIS_lead_comparison") or {}
    require(lead.get("source_identity_approved") is False if lead else True,
            "palindromic lead cannot authenticate official GWAS")
    lead_cmp=lead.get("status","NOT_PRESENT_OR_UNAVAILABLE")
    lead_pass=lead_cmp=="PALINDROMIC_LEAD_NUMERIC_CONCORDANCE_STRAND_HOLD"
    require(lead.get("n_independent_original_variant_matches_for_source_approval",0)==0,
            "palindromic single lead falsely counted for source identity")
    return {
        "schema":"is-p3a-chr10-old-new-actual-GWAS-descriptive-match-v1",
        "technical_status":"PASS_GIT_AND_FULL_PUBLISHER_SOURCE_AUDIT",
        "scientific_status":"HOLD",
        "original_source_commit":"0bf646fd6f29322d17da878253f37def24c01c74",
        "pinned_30_region_Git_blob_sha1":SOURCE_SHA,
        "official_EAS_GWAS_accession":eas["source"]["accession"],
        "legacy_region":"IS_XDATA_G0015",
        "original_GWAS_EAS_context_valid_rows":rcount,
        "archived_GWS_EAS_AIS_context_valid_rows":source_rows,
        "regional_original_row_count_equal":count_ok,
        "original_GWAS_EAS_context_p_lt_5e8_rows":sig,
        "archived_GWS_EAS_AIS_context_p_lt_5e8_rows":source_sig,
        "regional_genomewide_significance_count_equal":sig_ok,
        "original_vs_archived_best_p_log10_within_0p05":minp_ok,
        "palindromic_single_lead_numeric_comparison":lead_cmp,
        "palindromic_single_lead_numerically_concordant":lead_pass,
        "descriptive_region_signature_match":count_ok and sig_ok and minp_ok,
        "twenty_independent_nonpal_source_SNPs_matched":False,
        "sample_size_and_BBJ_study_source_approved":False,
        "matched_fullcis_QTL_GWAS_effects_verified":False,
        "new_P3A_coloc_executed":0,
        "groups":36,"genes":20,
        "interpretation":"Matching aggregate counts/p-values across one region support source consistency but NOT per-variant signed allele provenance, multiple independent loci, cohort N or functional causality",
    }

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--eur",type=Path,required=True)
    a.add_argument("--eas",type=Path,required=True)
    a.add_argument("--legacy30",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    try:
        d=reconcile(json.loads(args.eur.read_text()),json.loads(args.eas.read_text()),
                    source_table(args.legacy30))
    except (ValueError,OSError,KeyError,TypeError) as e:
        d={"schema":"is-p3a-chr10-old-new-actual-GWAS-descriptive-match-v1",
           "technical_status":"HOLD_SOURCE_OR_REGION_PROVENANCE_ERROR",
           "scientific_status":"HOLD","new_P3A_coloc_executed":0,"error":str(e)[:240]}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in d.items() if k in (
        "technical_status","scientific_status","descriptive_region_signature_match",
        "regional_original_row_count_equal","regional_genomewide_significance_count_equal",
        "original_vs_archived_best_p_log10_within_0p05",
        "palindromic_single_lead_numeric_comparison","error")},ensure_ascii=False))
    return 0 if d["technical_status"].startswith("PASS_") else 1

if __name__=="__main__":
    raise SystemExit(main())
