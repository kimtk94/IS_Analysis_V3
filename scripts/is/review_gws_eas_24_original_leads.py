#!/usr/bin/env python3
"""Independent readback of 24 original EAS GWAS nonpal lead comparisons.

Only interprets results actually computed using the full original publisher
compressed GWAS. If effects, alleles, source digest, legacy IDs, case fraction
or 36-group coloc status drift, report HOLD instead of inventing an approval.
"""
from __future__ import annotations
import argparse
import collections
import importlib
import json
from pathlib import Path

M=importlib.import_module("scripts.is.audit_gws_eas_24_leads_original")
GATE=importlib.import_module("scripts.is.p3a_publisher_md5_gate")

def require(v,s):
    if not v:raise ValueError(s)

def inspect(source, old):
    manifest=M.read_legacy_180(old)
    exp={x["region_id"]:x for x in manifest["selected"]}
    require(source.get("schema")=="is-gws-eas-30-provisional-region-24-nonpal-original-variant-audit-v1"
            and source.get("technical_status")=="PASS_FULL_PUBLISHER_FILE_AND_24_LEAD_SOURCE_SCREEN"
            and source.get("scientific_status")=="HOLD"
            and source.get("historical_GWAS_table_blob_sha1")==M.EXPECTED_LEGACY_BLOB,
            "source schema/old SHA/source scientific status changed")
    pub=source.get("verified_whole_file_original_source",{})
    require(pub.get("accession")=="GCST90104545"
            and pub.get("compressed_original_file_bytes")==GATE.SOURCES["EAS"]["bytes"]
            and pub.get("official_publisher_md5_and_pinned_sha256") is True
            and source.get("original_GWAS_full_variant_rows_scanned")==6725216,
            "actual complete original GWAS publisher source cannot be trusted")
    rows=source.get("tested_leads",[])
    require(len(rows)==24
            and len({x.get("region_id") for x in rows})==24
            and {x.get("region_id") for x in rows}==set(exp),
            "not the exact 24 true original nonpal lead source identities")
    nreproduced=0;nallele=0;summary=[]
    for row in rows:
        expected=exp[row["region_id"]]
        require(row.get("GRCh37_variant_key")==expected["GRCh37_variant_key"]
                and row.get("nonpalindromic") is True,
                "variant key or nonpalindromic flag was changed")
        require(row.get("original_GWAS_SNP_effect_valid_for_P3A_coloc") is False
                if "original_GWAS_SNP_effect_valid_for_P3A_coloc" in row
                else True,
                "single GWAS lead not eligible for P3A coloc")
        accept=row.get("joint_effect_p_eaf_reproduced") is True
        if accept:
            require(row.get("status")=="SIGNED_ORIGINAL_EAS_LEAD_REPRODUCED_ARCHIVE"
                    and row.get("ALT_oriented_beta_match") is True
                    and row.get("p_value_log10_tolerance_match") is True
                    and row.get("ALT_oriented_eaf_match") is True
                    and row.get("allele_pair_matches")==1
                    and row.get("original_effect_se_finite_and_positive") is True,
                    "invalid original source joint agreement status")
            nreproduced+=1
        if row.get("allele_pair_matches")==1:nallele+=1
        summary.append({"region_id":row["region_id"],
                        "source_variant":expected["GRCh37_variant_key"],
                        "chromosome":expected["chromosome"],
                        "source_reproduction_status":row.get("status"),
                        "allele_pair_ok":row.get("allele_pair_matches")==1,
                        "signed_beta_p_eaf_reproduced":accept,
                        "scientific_status":"HOLD_NOT_INDEPENDENT_REPLICATION"})
    require(nreproduced==source.get("original_signed_beta_p_eaf_all_three_match_count")
            and nallele==source.get("original_same_position_allele_pair_match_count")
            and source.get("original_nonpalindromic_archived_leads_tested")==24,
            "posterior edited aggregate or per-SNP original technical result")
    chrom={x["chromosome"] for x in summary if x["signed_beta_p_eaf_reproduced"]}
    require(len(chrom)==source.get("unique_chromosomes_with_concordant_leads"),
            "multiple chromosome aggregation mismatch")
    pal=source.get("palindromic_leads_excluded")
    miss=source.get("historical_missing_leads")
    require({x["region_id"] for x in pal}==
            {x["region_id"] for x in manifest["withheld"]}
            and {x["region_id"] for x in miss}==
            {x["region_id"] for x in manifest["missing"]}
            and source.get("independently_excluded_palindromic_source_leads")==4
            and source.get("source_missing_lead_regions")==2,
            "archive palindromic/missing-source lead records drift")
    require(source.get("published_original_study_Ncase_Ncontrol_independently_approved") is False
            and source.get("BBJ_stroke_source_provenance_verified") is False
            and source.get("matched_QTL_cis_effects_verified") is False
            and source.get("coloc_abf_executed")==0
            and source.get("coloc_susie_executed")==0,
            "GWAS single marker technical check cannot grant scientific eligibility")
    tech=source.get("technical_source_archive_signed_reproduction_screen")
    require(tech==("PASS_REVIEW_REQUIRED" if nreproduced>=23
                   and len(chrom)>=10 and nallele>=23
                   else "HOLD_SOURCE_REPRODUCTION_INSUFFICIENT"),
            "source matching acceptance screen inconsistent with original rows")
    return {
        "schema":"is-gws-eas-original-24-variant-cross-region-readback-v1",
        "source_audit_status":"PASS_ORIGINAL_ROWS_RECONCILED_TO_PINNED_SOURCE_JSON",
        "source_match_screen":tech,
        "scientific_status":"HOLD",
        "original_source":pub["accession"],
        "historical_source_sha1":M.EXPECTED_LEGACY_BLOB,
        "original_ancestry":"EAS",
        "legacy_regions_reviewed":30,
        "nonpalindromic_leads_tested":24,
        "palindromic_leads_scientifically_held":4,
        "original_leads_missing":2,
        "signed_beta_p_eaf_reproduced":nreproduced,
        "signed_beta_p_eaf_discordant_or_missing":24-nreproduced,
        "allele_pair_matched_leads":nallele,
        "concordant_chromosome_count":len(chrom),
        "not_independent_BBJ_excluded_stroke_replication":True,
        "P3A_QTL_gene_assay_count":36,
        "P3A_unique_candidate_gene_count":20,
        "P3A_all_on_single_chromosome10_region":True,
        "P3A_exact_QTL_GWAS_source_rows_paired":0,
        "P3A_coloc_abf":0,"P3A_coloc_susie":0,
        "source_GWAS_case_fraction_N_approved":False,
        "signed_effect_variants":summary,
        "held_palindromic_leads":pal,"unavailable_source_leads":miss,
        "next_source_gate":"Archive FULL cis QTL/GWAS signed beta, SE, actual per-cell N and case fraction, plus source matched LD"
    }

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--from-actual-run",type=Path,required=True)
    a.add_argument("--legacy",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    try:
        x=inspect(json.loads(args.from_actual_run.read_text(encoding="utf-8")),
                  args.legacy)
    except (OSError,TypeError,ValueError,KeyError) as err:
        x={"schema":"is-gws-eas-original-24-variant-cross-region-readback-v1",
           "source_audit_status":"HOLD_SOURCE_ORIGINAL_READBACK_INVALID",
           "scientific_status":"HOLD","P3A_coloc_abf":0,"error":str(err)[:240]}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(x,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    for field in ("source_audit_status","source_match_screen","signed_beta_p_eaf_reproduced",
                  "signed_beta_p_eaf_discordant_or_missing","concordant_chromosome_count",
                  "P3A_coloc_abf","error"):
        if field in x:print(field,x[field])
    if "signed_effect_variants" in x:
        print("region_id\\tsource_variant\\tstatus")
        for row in sorted(x["signed_effect_variants"],key=lambda x:x["region_id"]):
            print(row["region_id"]+"\\t"+row["source_variant"]+"\\t"+
                  str(row["source_reproduction_status"]))
    return 0 if x["source_audit_status"].startswith("PASS_") else 1

if __name__=="__main__":
    raise SystemExit(main())
