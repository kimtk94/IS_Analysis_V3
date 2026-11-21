#!/usr/bin/env python3
"""Source-backed per-chromosome SNP attrition QC; no LDSC or summary-stat export."""
from __future__ import annotations
import argparse
import json
from pathlib import Path

from scripts.ckm_stroke.review_gigastroke_full_evidence import (
    EXPECTED_BMI_CANDIDATE, PINNED, review,
)

COLUMNS=(
    "candidate_reference_bmi_snps",
    "unique_outcome_coordinate_snps",
    "outcome_duplicate_coordinate_exclusions",
    "candidate_not_in_outcome",
    "allele_aligned_numeric_valid_pre_munge",
    "nonaligned_or_numeric_invalid",
)

def checked_num(row, field):
    x=row.get(field)
    if type(x) is not int or x < 0:
        raise ValueError("invalid per-chromosome SNP count: "+field)
    return x

def audit_attrition(source:dict):
    """Requires full original publisher MD5 + prior independent global counts."""
    validated=review(source)
    source_outcomes={o["phenotype"]:o for o in source["outcomes"]}
    reference_profile=None
    summaries=[]
    for phenotype in ("AIS","LAA","CE","SVD"):
        raw=source_outcomes[phenotype]
        original_counts=raw["counts"]
        per_chr=raw.get("chr_attrition_pre_munge")
        if not isinstance(per_chr,list) or len(per_chr)!=22:
            raise ValueError("missing 22 per-chromosome rows: "+phenotype)
        chroms=[r.get("chromosome") for r in per_chr if isinstance(r,dict)]
        if chroms!=list(range(1,23)):
            raise ValueError("chromosome order, identity or number invalid: "+phenotype)
        agg={key:0 for key in COLUMNS}
        profile=[]
        summary_chr=[]
        for row in per_chr:
            counts={k:checked_num(row,k) for k in COLUMNS}
            c=counts["candidate_reference_bmi_snps"]
            u=counts["unique_outcome_coordinate_snps"]
            d=counts["outcome_duplicate_coordinate_exclusions"]
            m=counts["candidate_not_in_outcome"]
            good=counts["allele_aligned_numeric_valid_pre_munge"]
            invalid=counts["nonaligned_or_numeric_invalid"]
            if c==0 or c!=u+d+m or u!=good+invalid:
                raise ValueError("chromosome SNP partition does not reconcile: "+phenotype)
            for k,v in counts.items():agg[k]+=v
            profile.append(c)
            summary_chr.append({"chr":row["chromosome"],**counts})
        if agg["candidate_reference_bmi_snps"]!=EXPECTED_BMI_CANDIDATE:
            raise ValueError("candidate SNP reference total does not match pinned source")
        if reference_profile is None:
            reference_profile=profile
        elif profile!=reference_profile:
            raise ValueError("inconsistent reference chromosome SNP universe among outcomes")
        expected={
            "unique_outcome_coordinate_snps":"unique_outcome_coordinate_matches",
            "outcome_duplicate_coordinate_exclusions":"outcome_collision_coordinates_excluded",
            "candidate_not_in_outcome":"candidate_coordinates_not_in_outcome",
            "allele_aligned_numeric_valid_pre_munge":"unique_harmonizable_BMI_positions",
        }
        if any(agg[k]!=original_counts.get(v) for k,v in expected.items()):
            raise ValueError("per-chromosome counts do not match original source: "+phenotype)
        if (agg["nonaligned_or_numeric_invalid"]+
            agg["allele_aligned_numeric_valid_pre_munge"]!=
            agg["unique_outcome_coordinate_snps"]):
            raise ValueError("not partitioned by valid allele/numeric evidence")
        if (raw.get("accession")!=PINNED[phenotype][0]
            or raw.get("outcome_per_variant_N_available") is not False
            or raw.get("raw_md5_verified") is not True):
            raise ValueError("outcome source accession, MD5 or N availability changed")
        summaries.append({
            "phenotype":phenotype,
            "accession":PINNED[phenotype][0],
            "source_publisher_MD5_verified":True,
            "chromosomes":summary_chr,
            "totals":agg,
            "pre_munge_reference_aligned_coverage_pct":
                round(100*agg["allele_aligned_numeric_valid_pre_munge"]/EXPECTED_BMI_CANDIDATE,5),
            "below_200k_prefilter_overlap":agg["allele_aligned_numeric_valid_pre_munge"]<200000,
            "final_LDSC_regression_SNP_count_verified":False,
            "scientific_status":"HOLD",
        })
    return {
        "schema":"ckm-stroke-eur-bmi-2018-GIGASTROKE-chr-attrition-pre-munge-v1",
        "technical_status":"PASS_ORIGINAL_MD5_PER_CHROMOSOME_PRE_MUNGE_COUNTS_ONLY",
        "ancestry":"EUR","genome_build":"GRCh37",
        "source_BMI_candidate_count":EXPECTED_BMI_CANDIDATE,
        "reference_source_note":"BMI candidate originally selected from pinned triple HM3 and 22 canonical EUR LD score rsIDs",
        "outcomes":summaries,
        "not_actual_ldsc_regression_counts":True,
        "scientific_status":"HOLD",
        "ready_for_ldsc":False,
        "approved_case_control_N_eff":False,
        "LDSC_h2_intercept_rg_MR_coloc_performed":False,
        "remaining_gates":[
            "Outcome case-control meta-analysis SNP-specific effective N unapproved",
            "SNP-specific munge filters and full Z/SE validation not completed",
            "Independent weight-LD/LD-score regression SNP merge and post-filter count not computed",
            "Genomic-control preprocessing provenance and participant overlap unresolved"
        ],
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--full-source-evidence",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    try:
        r=audit_attrition(json.loads(args.full_source_evidence.read_text("utf-8")))
    except Exception as exc:
        r={"technical_status":"HOLD_CHR_ATTRITION_EVIDENCE_INVALID",
           "scientific_status":"HOLD","ready_for_ldsc":False,
           "LDSC_h2_intercept_rg_MR_coloc_performed":False,
           "error_type":type(exc).__name__,"error_detail":str(exc)[:250]}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":r["technical_status"],
                      "outcome_counts":{x["phenotype"]:x["totals"]["allele_aligned_numeric_valid_pre_munge"] for x in r.get("outcomes",[])},
                      "scientific_status":r["scientific_status"],
                      "error":r.get("error_detail")},ensure_ascii=False))
    if r["technical_status"]!="PASS_ORIGINAL_MD5_PER_CHROMOSOME_PRE_MUNGE_COUNTS_ONLY":
        raise SystemExit(1)

if __name__=="__main__":
    main()
