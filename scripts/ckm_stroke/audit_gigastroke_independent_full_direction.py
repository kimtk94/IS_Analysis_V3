#!/usr/bin/env python3
"""Fail-closed independent Locke-2015 effect allele orientation review on full GWAS."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from scripts.ckm_stroke.review_gigastroke_full_evidence import PINNED,review
from scripts.ckm_stroke.audit_ldsc_chr_attrition import audit_attrition
from scripts.ckm_stroke.audit_ldsc_preregression_readiness import EXPECTED

VALID=("2015_same","2015_swapped","2015_complement_same","2015_complement_swapped",
       "2015_allele_incompatible","2015_non_ACGT_alleles","2015_same_allele_invalid",
       "2015_palindromic_ambiguous","not_in_locke_2015","yengo_numeric_or_allele_invalid")
MATCHED=VALID[:4]

def nn(x,label):
    if type(x) is not int or x<0:
        raise ValueError("invalid full-source direction count: "+label)
    return x

def direction_counts(source):
    review(source)
    audit_attrition(source)
    old=source.get("independent_locke_panel_evidence",{})
    if (old.get("candidate_rows")!=EXPECTED["locke_2015"]["rows"] or
        old.get("plain_sha256")!=EXPECTED["locke_2015"]["sha256"] or
        old.get("fractional_N_rows")!=11506 or
        old.get("selected_without_Yengo_2018") is not True):
        raise ValueError("independent Locke 2015 pinned source mismatch")
    cross=source.get("independent_2015_vs_2018_order_evidence",{})
    cat=cross.get("categories",{})
    if (cross.get("total_Yengo_candidate")!=1007899 or
        cross.get("both_original_BMI_order_comparable")!=1007858 or
        cross.get("2015_order_reversed_vs_2018")!=491428 or
        cross.get("individual_participant_independence_confirmed") is not False or
        cat.get("same")!=516430 or cat.get("swapped")!=491428 or
        cat.get("not_in_independent_2015")!=41 or
        sum(v for v in cat.values() if type(v) is int)!=1007899):
        raise ValueError("pinned independent 2015/2018 allele orders differ")
    result=[]
    for row in source["outcomes"]:
        ph=row["phenotype"]
        check=row.get("locke_2015_independent_full_orientation",{})
        if (row.get("raw_md5_verified") is not True or
            check.get("complete_original_source_scanned") is not True or
            check.get("source_MD5_validated_by_caller") is not True or
            check.get("approved_LDSC_analysis") is not False):
            raise ValueError("source full MD5 or 2015 direction verification missing")
        original_n=row["counts"]["unique_outcome_coordinate_matches"]
        if check.get("unique_original_outcome_candidate_positions")!=original_n:
            raise ValueError("original SNP count versus independent direction drift")
        raw=check.get("categories")
        if not isinstance(raw,dict) or set(raw)-set(VALID):
            raise ValueError("unrecognized original independent direction category")
        c={k:nn(raw.get(k,0),k) for k in VALID}
        if sum(c.values())!=original_n:
            raise ValueError("not all full original matched SNPs classified")
        aligned=sum(c[k] for k in MATCHED)
        reversed_n=c["2015_swapped"]+c["2015_complement_swapped"]
        if (aligned!=check.get("aligned_to_both_2015_and_2018") or
            reversed_n!=check.get("2015_effect_sign_reversals") or aligned<900000 or
            c["2015_same"]==0 or reversed_n==0):
            raise ValueError("independent full outcome allele reversal ledger mismatch")
        if any(c[k] for k in VALID[4:8]):
            raise ValueError("independent original allele incompatibility requires review")
        by_chr=check.get("by_chr")
        if (not isinstance(by_chr,list) or len(by_chr)!=22 or
            [item.get("chr") for item in by_chr]!=list(range(1,23))):
            raise ValueError("22-chromosome independent direction ledger invalid")
        c22={k:0 for k in VALID}
        output_chr=[]
        for item in by_chr:
            vals=item.get("categories")
            if not isinstance(vals,dict) or set(vals)-set(VALID):
                raise ValueError("unrecognized per-chromosome classification")
            cur={k:nn(vals.get(k,0),k) for k in VALID}
            if item.get("candidate_positions")!=sum(cur.values()):
                raise ValueError("chromosome independent direction count mismatch")
            for k,v in cur.items():c22[k]+=v
            output_chr.append({"chr":item["chr"],"aligned":sum(cur[k] for k in MATCHED),
                               "sign_reversed":cur["2015_swapped"]+cur["2015_complement_swapped"],
                               "2015_not_present":cur["not_in_locke_2015"]})
        if c22!=c:
            raise ValueError("22-chromosome direction totals differ from full original")
        result.append({"phenotype":ph,"accession":PINNED[ph][0],
                       "original_full_publisher_MD5_verified":True,
                       "yengo_reference_matched_source_coordinates":original_n,
                       "also_aligned_to_independent_locke_2015":aligned,
                       "2015_effect_allele_order_reversals":reversed_n,
                       "2015_allele_order_reversal_pct":round(100*reversed_n/aligned,4),
                       "same_count":c["2015_same"],"swapped_count":c["2015_swapped"],
                       "strand_complement_same":c["2015_complement_same"],
                       "strand_complement_swapped":c["2015_complement_swapped"],
                       "by_chr":output_chr})
    if {x["phenotype"] for x in result}!=set(PINNED):
        raise ValueError("four full original GIGASTROKE outcomes missing")
    return {"schema":"ckm-gigastroke-four-full-locke2015-original-effect-order-v1",
            "technical_status":"PASS_FULL_ORIGINAL_MD5_LOCKE_2015_ALLELE_ORDER",
            "Locke2015_Yengo2018_independently_shared":1007858,
            "Locke2015_Yengo2018_independent_reversed":491428,
            "outcomes":result,
            "source_effect_allele_order_only_not_causal_effect":True,
            "cross_study_frequency_equivalence_verified":False,
            "source_EUR_GWAS_original_effects_exported":False,
            "approved_outcome_N_eff":False,
            "LDSC_h2_rg_MR_coloc_executed":False,
            "scientific_status":"HOLD","ready_for_ldsc":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--evidence",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    try:r=direction_counts(json.loads(a.evidence.read_text("utf-8")))
    except Exception as exc:
        r={"technical_status":"HOLD_UNVERIFIED_INDEPENDENT_FULL_DIRECTION",
           "error_type":type(exc).__name__,"error_detail":str(exc)[:240],
           "scientific_status":"HOLD","ready_for_ldsc":False,
           "approved_outcome_N_eff":False,"LDSC_h2_rg_MR_coloc_executed":False}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":r["technical_status"],
                      "reversed":{x["phenotype"]:x["2015_effect_allele_order_reversals"] for x in r.get("outcomes",[])},
                      "science":r["scientific_status"],"error":r.get("error_detail")}))
    if r["technical_status"]!="PASS_FULL_ORIGINAL_MD5_LOCKE_2015_ALLELE_ORDER":
        raise SystemExit(1)
if __name__=="__main__":
    main()
