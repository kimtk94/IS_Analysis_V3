#!/usr/bin/env python3
"""Consolidate GIGASTROKE original MD5-verified TPMI lead probe shards."""
import argparse
import json
from pathlib import Path
from scripts.ckm_stroke.review_gigastroke_full_evidence import PINNED,OBSERVED_COUNTS
TRAITS=("AIS","LAA","CE","SVD")
def consolidate(records):
    if set(records)!=set(TRAITS):
        raise ValueError("all four original subtype source results required")
    comparisons=[]
    for ph in TRAITS:
        r=records[ph]
        if (r.get("schema")!="ckm-TPMI43321-GIGASTROKE-four-subtype-original-variant-probe-v1"
            or r.get("technical_status")!="PASS_ORIGINAL_PUBLISHER_MD5_AND_EXACT_TARGET_LOCUS_PROBE"
            or r.get("trait")!=ph
            or r.get("original_source_accession")!=PINNED[ph][0]
            or r.get("original_publisher_MD5")!=PINNED[ph][1]
            or r.get("original_compressed_bytes")!=PINNED[ph][2]
            or r.get("original_all_rows_streamed")!=OBSERVED_COUNTS[ph][0]
            or r.get("original_MD5_verified") is not True
            or r.get("scientific_status")!="HOLD"
            or r.get("ready_for_ldsc") is not False
            or r.get("LDSC_MR_coloc_not_executed") is not True):
            raise ValueError("unverified original whole GWAS publisher MD5 or scientific gate for "+ph)
        count=r.get("exact_locus_original_row_count")
        assoc=r.get("one_locus_original_source_rows")
        if type(count) is not int or count<0 or count!=len(assoc):
            raise ValueError("invalid targeted locus row count")
        if count==0:
            if r.get("target_source_status")!="NOT_IN_ORIGINAL_SOURCE_GWAS":
                raise ValueError("missing variant absence was misclassified")
            category="NOT_OBSERVED_IN_EUR_SOURCE"
        elif count==1 and assoc[0].get("allele_aligned_effects_verified") is True:
            if r.get("target_source_status")!="PASS_ORIGINAL_ALLELE_DIRECTION_AND_NUMERIC":
                raise ValueError("source allele effect direction drift")
            category="OBSERVED_HARMONIZABLE_SOURCE_ONLY"
        else:
            category="HOLD_DUPLICATE_OR_UNHARMONIZABLE"
        comparisons.append({
          "phenotype":ph,"original_accession":r["original_source_accession"],
          "original_MD5_verified":True,
          "target_locus":r["searched_GRCh37_chrom"]+":"+str(r["searched_GRCh37_position"]),
          "original_row_count":count,
          "source_detected_category":category,
          "original_source_target_associations":assoc,
          "independent_replication_approved":False
        })
    return {
       "schema":"ckm-TPMI43321-GIGASTROKE-EUR-original-four-locus-compare-v1",
       "technical_status":"PASS_FOUR_GIGASTROKE_ORIGINAL_MD5_TARGET_LOCUS_DESCRIPTIVE_ONLY",
       "TPMI_source_rsID":"rs187759247",
       "TPMI_allele_C_GRCh37_chr19_pos":15573379,
       "outcomes":comparisons,
       "number_original_GIGASTROKE_subtypes_with_source_variant":sum(x["original_row_count"]>0 for x in comparisons),
       "any_actual_EUR_replication_test_executed":False,
       "transancestry_effect_beta_homogeneity_not_tested":True,
       "GIGASTROKE_EUR_and_TPMI_EAS_ancestry_LD_contexts_differ":True,
       "case_control_N_eff_verified":False,
       "scientific_replication_conclusion":"HOLD",
       "scientific_status":"HOLD","ready_for_ldsc":False,
       "LDSC_h2_rg_MR_coloc_executed":False
    }
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--folder",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    try:
        items={}
        for ph in TRAITS:
            matches=list(a.folder.glob(f"**/probe_{ph}.json"))
            if len(matches)!=1:raise ValueError("required original trait source evidence missing "+ph)
            items[ph]=json.loads(matches[0].read_text("utf-8"))
        result=consolidate(items)
    except Exception as exc:
        result={"technical_status":"HOLD_TARGETED_CROSSANCESTRY_FOUR_SOURCE_REVIEW",
                "error_type":type(exc).__name__,"error_detail":str(exc)[:250],
                "scientific_status":"HOLD","ready_for_ldsc":False}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n","utf-8")
    print(json.dumps({
        "technical":result["technical_status"],
        "traits":{x["phenotype"]:x["source_detected_category"] for x in result.get("outcomes",[])},
        "number_present":result.get("number_original_GIGASTROKE_subtypes_with_source_variant"),
        "error":result.get("error_detail"),
        "science":result["scientific_status"]}))
    if result["technical_status"]!="PASS_FOUR_GIGASTROKE_ORIGINAL_MD5_TARGET_LOCUS_DESCRIPTIVE_ONLY":
        raise SystemExit(1)
if __name__=="__main__":
    main()
