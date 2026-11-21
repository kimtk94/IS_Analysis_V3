#!/usr/bin/env python3
"""Create a strictly non-inferential GIGASTROKE lookup request from validated TPMI liftover.
Does not claim a GIGASTROKE effect or independence, and never fabricates N_eff.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

EXPECTED_SHA="ce14daef23f9fadf921030d9ed2429383d927f6da391c620d9abf870045b2e62"

def prepare(data):
    if (data.get("schema")!="ckm-TPMI43321-lead-variant-Ensembl-GRCh38-GRCh37-v1"
        or data.get("technical_status")!="PASS_BIDIRECTIONAL_ENSEMBL_BUILD_MAP_AND_REF_ALLELE"
        or data.get("original_source_SHA256_verified")!=EXPECTED_SHA
        or data.get("liftover_reverse_verified") is not True
        or data.get("scientific_status")!="HOLD"
        or data.get("ready_for_ldsc") is not False
        or data.get("GIGASTROKE_EUR_source_effect_at_target_observed") is not False):
        raise ValueError("mapped original TPMI lead has no trustworthy verified provenance")
    source=data.get("lead_GRCh38",{})
    mapped=data.get("mapped_GRCh37",{})
    chrom=str(mapped.get("chrom_GRCh37",""))
    pos=mapped.get("pos_GRCh37")
    ref=mapped.get("REF_GRCh37")
    alt=mapped.get("ALT_GRCh37")
    a1=mapped.get("A1_effect_allele_GRCh37")
    if (not chrom.isdecimal() or not 1<=int(chrom)<=22
        or type(pos) is not int or pos<1
        or not all(isinstance(a,str) and len(a)==1 and a in "ACGT" for a in (ref,alt,a1))
        or ref==alt or a1 not in (ref,alt)
        or mapped.get("mapped_GRCh37_reference_verified") is not True
        or mapped.get("allele_source_GRCh38_reference_verified") is not True):
        raise ValueError("mapped target locus, reference alleles or A1 not validated")
    if source.get("original_OBS_CT")!=204099:
        raise ValueError("original PLINK2 observation count changed")
    beta=source.get("logOR_original_A1")
    se=source.get("LOGOR_SE_original")
    p=source.get("P_original")
    if (not all(type(x) in (int,float) and math.isfinite(x) for x in (beta,se,p))
        or se<=0 or not 0<=p<=5e-8):
        raise ValueError("original TPMI A1 logOR, SE or significant P unavailable")
    return {
       "schema":"ckm-TPMI43321-GIGASTROKE-GRCh37-target-request-v1",
       "technical_status":"PASS_VERIFIED_TARGET_ONLY",
       "target_build":"GRCh37",
       "target_chr":chrom,"target_pos":pos,
       "target_REF":ref,"target_ALT":alt,"target_A1":a1,
       "target_variant_key":f"{chrom}:{pos}:{ref}:{alt}",
       "TPMI_reported_rsid":source.get("reported_rsid"),
       "TPMI_original_P":p,
       "TPMI_A1_logOR":beta,
       "TPMI_A1_logOR_SE":se,
       "TPMI_original_OBS_CT_not_Neff":204099,
       "GIGASTROKE_original_row_matched":False,
       "GIGASTROKE_EUR_effect_allele_harmonized":False,
       "population_independence_confirmed":False,
       "binary_case_control_Neff_approved":False,
       "scientific_status":"HOLD",
       "ready_for_ldsc":False,
       "mr_coloc_ldsc_executed":False
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--mapped",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    try:
        result=prepare(json.loads(a.mapped.read_text("utf-8")))
    except Exception as e:
        result={"technical_status":"HOLD_TPMI_GIGASTROKE_TARGET_UNVERIFIED",
                "error_type":type(e).__name__,"error_detail":str(e)[:240],
                "scientific_status":"HOLD","ready_for_ldsc":False}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n","utf-8")
    print(json.dumps({"status":result["technical_status"],"target":result.get("target_variant_key"),
                      "scientific_status":result["scientific_status"]}))
    if result["technical_status"]!="PASS_VERIFIED_TARGET_ONLY":
        raise SystemExit(1)

if __name__=="__main__":
    main()
