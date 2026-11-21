#!/usr/bin/env python3
"""Independent GRCh37 rs198472 EUR frequency lookup; GIGASTROKE Neff HOLD.

Ensembl GRCh37 variant REST endpoint provides 1000 Genomes Phase 3 EUR
population allele frequencies. The original GIANT 2018 A1/tested allele
is independently recovered by a whole-publisher-SHA256-verified streaming
read before comparing EUR frequencies to the original GIGASTROKE source
diagnostic four phenotype variants. No inferred GC correction or N_eff use.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
from urllib.parse import urlsplit
from urllib.request import Request,urlopen

from scripts.ckm_stroke.audit_bmi_duplicate_zero_p import iter_original_rows

VARIANT="rs198472"
CHROM="11"
BP=61527266
ALLOWED_ALLELES={"A","G"}
POPULATION="1000GENOMES:phase_3:EUR"
API="https://grch37.rest.ensembl.org/variation/homo_sapiens/rs198472?pops=1"
MAX_API_BYTES=4_000_000
OUTCOMES=("AIS","LAA","CE","SVD")
CASES={"AIS":62100,"LAA":6399,"CE":10804,"SVD":6811}
CONTROLS=1234808
TOTALS={"AIS":1296908,"LAA":1241207,"CE":1245612,"SVD":1241619}
EXPECTED_SOURCE_REF=0.5341
EXPECTED_EAF={"AIS":0.7749,"LAA":0.7534,"CE":0.7575,"SVD":0.7588}

def verified_target_artifact(data):
    if (data.get("schema")!="ckm-GIGASTROKE-EUR-chr11-variant-identity-chr4-P-Wald-v1"
        or data.get("technical_status")!="PASS_PINNED_FULL_MD5_TARGET_CHR4_11_VARIANT_IDENTITY"
        or data.get("all_four_chr11_outlier_same_genomic_variant") is not True
        or data.get("full_source_reverified_four_original_MD5_and_SHA256") is not True
        or data.get("scientific_status")!="HOLD"
        or data.get("outcome_per_variant_N_eff_approved") is not False):
        raise ValueError("previous full original source MD5 rs198472 identity not proven")
    records=data.get("chr11_source_outlier_frequency_cross_trait_details")
    if (not isinstance(records,list) or len(records)!=4
        or {x.get("chromosome") for x in records}!={11}
        or {x.get("position_GRCh37") for x in records}!={BP}
        or {x.get("rsid_from_pinned_GIANT") for x in records}!={VARIANT}
        or {x.get("source_reference_alleles_sorted") for x in records}!={"AG"}
        or {round(float(x.get("GIANT_HRS_tested_allele_frequency")),4) for x in records}!={EXPECTED_SOURCE_REF}):
        raise ValueError("rs198472 phenotype source allele/frequency drift")
    vals=[float(x["GIGASTROKE_frequency_oriented_to_GIANT"]) for x in records]
    if sorted(round(x,4) for x in vals)!=sorted(EXPECTED_EAF.values()):
        raise ValueError("original four outcome EAF reference values drift")
    if data.get("chr11_outlier_variant_identity")!=records[0]:
        raise ValueError("source subtype original ID inconsistency")
    if len({(r["rsid_from_pinned_GIANT"],r["position_GRCh37"]) for r in records})!=1:
        raise ValueError("source rsID/coordinate inconsistency")
    return records

def verify_giant_tested_allele(original=iter_original_rows,expected_source_rows=2336269):
    hits=[]
    count=0
    verified=0
    for lineno,row in original():
        if lineno==0:
            if row[0]!="__SOURCE_COMPLETE_SHA256_PASS__":
                raise ValueError("GIANT 2018 original independent publisher SHA sentinel drift")
            verified+=1
            break
        count+=1
        if row[2]==VARIANT or (row[0]==CHROM and row[1]==str(BP)):
            if row[2]!=VARIANT or row[0]!=CHROM or row[1]!=str(BP):
                raise ValueError("GIANT rs198472 identity conflicts with GRCh37 coordinate")
            if {row[3],row[4]}!=ALLOWED_ALLELES:
                raise ValueError("GIANT original rs198472 expected A/G alleles changed")
            freq=float(row[5])
            if not math.isfinite(freq) or abs(freq-EXPECTED_SOURCE_REF)>1e-6:
                raise ValueError("GIANT original HRS reference frequency differs from source audit")
            hits.append({"source":"original_GIANT_Yengo2018_compressed_SHA256_verified",
                         "GRCh37_chromosome":int(CHROM),"GRCh37_position":BP,
                         "rsid":VARIANT,"tested_allele":row[3],
                         "other_allele":row[4],"HRS_tested_allele_frequency":freq})
    else:
        raise ValueError("GIANT original SHA sentinel absent: source stream incomplete")
    if verified!=1 or count!=expected_source_rows or len(hits)!=1:
        raise ValueError("GIANT publisher original source row count or rsID uniqueness drift")
    return hits[0]

def parse_ensembl_variant(data,giant):
    if data.get("name")!=VARIANT:
        raise ValueError("Ensembl exact rsID not returned")
    mappings=data.get("mappings")
    if not isinstance(mappings,list):
        raise ValueError("Ensembl GRCh37 mapping list absent")
    confirmed=[x for x in mappings if (x.get("assembly_name")=="GRCh37"
               and str(x.get("seq_region_name"))==CHROM
               and x.get("start")==BP and x.get("end")==BP
               and set(str(x.get("allele_string","")).split("/"))==ALLOWED_ALLELES)]
    if len(confirmed)!=1:
        raise ValueError("Ensembl rs198472 GRCh37 coordinate/allele validation failed")
    populations=data.get("populations")
    if not isinstance(populations,list):
        raise ValueError("Ensembl population data absent")
    observed={}
    for p in populations:
        if p.get("population")!=POPULATION:
            continue
        allele=p.get("allele")
        value=p.get("frequency")
        if allele not in ALLOWED_ALLELES or allele in observed or (
            type(value) not in (int,float) or not math.isfinite(value)
            or not 0<=value<=1):
            raise ValueError("invalid or duplicated European 1000 Genomes allele population frequency")
        observed[allele]=value
    if set(observed)!=ALLOWED_ALLELES or abs(sum(observed.values())-1)>0.02:
        raise ValueError("missing/inconsistent biallelic 1000 Genomes EUR A/G data")
    tested=giant["tested_allele"]
    if tested not in observed or giant["other_allele"] not in observed:
        raise ValueError("GIANT A1 not represented in independent EUR populations")
    return {
        "population":POPULATION,
        "original_GIANT_tested_allele":tested,
        "eur_AF_of_GIANT_tested_allele":round(observed[tested],7),
        "EUR_AF_allele_A":round(observed["A"],7),
        "EUR_AF_allele_G":round(observed["G"],7),
        "European_1000G_frequency_is_independent_population_reference":True,
        "GRCh37_alleles_and_coordinate_verified":True,
        "source":"Ensembl_GRCh37_REST_1000GENOMES_phase3_EUR",
    }

def retrieve_ensembl(opener=urlopen):
    req=Request(API,headers={
        "Accept":"application/json",
        "User-Agent":"CKM-Stroke-GIGASTROKE-independent-1000G-EUR/1.0"})
    with opener(req,timeout=45) as response:
        loc=urlsplit(response.geturl())
        if getattr(response,"status",None)!=200 or loc.scheme!="https" or loc.hostname!="grch37.rest.ensembl.org":
            raise ValueError("Ensembl GRCh37 reference source host or HTTP status invalid")
        if response.headers.get("Content-Type","").lower().find("json")<0:
            raise ValueError("Ensembl variation endpoint did not return JSON")
        payload=response.read(MAX_API_BYTES+1)
        if len(payload)>MAX_API_BYTES:
            raise ValueError("Ensembl reference API exceeds source size cap")
        return json.loads(payload.decode("utf-8"))

def effective_n_scenario():
    """Population sample-count illustration only; not GIGASTROKE SNP-wise N."""
    samples=[]
    for ph in OUTCOMES:
        cases,controls=CASES[ph],CONTROLS
        if cases+controls!=TOTALS[ph]:
            raise ValueError("reported case/control and total N do not reconcile")
        demo=4*cases*controls/(cases+controls)
        samples.append({
            "phenotype":ph,
            "reported_cases":cases,
            "reported_controls":controls,
            "publisher_reported_total_N":TOTALS[ph],
            "illustrative_4_case_control_harmonic_N_eff":round(demo,2),
            "illustrative_value_NOT_original_SNP_level_N_eff":True,
        })
    return {
        "source":"Mishra2022_GIGASTROKE_EUR_44_contributing_cohort_meta_analysis",
        "meta_analysis":"inverse_variance_weighted_multi_cohort",
        "original_public_outcome_source_column_count":11,
        "per_variant_case_N_available":False,
        "per_variant_control_N_available":False,
        "per_variant_N_eff_available":False,
        "input_cohort_allele_availability_heterogeneity_resolved":False,
        "constant_total_N_approvable_for_LDSC":False,
        "constant_illustrative_N_eff_approvable_for_LDSC":False,
        "original_cov_LDSC_reference_method_matches_our_pipeline_unverified":True,
        "illustrative_case_control_harmonic_N_eff_by_trait":samples,
        "N_eff_scientific_status":"HOLD",
    }

def audit(original_target,reference,giant_original):
    records=verified_target_artifact(original_target)
    giant=giant_original
    if giant.get("rsid")!=VARIANT or giant.get("tested_allele") not in ALLOWED_ALLELES or giant.get("other_allele") not in ALLOWED_ALLELES:
        raise ValueError("GIANT tested allele identity missing")
    independent=parse_ensembl_variant(reference,giant)
    freq=independent["eur_AF_of_GIANT_tested_allele"]
    comparisons=[]
    for trait,rec in zip(OUTCOMES,records):
        result=round(float(rec["GIGASTROKE_frequency_oriented_to_GIANT"])-freq,7)
        comparisons.append({
            "phenotype":trait,
            "GIGASTROKE_EAF_oriented_to_GIANT_tested":float(rec["GIGASTROKE_frequency_oriented_to_GIANT"]),
            "independent_1000G_EUR_AF_same_tested_allele":freq,
            "absolute_frequency_difference_GIGASTROKE_vs_1000G_EUR":round(abs(result),7),
            "signed_difference_GIGASTROKE_minus_1000G_EUR":result,
        })
    return {
        "schema":"ckm-GIGASTROKE-rs198472-independent-1000G-EUR-Neff-evidence-v1",
        "technical_status":"PASS_GRCh37_1000G_EUR_SOURCE_BIALLELIC_VARIANT_MATCH",
        "GRCh37_reference_rsID":VARIANT,
        "GIANT_original_source":giant,
        "independent_EUR_population_frequency":independent,
        "GIGASTROKE_allele_aligned_subtype_comparisons":comparisons,
        "GIGASTROKE_source_chr11_exact_identity_from_publisher_MD5":True,
        "frequency_difference_implies_genotyping_or_ancestry_error":False,
        "source_allele_effect_data_uploaded":False,
        "LDSC_effective_sample_size_evidence":effective_n_scenario(),
        "scientific_status":"HOLD","ready_for_ldsc":False,
        "LDSC_h2_rg_MR_coloc_executed":False,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--verified-original-target",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    try:
        src=json.loads(args.verified_original_target.read_text("utf-8"))
        giant=verify_giant_tested_allele()
        reference=retrieve_ensembl()
        result=audit(src,reference,giant)
    except Exception as exc:
        result={
            "technical_status":"HOLD_GRCh37_1000G_EUR_ORIGINAL_GIANT_SOURCE_UNVERIFIED",
            "error_type":type(exc).__name__,
            "error_detail":str(exc)[:300],
            "scientific_status":"HOLD","ready_for_ldsc":False,
            "LDSC_h2_rg_MR_coloc_executed":False,
            "LDSC_effective_sample_size_evidence":effective_n_scenario(),
        }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "technical_status":result["technical_status"],
        "GIANT_tested_allele":result.get("GIANT_original_source",{}).get("tested_allele"),
        "EUR_AF_same_allele":result.get("independent_EUR_population_frequency",{}).get("eur_AF_of_GIANT_tested_allele"),
        "N_eff_status":result["LDSC_effective_sample_size_evidence"]["N_eff_scientific_status"],
        "science":result["scientific_status"],
        "error":result.get("error_detail"),
    }))
    if result["technical_status"]!="PASS_GRCh37_1000G_EUR_SOURCE_BIALLELIC_VARIANT_MATCH":
        raise SystemExit(1)

if __name__=="__main__":
    main()
