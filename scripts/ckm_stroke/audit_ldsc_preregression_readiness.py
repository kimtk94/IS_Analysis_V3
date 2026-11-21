#!/usr/bin/env python3
"""Checksum-gated descriptive Z^2/mean-chi2 BMI audit and GIGASTROKE h2 input HOLD.

This is not LD Score regression. Never writes per-SNP effects, estimates h2,
fits an intercept, chooses outcome N_eff or computes rg/MR.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
from statistics import median

EXPECTED={
 "locke_2015":{
  "columns":("SNP","N_REPORTED","Z_FROM_BETA_SE","A1","A2"),
  "rows":1070753,
  "sha256":"a1cc5ec06cc81d2bb0a583eb76c0f7e81784cd5dc3f622270de61023f95de4bf",
  "status":"PASS_2015_ORIGINAL_INDEPENDENT_SNP_PANEL_ONLY",
  "rows_key":"selected_rows","sha_key":"2015_candidate_plaintext_sha256",
 },
 "yengo_2018":{
  "columns":("SNP","N","Z","A1","A2"),
  "rows":1007899,
  "sha256":"ed8f713827cd3281bea5f88a9cb42a3eabfecb1ee67c912ca7e86f151e4399bf",
  "status":"PASS_CONSERVATIVE_CANDIDATE_CREATED_NOT_LDSC_APPROVED",
  "rows_key":"candidate_rows","sha_key":"candidate_plaintext_sha256",
 },
}
NUL_CHI2_MEDIAN=0.454936423119572
EXPECTED_OUTCOMES={
 "AIS":("GCST90104540",62100),
 "LAA":("GCST90104542",6399),
 "CE":("GCST90104541",10804),
 "SVD":("GCST90104543",6811),
}
EXPECTED_CONTROLS=1234808

# Both BMI pipelines must trace back to the same pinned, ancestry-matched
# reference bytes. This does NOT validate the final LDSC regression SNP count.
MRCIEU_HM3_SHA256 = "4ed5a031f88edff7f94663841160d0fff750c45d4d55a630515b4d6de082c49e"
ZENODO_EUR_LDSCORE_MD5 = "e2f16343c4cfaa76caa7d0c03d26b489"
ZENODO_EUR_LDSCORE_SHA256 = "9537f00eb0d163a935aaa2cf04b358b7cf21852279b9c7925802526f6060b069"
EXPECTED_THREE_WAY_REF_SNPS = 1172984


def reference_evidence_contract(name: str, provenance: dict):
    """Reject missing, swapped or mutated HM3/EUR LD reference provenance."""
    if name not in EXPECTED:
        raise ValueError("unknown BMI reference profile")
    sources = provenance.get("reference_sources" if name == "locke_2015" else "sources")
    if not isinstance(sources, dict):
        raise ValueError("missing pinned EUR reference/HM3 provenance")
    sha_key = "original_mrcieu_ref_sha256" if name == "locke_2015" else "mrcieu_sha256"
    pins = {
        sha_key: MRCIEU_HM3_SHA256,
        "zenodo_archive_md5": ZENODO_EUR_LDSCORE_MD5,
        "zenodo_archive_sha256": ZENODO_EUR_LDSCORE_SHA256,
    }
    if any(sources.get(k) != value for k, value in pins.items()):
        raise ValueError("pinned EUR reference/HM3 checksum mismatch")
    if name == "locke_2015":
        if (sources.get("strict_three_reference_rsid_count") != EXPECTED_THREE_WAY_REF_SNPS
            or sources.get("direct_unordered_allele_pair_agreement") != EXPECTED_THREE_WAY_REF_SNPS
            or sources.get("2018_BMI_candidate_used") is not False):
            raise ValueError("Locke HM3 reference SNP agreement or independence invalid")
    else:
        if (provenance.get("reference_triple_rsids") != EXPECTED_THREE_WAY_REF_SNPS
            or provenance.get("reference_equal_pair_rsids") != EXPECTED_THREE_WAY_REF_SNPS):
            raise ValueError("Yengo HM3 reference SNP agreement invalid")
    return {
        "ancestry": "EUR",
        "reference_build": "GRCh37",
        "reference_triple_direct_agreement_snps": EXPECTED_THREE_WAY_REF_SNPS,
        "M_5_50_source_covered": "22_chr_in_pinned_archive",
        "reference_sha256_match": True,
        "final_LDSC_regression_SNP_count_verified": False,
        "scientific_status": "HOLD",
    }


MAX_LINE=4096
MAX_BYTES=60000000
MAX_ROWS=1500000


def examine_candidate(name:str,candidate:Path,provenance:dict,*,strict=True):
    spec=EXPECTED[name]
    wanted_rows=spec["rows"] if strict else provenance["rows"]
    wanted_sha=spec["sha256"] if strict else provenance["plaintext_sha256"]
    if strict:
        source_count=provenance.get("counts",{}).get(spec["rows_key"])
        if (provenance.get("technical_status")!=spec["status"] or
            provenance.get("scientific_status")!="HOLD" or
            provenance.get("ready_for_ldsc") is not False or
            source_count!=spec["rows"] or
            provenance.get(spec["sha_key"])!=spec["sha256"]):
            raise ValueError("candidate provenance invalid for "+name)
        if name=="locke_2015" and provenance.get("2018_exposure_candidates_used") is not False:
            raise ValueError("Locke panel must be independently selected")
    reference_contract = reference_evidence_contract(name, provenance) if strict else None
    if candidate.stat().st_size>MAX_BYTES:raise ValueError("candidate compressed size exceeded cap")
    seen=set();nvals=[];chisq=[];sha=hashlib.sha256()
    count=0;fractional=0;sumchi=0.0;comp=0.0;large=0
    with gzip.open(candidate,"rb") as f:
        head=f.readline(MAX_LINE)
        sha.update(head)
        if tuple(head.rstrip(b"\r\n").decode("ascii").split("\t"))!=spec["columns"]:
            raise ValueError("candidate headers wrong / accidental LDSC-munging")
        for line in f:
            sha.update(line)
            if len(line)>MAX_LINE:raise ValueError("candidate line too large")
            row=line.rstrip(b"\r\n").decode("ascii").split("\t")
            if len(row)!=5:raise ValueError("candidate row width invalid")
            rsid,n,z,a1,a2=row
            if not rsid or rsid in seen:raise ValueError("duplicate candidate rsID")
            seen.add(rsid)
            if (len(a1)!=1 or len(a2)!=1 or a1 not in "ACGT" or a2 not in "ACGT" or
                a1==a2 or frozenset((a1,a2)) in (frozenset("AT"),frozenset("CG"))):
                raise ValueError("invalid candidate alleles")
            try:
                sample=float(n);z_value=float(z)
            except ValueError as exc:raise ValueError("invalid candidate N/Z") from exc
            if not math.isfinite(sample) or sample<=0 or not math.isfinite(z_value):
                raise ValueError("nonfinite or invalid N/Z")
            if name=="yengo_2018" and (not n.isdecimal() or not sample.is_integer()):
                raise ValueError("2018 candidate expected integral source N")
            chi=z_value*z_value
            if not math.isfinite(chi):raise ValueError("Z^2 overflow")
            count+=1
            if count>MAX_ROWS:raise ValueError("candidate row safety cap")
            if not sample.is_integer():fractional+=1
            if chi>80:large+=1
            nvals.append(sample);chisq.append(chi)
            y=chi-comp;t=sumchi+y;comp=(t-sumchi)-y;sumchi=t
    digest=sha.hexdigest()
    if count!=wanted_rows or digest!=wanted_sha:
        raise ValueError("candidate SHA/row pin mismatch for "+name)
    if not count:raise ValueError("empty candidate")
    chisq.sort();nvals.sort()
    def quantile(values,f):return values[int((len(values)-1)*f)]
    mean=sumchi/count
    chmedian=median(chisq)
    return {
        "name":name,"SNPs":count,"canonical_plaintext_sha256":digest,
        "reference_evidence": reference_contract or {"verification":"FIXTURE_ONLY_NOT_SOURCE_VERIFIED"},
        "reported_N":{"min":nvals[0],"p10":quantile(nvals,.1),"median":median(nvals),
                      "p90":quantile(nvals,.9),"max":nvals[-1],
                      "fractional_N_rows":fractional},
        "Z_definition":"2015_rounded_beta_over_se_GC_P_not_equivalent" if name=="locke_2015" else
                       "2018_GIANT_source_beta_over_se",
        "descriptive_only":{
            "mean_chi2":round(mean,8),"median_chi2":round(chmedian,8),
            "lambda_gc":round(chmedian/NUL_CHI2_MEDIAN,8),
            "chi2_gt80":large,
            "mean_chi2_gt_1_02":mean>1.02,
            "not_an_LDSC_h2_or_intercept":True
        },
        "regression_SNPs_gt_200000":count>200000,
        "regression_SNPs_gt_600000":count>600000,
        "scientific_status":"HOLD","ready_for_ldsc":False,
    }


def outcome_gates(cohort:dict,effects:dict,overlap:dict):
    if (cohort.get("technical_status")!="PASS_PRIMARY_CORRECTED_COHORT_AND_TABLE3_ONLY" or
        cohort.get("ready_for_ldsc") is not False or
        effects.get("technical_status")!="PASS_FOUR_BOUNDED_LOGOR_EFFECT_SCALE_ONLY" or
        effects.get("ready_for_ldsc") is not False or
        overlap.get("technical_status")!="PASS_FIRSTPARTY_TWO_SUPPLEMENT_COHORT_NAME_RECONCILIATION" or
        overlap.get("exact_individual_overlap_count_identified") is not False or
        overlap.get("direct_cohort_label_matches")!=10):
        raise ValueError("original phenotype, logOR or cohort-overlap evidence insufficient")
    outs=cohort.get("outcomes",[])
    ps=effects.get("prefixes",[])
    if len(outs)!=4 or len(ps)!=4:
        raise ValueError("4 corrected outcomes and effect probes required")
    mapped={o["phenotype"]:o for o in outs}
    probed={p["phenotype"]:p for p in ps}
    if set(mapped)!=set(EXPECTED_OUTCOMES) or set(probed)!=set(EXPECTED_OUTCOMES):
        raise ValueError("duplicate or missing subtype")
    output=[]
    for name in EXPECTED_OUTCOMES:
        o=mapped[name];p=probed[name];accession,cases=EXPECTED_OUTCOMES[name]
        ratio=o.get("original_CovLDSC_h2_SE_ratio")
        if (o["accession"]!=accession or p["accession"]!=accession or
            o["EUR_case_count_from_original_table1"]!=cases or
            o["EUR_control_count_from_original_table1"]!=EXPECTED_CONTROLS or
            p.get("source_per_variant_N_available") is not False or
            p.get("observed",{}).get("parsed_effect_rows",0)<500 or
            p.get("maximum_abs_beta_minus_lnOR",100)>1e-10 or
            not isinstance(ratio,(float,int)) or ratio<=0):
            raise ValueError("source phenotype/scale/N match invalid "+name)
        output.append({
            "phenotype":name,"accession":accession,
            "publisher_binary_beta_scale":"log_odds_ratio_prefix_checked",
            "original_prefix_variants_tested":p["observed"]["parsed_effect_rows"],
            "outcome_per_variant_N_in_source":False,
            "original_EUR_cases":cases,"original_EUR_controls":EXPECTED_CONTROLS,
            "published_CovLDSC_h2_over_SE":ratio,
            "published_is_not_new_LDSC":True,
            "published_h2_precision_below_precaution_ratio4":ratio<4.0,
            "total_first_Neff_not_approved":True,
            "cohortwise_Neff_not_approved":True,
            "individual_sample_overlap_unknown":True,
            "allow_constrain_cross_trait_intercept_to_zero":False,
            "input_gate":"HOLD_CASE_CONTROL_META_ANALYSIS_EFFECTIVE_N_NOT_VERIFIED",
            "LDSC_h2_intercept_rg_NOT_RUN":True
        })
    return output


def audit(locke:Path,locke_evidence:Path,yengo:Path,yengo_evidence:Path,
          cohort_file:Path,probe_file:Path,overlap_file:Path):
    report={
        "schema":"ckm-2015-2018-BMI-and-stroke-h2-preregression-preflight-v1",
        "checked_utc":datetime.now(timezone.utc).isoformat(),
        "technical_status":"HOLD_SOURCE_OR_SCIENTIFIC_INPUT_EVIDENCE",
        "scientific_status":"HOLD","ready_for_ldsc":False,
        "LDSC_h2_or_intercept_run":False,
        "LDSC_rg_MR_coloc_run":False,
        "cross_trait_intercept_zero_fixed":False,
    }
    try:
        p=json.loads(locke_evidence.read_text("utf-8"))
        q=json.loads(yengo_evidence.read_text("utf-8"))
        a=examine_candidate("locke_2015",locke,p)
        b=examine_candidate("yengo_2018",yengo,q)
        cohorts=json.loads(cohort_file.read_text("utf-8"))
        prefixes=json.loads(probe_file.read_text("utf-8"))
        overlap=json.loads(overlap_file.read_text("utf-8"))
        outcomes=outcome_gates(cohorts,prefixes,overlap)
        report.update({
            "technical_status":"PASS_INPUT_DESCRIPTIVE_BMI_CHI2_OUTCOME_N_GATE_ONLY",
            "BMI":[a,b],"stroke_outcomes":outcomes,
            "diagnostic_rules":{
                "LDSC_SNP_count_warning_below":200000,
                "preferred_regression_SNP_count_at_least":600000,
                "review_mean_chi2_floor":1.02,
                "published_h2_over_SE_precaution":4.0,
                "source":"LDSC FAQ; these are review thresholds, not scientific GO"
            },
            "remaining_scientific_gates":[
                "2015 BMI rounded beta/se, GC P discrepancy and fractional SNP N must be treated explicitly",
                "GIGASTROKE per-SNP N unavailable and meta-analysis effective N/weights unsettled",
                "Locke and GIGASTROKE share study-cohort names; exact participants unquantified",
                "No finalized LDSC regression-specific allele, chi2 and SNP attrition protocol",
                "Run single-trait h2/mean-chi2/intercept with approved inputs before any rg"
            ],
        })
    except Exception as e:
        report["error_type"]=type(e).__name__
        report["error_detail"]=str(e)[:240]
    return report


def main():
    ap=argparse.ArgumentParser()
    for key in ("locke","locke_evidence","yengo","yengo_evidence",
                "cohort_file","probe_file","overlap_file"):
        ap.add_argument("--"+key.replace("_","-"),type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    report=audit(*(getattr(args,k) for k in (
        "locke","locke_evidence","yengo","yengo_evidence",
        "cohort_file","probe_file","overlap_file")))
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":report["technical_status"],
        "BMI":[{"name":x["name"],"N":x["SNPs"],
                "mean_chi2":x["descriptive_only"]["mean_chi2"]}
               for x in report.get("BMI",[])],
        "subtype_N_hold":len(report.get("stroke_outcomes",[])),
        "scientific":report["scientific_status"],
        "error":report.get("error_detail"),
    },ensure_ascii=False))
    if report["technical_status"]!="PASS_INPUT_DESCRIPTIVE_BMI_CHI2_OUTCOME_N_GATE_ONLY":
        raise SystemExit(1)


if __name__=="__main__":main()
