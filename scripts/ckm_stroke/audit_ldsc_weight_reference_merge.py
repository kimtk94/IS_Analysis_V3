#!/usr/bin/env python3
"""Actual EUR LDSC reference-score / weight-file SNP intersection preflight.

Rebuild Yengo2018 from pinned GIANT/HM3/Zenodo archive in its own workflow step.
Independently stream the exact canonical 22 chromosome original EUR LD score
files, both M and M_5_50, and w_hm3.snplist from SHA+MD5 pinned Zenodo. Assert
candidate rsID overlap *per chromosome* agrees with previous MD5-verified
GIGASTROKE original 4-outcome chromosome pre-munge SNP reference counts.

These are available reference-weight SNPs, not actual post-munge LDSC outcome
regression SNPs. Full GIGASTROKE per-SNP N_eff remains unverified.
"""
from __future__ import annotations
import argparse
from collections import Counter
import gzip
import hashlib
import json
import math
from pathlib import Path
import tarfile
import urllib.request
from urllib.parse import urlsplit

from scripts.ckm_stroke.verify_bmi_ldsc_candidate import inspect
from scripts.ckm_stroke.audit_ldsc_chr_attrition import audit_attrition
from scripts.ckm_stroke.audit_ldsc_preregression_readiness import (
    EXPECTED, EXPECTED_OUTCOMES, EXPECTED_CONTROLS
)
from scripts.ckm_stroke.audit_eur_ldscore_zenodo import (
    FILE_URL,SIZE,EXPECTED_MD5,MAX_ARCHIVE_BYTES,DigestReader,safe_member,
    LD_SCORE_PATTERN,M_PATTERN,M_5_50_PATTERN,read_record
)
from scripts.ckm_stroke.audit_hm3_zenodo_catalog_diff import (
    BUNDLED_NAME,BUNDLED_SIZE,parse_ldscore_member,read_bundle_hm3
)

EXPECTED_SHA256="9537f00eb0d163a935aaa2cf04b358b7cf21852279b9c7925802526f6060b069"
EXPECTED_SUPPLEMENT_SHA256="7c688f864c73ec2e9624bde9f53ac48a9880a6a08e59f5e1266690e243618ba4"
EXPECTED_SNP_COUNT=1007899
TRAITS=("AIS","LAA","CE","SVD")
CHROMOSOMES=set(range(1,23))

def verified_candidate(candidate:Path,provenance:dict):
    r=inspect(candidate,provenance)
    if (r["candidate_rows_verified"]!=EXPECTED_SNP_COUNT or
        r["candidate_plaintext_sha256"]!=EXPECTED["yengo_2018"]["sha256"] or
        provenance.get("sources",{}).get("zenodo_archive_sha256")!=EXPECTED_SHA256 or
        provenance.get("ready_for_ldsc") is not False):
        raise ValueError("original independently pinned GIANT BMI candidate missing")
    seen=set()
    with gzip.open(candidate,"rt",encoding="ascii",newline="") as fd:
        if fd.readline().rstrip("\r\n")!="SNP\tN\tZ\tA1\tA2":
            raise ValueError("candidate GWAS column order mismatch")
        for line in fd:
            name=line.split("\t",1)[0]
            if name in seen:raise ValueError("candidate duplicate rsid")
            seen.add(name)
    if len(seen)!=EXPECTED_SNP_COUNT:
        raise ValueError("candidate original SNP count wrong")
    return seen

def archive_reference(opener=urllib.request.urlopen,*,required_meta=True):
    """Strict full-archive scan, no disk extraction and no reused SNP contents."""
    if required_meta:read_record(opener)
    req=urllib.request.Request(FILE_URL,headers={"User-Agent":"CKM-Stroke-LDSC-exact-weight-reference/1.0"})
    with opener(req,timeout=60) as http:
        url=urlsplit(http.geturl())
        if getattr(http,"status",None)!=200 or url.scheme!="https" or url.hostname!="zenodo.org":
            raise ValueError("pinned original LD-score source host/HTTP status drift")
        if http.headers.get("Content-Length")!=str(SIZE):
            raise ValueError("pinned archive Content-Length changed")
        reader=DigestReader(http,MAX_ARCHIVE_BYTES)
        by_chr={}
        hm3=None
        counts={"M":{},"M_5_50":{}}
        members=0
        with tarfile.open(fileobj=reader,mode="r|gz") as tar:
            for item in tar:
                members+=1
                if members>130:raise ValueError("pinned LDSC archive member count cap exceeded")
                safe_member(item.name)
                if not(item.isfile() or item.isdir()):
                    raise ValueError("unexpected special archive member")
                if item.isdir():continue
                source=tar.extractfile(item)
                if source is None:raise ValueError("missing source tar member")
                if item.name==BUNDLED_NAME:
                    if hm3 is not None or item.size!=BUNDLED_SIZE:
                        raise ValueError("mutated HM3 SNP-list reference")
                    hm3,meta=read_bundle_hm3(source)
                elif LD_SCORE_PATTERN.fullmatch(item.name):
                    ch=int(LD_SCORE_PATTERN.fullmatch(item.name).group(1))
                    if ch in by_chr:raise ValueError("duplicated chromosome LD scores")
                    ids,qc=parse_ldscore_member(source,ch)
                    if qc["duplicates"]!=0:raise ValueError("duplicated LDscore rsID")
                    by_chr[ch]=ids
                else:
                    kind=None;ch=None
                    for label,pat in (("M",M_PATTERN),("M_5_50",M_5_50_PATTERN)):
                        match=pat.fullmatch(item.name)
                        if match:
                            kind=label
                            ch=int(match.group(1))
                            break
                    if kind:
                        if ch in counts[kind]:raise ValueError("duplicate LD-score M normalizer")
                        values=source.read(256).strip().split()
                        if len(values)!=1:raise ValueError("invalid M reference scalar count")
                        scalar=float(values[0])
                        if not math.isfinite(scalar) or scalar<=0:
                            raise ValueError("invalid positive LDSC M normalizer")
                        counts[kind][ch]=scalar
        while reader.read(1048576):pass
    diagnosis={
        "actual_archive_bytes":reader.nbytes,
        "expected_archive_bytes":SIZE,
        "publisher_MD5_match":reader.md5.hexdigest()==EXPECTED_MD5,
        "publisher_SHA256_match":reader.sha256.hexdigest()==EXPECTED_SHA256,
        "LDscore_chromosomes":sorted(by_chr),
        "M_chromosomes":sorted(counts["M"]),
        "M_5_50_chromosomes":sorted(counts["M_5_50"]),
        "HM3_present":hm3 is not None,
        "archive_tar_member_count":members,
    }
    if (reader.nbytes!=SIZE or reader.md5.hexdigest()!=EXPECTED_MD5 or
        reader.sha256.hexdigest()!=EXPECTED_SHA256 or set(by_chr)!=CHROMOSOMES or
        set(counts["M"])!=CHROMOSOMES or set(counts["M_5_50"])!=CHROMOSOMES or
        hm3 is None):
        raise ValueError("publisher 22-chromosome LD reference checksum/member completeness drift: "+
                         json.dumps(diagnosis,sort_keys=True))
    all_ids=set()
    for ch,ids in by_chr.items():
        if all_ids & ids:
            raise ValueError("reference rsID assigned to multiple chromosomes")
        all_ids.update(ids)
    return by_chr,set(hm3),{"record":"Zenodo 8182036","publisher_archive_bytes":reader.nbytes,
        "publisher_archive_md5_verified":True,"publisher_archive_sha256_verified":True,
        "reference_LD_score_chromosomes":22,"M_and_M_5_50_chromosomes":22,
        "total_reference_LD_score_unique_SNPs":len(all_ids),
        "HM3_reference_unique_SNPs":len(hm3),
        "M_5_50_per_chromosome":{str(k):counts["M_5_50"][k] for k in sorted(counts["M_5_50"])},
        "non_partitioned_weights_share_LD_score_file":True}

def verify_neff_source(cohort):
    if (cohort.get("technical_status")!="PASS_PRIMARY_CORRECTED_COHORT_AND_TABLE3_ONLY"
        or cohort.get("publisher_source_sha256")!=EXPECTED_SUPPLEMENT_SHA256 or
        cohort.get("publisher_original_table1_eur_cohort_count")!=44 or
        cohort.get("ready_for_ldsc") is not False or
        cohort.get("LDSC_h2_rg_executed") is not False):
        raise ValueError("corrected original stroke 44-cohort supplement provenance invalid")
    rows=cohort.get("outcomes")
    if not isinstance(rows,list) or len(rows)!=4 or [x.get("phenotype") for x in rows]!=list(TRAITS):
        raise ValueError("missing outcome in corrected original per-cohort source")
    out=[]
    for x in rows:
        ph=x["phenotype"]; accession,cases=EXPECTED_OUTCOMES[ph]
        c=x.get("illustrative_cohortwise_neff",{})
        pooled=x.get("illustrative_harmonic_neff")
        cohortwise=c.get("cohortwise_sum_harmonic_neff_illustrative")
        missing=c.get("cases_missing_cohort_control_count")
        if (x.get("accession")!=accession or x.get("EUR_case_count_from_original_table1")!=cases or
            x.get("EUR_control_count_from_original_table1")!=EXPECTED_CONTROLS or
            x.get("EUR_total_N_from_original_table1")!=cases+EXPECTED_CONTROLS or
            x.get("illustrative_neff_not_per_snp_or_approved") is not True or
            c.get("not_per_variant_N_and_not_approved_for_LDSC") is not True or
            type(missing) is not int or missing<0 or
            type(pooled) not in (int,float) or type(cohortwise) not in (int,float) or
            not (0<cohortwise<=pooled)):
            raise ValueError("original cohortwise meta-analysis N evidence mismatch: "+ph)
        ratio=(pooled-cohortwise)/pooled*100
        out.append({
            "phenotype":ph,"accession":accession,"N_total":cases+EXPECTED_CONTROLS,
            "pooled_case_control_harmonic_Neff_NOT_APPROVED":pooled,
            "cohortwise_known_controls_harmonic_Neff_NOT_APPROVED":cohortwise,
            "difference_percent":round(ratio,4),
            "number_original_cohorts_with_case_and_control_counts":c["cohort_count_with_cases_and_reported_controls"],
            "cases_in_cohorts_without_published_control_count":missing,
            "N_SNP_specific_observed":False,
            "scientific_status":"HOLD",
        })
    return out

def reconcile_reference(ids,by_chr,hm3,source,cohort):
    checked=audit_attrition(source)
    missing_hm3=ids-hm3
    score_union=set().union(*by_chr.values())
    missing_score=ids-score_union
    if missing_hm3 or missing_score:
        raise ValueError("pinned BMI candidate absent from independently verified HM3/22 chr LDscore")
    matches={i:ids & by_chr[i] for i in range(1,23)}
    if sum(len(v) for v in matches.values())!=len(ids):
        raise ValueError("reference weight-LD chromosome assignment does not reconcile")
    reference_chr=[r["candidate_reference_bmi_snps"] for r in source["outcomes"][0]["chr_attrition_pre_munge"]]
    bychr=[]
    for i in range(1,23):
        item={"chr":i,"pinned_reference_BMI_candidate_SNPs":len(matches[i]),
              "reference_LD_score_total_SNPs":len(by_chr[i]),
              "BMI_candidate_also_in_HM3_and_LD_score":len(matches[i])}
        if len(matches[i])!=reference_chr[i-1]:
            raise ValueError("original GIGASTROKE outcome reference BMI chr attrition disagrees with canonical LDscore archive")
        bychr.append(item)
    if (checked.get("source_BMI_candidate_count")!=EXPECTED_SNP_COUNT
        or len(checked.get("outcomes",[]))!=4):
        raise ValueError("prior source GIGASTROKE chromosome attrition not fully proven")
    return {"schema":"ckm-EUR-LDSC-reference-weight-and-outcome-N-source-v1",
        "technical_status":"PASS_22CHR_LD_REFERENCE_AND_SOURCE_NEFF_HOLD_AUDITED",
        "reference_BMI_candidate_count":len(ids),
        "reference_BMI_candidate_22chr_merge_count":sum(len(v) for v in matches.values()),
        "reference_BMI_candidate_also_in_HM3":len(ids),
        "chr_reference":bychr,
        "outcomes_pre_munge_counts":[
            {"phenotype":x["phenotype"],"valid_pre_munge_SNP_count":x["totals"]["allele_aligned_numeric_valid_pre_munge"]}
            for x in checked["outcomes"]],
        "cohort_meta_analysis_neff_illustrations":verify_neff_source(cohort),
        "case_control_N_eff_per_variant_known":False,
        "outcome_post_munge_SNP_selection_executed":False,
        "weight_file_also_available_as_ref_ld_file":True,
        "same_archive_weight_reference_not_independent_weight_validation":True,
        "reference_must_NOT_be_called_final_regression_SNP_count":True,
        "GC_intercept_and_h2_NOT_RUN":True,
        "rg_MR_coloc_NOT_RUN":True,
        "scientific_status":"HOLD","ready_for_ldsc":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--candidate",required=True,type=Path)
    p.add_argument("--provenance",required=True,type=Path)
    p.add_argument("--original-four-evidence",required=True,type=Path)
    p.add_argument("--corrected-cohort-evidence",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    try:
        orig=json.loads(a.original_four_evidence.read_text("utf-8"))
        cohort=json.loads(a.corrected_cohort_evidence.read_text("utf-8"))
        ids=verified_candidate(a.candidate,json.loads(a.provenance.read_text("utf-8")))
        by_chr,hm3,meta=archive_reference()
        result=reconcile_reference(ids,by_chr,hm3,orig,cohort)
        result["LD_score_source"]=meta
    except Exception as exc:
        result={"technical_status":"HOLD_LDSC_REFERENCE_OR_SAMPLE_COUNT_EVIDENCE",
          "error_type":type(exc).__name__,"error_detail":str(exc)[:240],
          "scientific_status":"HOLD","ready_for_ldsc":False,
          "case_control_N_eff_per_variant_known":False,"GC_intercept_and_h2_NOT_RUN":True}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["technical_status"],
          "LD_reference_intersection":result.get("reference_BMI_candidate_22chr_merge_count"),
          "outcomes":result.get("outcomes_pre_munge_counts"),
          "cohort_Neff_differences":{x["phenotype"]:x["difference_percent"]
                                   for x in result.get("cohort_meta_analysis_neff_illustrations",[])},
          "science":result["scientific_status"],"error":result.get("error_detail")}))
    if result["technical_status"]!="PASS_22CHR_LD_REFERENCE_AND_SOURCE_NEFF_HOLD_AUDITED":
        raise SystemExit(1)

if __name__=="__main__":
    main()
