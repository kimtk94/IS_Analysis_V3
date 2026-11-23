#!/usr/bin/env python3
"""Audit pinned original coloc R source contract for binary GWAS beta+SE inputs.

The exact upstream source code is executed using base R in GitHub Actions,
without requiring any coloc package installation or simulated implementation.
- Case-control beta+varbeta path: 's' NOT REQUIRED by current check_dataset.
- Case-control pvalues+MAF+N approximation path: 's' REQUIRED.
Scientific status remains HOLD until QTL sdY/N, harmonized complete cis SNPs,
prior sensitivity and study comparability are reviewed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import urllib.parse
import urllib.request
from pathlib import Path

UPSTREAM_REPO = "chr1swallace/coloc"
UPSTREAM_COMMIT = "27776b1bb14e19d1542798ac7477fa378fc2fc0d"
SOURCE_BLOB_SHA1 = {
    "R/check.R": "7d244c2826ed36f363b986a26e4e0c612de9af1e",
    "R/claudia.R": "043783ac5ddec3477ff3b91e3e1faf271c51ce09",
}
UPSTREAM_RAW_BASE = f"https://raw.githubusercontent.com/{UPSTREAM_REPO}/{UPSTREAM_COMMIT}"
TEST_R = Path("tests/is/test_coloc_cc_beta_no_s_original_source.R")
EXPECTED_KEYS = (
    "cc_beta_varbeta_without_s",
    "cc_beta_varbeta_with_s_same_abf",
    "cc_pvalue_only_missing_s_rejected",
    "cc_pvalue_only_supplied_s_allowed",
    "quant_beta_varbeta_sdY_explicit_allowed",
    "quant_beta_varbeta_missing_sdY_N_MAF_rejected",
)
MAX_SOURCE_BYTES = 150_000


def require(ok, explanation):
    if not ok:
        raise ValueError(explanation)


def git_blob_sha1(content):
    blob = b"blob " + str(len(content)).encode("ascii") + b"\x00" + content
    return hashlib.sha1(blob, usedforsecurity=False).hexdigest()


def verify_source_code(source, expected_path):
    require(expected_path in SOURCE_BLOB_SHA1, "Unapproved coloc source path")
    raw=source.encode("utf-8")
    require(git_blob_sha1(raw)==SOURCE_BLOB_SHA1[expected_path],
            "Pinned original coloc R file SHA1 mismatch")
    stripped="\n".join(line for line in source.splitlines()
                       if not line.lstrip().startswith("#"))
    normalized=re.sub(r"\s+", "", stripped)
    if expected_path=="R/check.R":
        # 's' is only demanded INSIDE the no-(beta,varbeta) branch.
        if_condition='if(!("beta"%in%nd)||!("varbeta"%in%nd)){'
        s_condition='if(d$type=="cc"&&!("s"%in%nd))'
        branch=normalized.find(if_condition)
        s=normalized.find(s_condition)
        require(branch>=0 and s>branch and
                'require,s,proportionofsamples' in normalized[s:s+230] and
                normalized.find("p=pnorm(",s)>s,
                "Original coloc case fraction rule not conditional on missing effect estimates")
        require('if(d$type=="quant"&&!("sdY"%in%nd))' in normalized,
                "Original quantitative phenotype sdY condition missing")
    else:
        require('if("beta"%in%nd&&"varbeta"%in%nd)' in normalized,
                "Original coloc beta-varbeta pathway missing")
        require('approx.bf.estimates(z=d$beta/sqrt(d$varbeta)' in normalized,
                "Original effect-beta ABF route missing")
        require('if("pvalues"%in%nd&"MAF"%in%nd&"N"%in%nd)' in normalized,
                "Original p-value approximation branch missing")
        require('sd.prior<-if(type=="quant"){' in normalized and
                'else{effect_priors["cc"]}' in normalized,
                "Case-control prior differs or depends on phenotype sdY")
    return {"upstream_file":expected_path,
            "git_blob_sha1":SOURCE_BLOB_SHA1[expected_path],
            "source_utf8_bytes":len(raw),
            "original_code_semantics_verified":True}


def fetch_pinned_r_source(path, opener=urllib.request.urlopen):
    require(path in SOURCE_BLOB_SHA1, "Unknown upstream script")
    url=f"{UPSTREAM_RAW_BASE}/{path}"
    req=urllib.request.Request(url,headers={
        "User-Agent":"IS-Analysis-V3-coloc-original-cc-beta-variance-contract/1.0"
    })
    with opener(req,timeout=55) as response:
        require(response.status==200, "Upstream coloc original source HTTP status != 200")
        actual=urllib.parse.urlsplit(response.geturl())
        require(actual.scheme=="https" and actual.hostname=="raw.githubusercontent.com"
                and actual.path==urllib.parse.urlsplit(url).path,
                "Pinned coloc original source redirected outside exact immutable URL")
        body=response.read(MAX_SOURCE_BYTES+1)
    require(0<len(body)<=MAX_SOURCE_BYTES,"Pinned original coloc source exceeded bound")
    source=body.decode("utf-8")
    meta=verify_source_code(source,path)
    meta["original_immutable_url"]=url
    meta["sha256"]=hashlib.sha256(body).hexdigest()
    return body,meta


def parse_runtime_check(lines):
    require(isinstance(lines,str), "No genuine base-R coloc contract output")
    results={}
    for line in lines.splitlines():
        if not line.startswith("ORIGINAL_COLOC_TEST_"):
            continue
        label,sep,val=line.partition("=")
        require(sep=="=" and val=="PASS","Original coloc R source test failed")
        key=label.removeprefix("ORIGINAL_COLOC_TEST_")
        require(key in EXPECTED_KEYS and key not in results,
                "Unknown/duplicate original coloc runtime test")
        results[key]=True
    require(set(results)==set(EXPECTED_KEYS),
            "Original coloc source R tests absent/incomplete")
    return results


def audit_original_contract(manifest, source_directory, rscript="Rscript"):
    groups=manifest.get("groups")
    require(isinstance(groups,list) and len(groups)==36
            and len({(x["QTL_dataset"],x["gene_id"]) for x in groups})==36
            and len({x["gene_id"] for x in groups})==20
            and len([x for x in groups if x["QTL_dataset"]=="QTD000434"])==15,
            "P3A 36 original groups not source authenticated")
    require(TEST_R.is_file(),"Original-source R runtime tests not found")
    source_directory.mkdir(parents=True,exist_ok=True)
    provenance={}
    for original_path in SOURCE_BLOB_SHA1:
        payload, metadata=fetch_pinned_r_source(original_path)
        local=source_directory/original_path.replace("/","_")
        local.write_bytes(payload)
        provenance[original_path]=metadata
    cmd=[rscript,"--vanilla",str(TEST_R),
         str(source_directory/"R_check.R"),
         str(source_directory/"R_claudia.R")]
    run=subprocess.run(cmd,capture_output=True,text=True,
                       timeout=120,check=False)
    require(run.returncode==0,"Original coloc R code runtime test nonzero")
    tests=parse_runtime_check(run.stdout)
    return {
        "schema":"is-p3a-pinned-coloc-original-cc-beta-SE-casefraction-contract-v1",
        "technical_status":"PASS_ORIGINAL_COLOC_R_CHECK_AND_PROCESS_DATASET_EXECUTED",
        "scientific_status":"HOLD_QTL_MODEL_N_SDY_HARMONIZED_CIS_PRIORS",
        "upstream_original_coloc_repo":UPSTREAM_REPO,
        "upstream_original_coloc_commit":UPSTREAM_COMMIT,
        "pinned_original_R_source":provenance,
        "original_coloc_R_runtime_tests":tests,
        "n_original_coloc_R_checks_passed":len(tests),
        "exact_GIGASTROKE_EAS_source":"GCST90104545",
        "exact_GIGASTROKE_EUR_source":"GCST90104540",
        "GIGASTROKE_chr10_beta_log_original_OR_verified_in_prior_source":True,
        "GIGASTROKE_EAS_EUR_complete_chr10_effect_rows_verified_prior":725707,
        "original_case_control_outcome_type":"cc",
        "direct_beta_plus_varbeta_case_fraction_s_required_by_original_coloc":False,
        "pvalue_MAF_N_only_case_fraction_s_required_by_original_coloc":True,
        "source_case_control_fraction_s_observed_in_exact_original_GWAS_YAML":False,
        "source_case_fraction_cannot_be_invented_or_p_value_approximated":True,
        "cc_type_beta_SE_path_no_approved_case_N_eff_needed_for_ABF":True,
        "real_BBJ_independent_replication_permitted":False,
        "original_source_QTL_sdY_estimate_only_15_genes":True,
        "original_source_QTL_gene_regression_N_validated":False,
        "QTL_estimated_sdY_approved_for_formal_coloc":False,
        "cross_build_GRCh37_GRCh38_complete_QTL_GWAS_cis_harmonisation_approved":False,
        "complete_cis_intersected_SNP_coverage_approved":False,
        "source_prior_sensitivity_and_single_causal_variant_assumption_reviewed":False,
        "n_original_gene_tissue_assay_groups":36,
        "n_current_formal_coloc_approved":0,
        "n_current_formal_coloc_executed":0,
        "no_QTL_GWAS_SNP_betas_or_individual_sample_IDs_exported":True,
        "decision":"REMOVE_INCORRECT_CC_CASEFRACTION_S_AS_STANDALONE_HARD_BLOCKER_ONLY",
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--manifest",type=Path,required=True)
    parser.add_argument("--source-dir",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    result={
        "schema":"is-p3a-pinned-coloc-original-cc-beta-SE-casefraction-contract-v1",
        "technical_status":"HOLD_ORIGINAL_R_SOURCE_OR_EXECUTION",
        "scientific_status":"HOLD",
        "n_current_formal_coloc_approved":0,
        "n_current_formal_coloc_executed":0,
    }
    try:
        manifest=json.loads(args.manifest.read_text(encoding="utf-8"))
        result=audit_original_contract(manifest,args.source_dir)
    except (OSError,ValueError,TypeError,KeyError,UnicodeError,TimeoutError,
            subprocess.TimeoutExpired) as exc:
        result["failure_class_only"]=type(exc).__name__
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",
                        encoding="utf-8")
    print(json.dumps({
        "technical_status":result["technical_status"],
        "scientific_status":result["scientific_status"],
        "n_real_original_coloc_R_runtime_tests_passed":
            result.get("n_original_coloc_R_checks_passed"),
        "cc_beta_SE_requires_s":result.get(
            "direct_beta_plus_varbeta_case_fraction_s_required_by_original_coloc"),
        "n_coloc_executed":0,
        "failure_class_only":result.get("failure_class_only")
    },ensure_ascii=False),flush=True)
    return 0 if result["technical_status"].startswith("PASS_") else 1


if __name__=="__main__":
    raise SystemExit(main())
