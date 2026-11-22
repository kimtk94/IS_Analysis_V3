#!/usr/bin/env python3
"""Verify genuine GIGASTROKE GRCh37 chr10 stroke-locus source GWAS rows.

Read the whole publisher gzip into ephemeral disk while computing full-file
publisher MD5 and pinned SHA256, then parse the original raw per-SNP beta/SE/
p/EAF/effect allele in a predeclared region. Emit ONLY group-level summaries;
no individual variant rows, no invented N/SE, no automatic coloc permissions.

All 36 P3A assay identities share one chr10 region. 36 assay rows therefore
must NOT be counted as 36 independent GWAS loci or replications.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import tempfile
import urllib.error
import urllib.request

GATE=importlib.import_module("scripts.is.p3a_publisher_md5_gate")
REGION={
    "id":"IS_XDATA_G0015",
    "chromosome":"10",
    "original_region_start":104663676,
    "original_region_end":105720351,
    "source_context_start":104163676,
    "source_context_end":106220351,
}
LEGACY_AIS_LEAD={
    "chromosome":"10","position":105596641,
    "legacy_alt_allele":"C","legacy_other_allele":"G",
    "legacy_ALT_beta":-0.0762,
    "legacy_p":6.214e-11,
    "previous_source_label":"GWS_EAS_AIS",
    "meaning":"single previously archived regional lead only; not a >20 SNP independent source match"
}
HEADER=("chromosome","base_pair_location","effect_allele","other_allele",
        "effect_allele_frequency","beta","standard_error","p_value")
MAX_COMPRESSED_BYTES=330_000_000
MAX_GWAS_ROWS=9_000_000

def require(ok, msg):
    if not ok:
        raise ValueError(msg)

def inspect_row(row):
    # Return scientific status, not raw rows, to avoid uploading SNP effects.
    chrom=row["chromosome"].strip().upper().removeprefix("CHR")
    try:
        pos=int(row["base_pair_location"])
    except (ValueError,TypeError):
        return None
    if chrom!=REGION["chromosome"] or not (
            REGION["source_context_start"]<=pos<=REGION["source_context_end"]):
        return None
    try:
        effect=float(row["beta"]); se=float(row["standard_error"])
        p=float(row["p_value"]); eaf=float(row["effect_allele_frequency"])
    except (ValueError,TypeError):
        return {"error":"NONNUMERIC","pos":pos}
    ea=row["effect_allele"].strip().upper()
    oa=row["other_allele"].strip().upper()
    if not (
        all(math.isfinite(v) for v in (effect,se,p,eaf))
        and se>0 and 0<=p<=1 and 0<=eaf<=1
        and ea in ("A","C","G","T") and oa in ("A","C","G","T") and ea!=oa
    ):
        return {"error":"INVALID_EFFECT_SE_ALLELES_OR_P","pos":pos}
    return {"pos":pos,"beta":effect,"se":se,"p":p,"eaf":eaf,
            "effect_allele":ea,"other_allele":oa,
            "palindrome":frozenset((ea,oa)) in (
                frozenset(("A","T")),frozenset(("C","G")))}

def compare_prior_lead(parsed):
    require(parsed["pos"]==LEGACY_AIS_LEAD["position"],"lead not same position")
    pair={parsed["effect_allele"],parsed["other_allele"]}
    if pair!={LEGACY_AIS_LEAD["legacy_alt_allele"],LEGACY_AIS_LEAD["legacy_other_allele"]}:
        return {"status":"ALLELES_DISCORDANT","sign_and_p_concordant":False}
    oriented=(parsed["beta"] if parsed["effect_allele"]==
              LEGACY_AIS_LEAD["legacy_alt_allele"] else -parsed["beta"])
    p=parsed["p"]
    p0=LEGACY_AIS_LEAD["legacy_p"]
    bp=LEGACY_AIS_LEAD["legacy_ALT_beta"]
    agree=(math.isclose(oriented,bp,rel_tol=.015,abs_tol=.00005)
            and p>0 and abs(math.log10(p)-math.log10(p0))<=.05)
    return {
        "status":"PRIOR_REGIONAL_LEAD_BETA_P_MATCH" if agree
                 else "PRIOR_REGIONAL_LEAD_NUMERIC_DISCORDANCE",
        "same_unordered_allele_pair":True,
        "signed_legacy_alt_beta_match":math.isclose(oriented,bp,rel_tol=.015,abs_tol=.00005),
        "p_value_log10_match":p>0 and abs(math.log10(p)-math.log10(p0))<=.05,
        "sign_and_p_concordant":agree,
        "n_independent_original_variant_matches_for_source_approval":1 if agree else 0,
        "source_identity_approved":False,
        "allele_ref_validated_with_independent_FASTA":False
    }

def scan_original(path,ancestry, min_total_rows=1_000_000, min_context_rows=100):
    count=0; valid=0; failed=0; pal=0; significant=0; underflow=0
    in_original=0; minp=None; lead=None; header=None
    positions=set();dupe=0
    with gzip.open(path,"rt",encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f,delimiter="\t")
        header=list(reader.fieldnames or [])
        require(set(HEADER).issubset(header),"original publisher GWAS missing beta/SE/alleles")
        for raw in reader:
            count+=1
            require(count<=MAX_GWAS_ROWS,"unbounded source SNP input")
            value=inspect_row(raw)
            if value is None:continue
            if "error" in value:
                failed+=1
                continue
            valid+=1
            in_original+=REGION["original_region_start"]<=value["pos"]<=REGION["original_region_end"]
            if value["pos"] in positions:dupe+=1
            positions.add(value["pos"])
            pal+=value["palindrome"]
            significant+=value["p"]<5e-8
            underflow+=value["p"]==0
            if minp is None or value["p"]<minp:minp=value["p"]
            if value["pos"]==LEGACY_AIS_LEAD["position"]:
                require(lead is None,"multiple original variants at lead position (ambiguity)")
                lead=compare_prior_lead(value)
    require(count>=min_total_rows,"not a complete GWAS variant file")
    require(valid>=min_context_rows,"insufficient original GIGASTROKE chr10 context coverage")
    return {
        "original_GWAS_rows_read":count,
        "original_GWAS_fields_present":list(HEADER),
        "source_context_GRCh37":REGION,
        "source_context_valid_variant_rows":valid,
        "original_region_valid_variant_rows":in_original,
        "source_context_failed_numeric_or_allele_rows":failed,
        "source_context_duplicate_positions":dupe,
        "source_context_palindromic_rows":pal,
        "source_context_p_lt_5e8_rows":significant,
        "source_context_p_zero_underflow_rows":underflow,
        "source_context_min_p_descriptive":minp,
        "legacy_GWS_EAS_AIS_lead_comparison":lead if ancestry=="EAS" else None,
        "source_GWAS_effect_and_SE_original_rows_checked":True,
        "source_GWAS_effects_reconciled_to_archived_fullcis_36_groups":False,
        "scientific_status":"HOLD"
    }

def checked_manifest(path):
    m=json.loads(path.read_text(encoding="utf-8"))
    groups=m.get("groups",[])
    require(m.get("schema")=="is-p3a-36-real-group-keys-from-legacy-393-v1"
            and m.get("scientific_status")=="HOLD"
            and m.get("origin",{}).get("commit_sha")==
                 "0bf646fd6f29322d17da878253f37def24c01c74"
            and m.get("region")=="IS_XDATA_G0015"
            and len(groups)==36
            and len({(z["region_id"],z["QTL_dataset"],z["gene_id"]) for z in groups})==36
            and len({z["gene_id"] for z in groups})==20
            and all(z["region_id"]=="IS_XDATA_G0015"
                    and z["approved_for_new_coloc"] is False
                    and z["original_P3A_source_stat_positive_proven"] is False
                    for z in groups),
            "not genuine recovered 36/20 group source gate")
    return groups

def stream_source(ancestry,dest):
    spec=GATE.SOURCES[ancestry]
    url,_=GATE.build_urls(spec)
    length=0;sha=hashlib.sha256();md5=hashlib.md5(usedforsecurity=False)
    req=urllib.request.Request(url,headers={
        "User-Agent":"IS-Master-Research-GIGASTROKE-Chr10-GeneGroup-Original/1.0"})
    with urllib.request.urlopen(req,timeout=80) as resp:
        require(resp.status==200,"GWAS original full source HTTP 200 required")
        with dest.open("wb") as f:
            while True:
                chunk=resp.read(4*1024*1024)
                if not chunk:break
                length+=len(chunk)
                require(length<=MAX_COMPRESSED_BYTES,"GWAS file size security bound exceeded")
                sha.update(chunk);md5.update(chunk);f.write(chunk)
    require(length==spec["bytes"] and sha.hexdigest()==spec["full_sha256_pinned"]
            and md5.hexdigest()==spec["publisher_md5_pinned"],
            "original published GWAS MD5/SHA256/size not matched to previous verified source")
    return {
        "accession":spec["accession"],
        "source_original_full_file_bytes":length,
        "publisher_md5_and_prior_sha256_verified":True,
        "source_release_build":"GRCh37",
        "source_effect_allele_orientation":"PUBLISHER_EFFECT_ALLELE_NOT_ASSUMED_ALT",
        "per_variant_effective_N_in_source":"NOT_PRESENT_DO_NOT_INFER"
    }

def build_group_rows(groups,whole,ancestry):
    require(whole["scientific_status"]=="HOLD"
            and whole["source_GWAS_effects_reconciled_to_archived_fullcis_36_groups"] is False,
            "GWAS source locus technical PASS cannot grant variant-wise coloc")
    return [{
        "region_id":x["region_id"],"QTL_dataset":x["QTL_dataset"],
        "gene_id":x["gene_id"],"gene_symbol":x["gene_symbol"],
        "GWAS_source_ancestry":ancestry,
        "GWAS_original_chr10_context_rows_QC":whole["source_context_valid_variant_rows"],
        "QTL_source_complete_cis_tested_variants_from_archived_BH05":
            x["official_full_cis_tested_variant_count"],
        "QTL_public_BH05_assay":x["official_studywide_eGene_BH05_assay"],
        "variant_level_QTL_GWAS_signed_effect_match_proven":False,
        "GWAS_exact_case_fraction_and_group_N_approved":False,
        "QTL_effect_SE_N_alleles_independently_approved":False,
        "source_cohort_matched_LD_approved":False,
        "new_coloc_or_replication_executed":False,
        "science_status":"HOLD_NO_APPROVED_SOURCE_PAIRED_SUMSTATS"
    } for x in groups]

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--ancestry",required=True,choices=("EAS","EUR"))
    ap.add_argument("--manifest",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    try:
        groups=checked_manifest(args.manifest)
        with tempfile.TemporaryDirectory(prefix="is_chr10_ephemeral_") as tmp:
            source=Path(tmp)/"original_gwas.tsv.gz"
            integrity=stream_source(args.ancestry,source)
            scanned=scan_original(source,args.ancestry)
        group_rows=build_group_rows(groups,scanned,args.ancestry)
        result={
            "schema":"is-p3a-chr10-original-publisher-gwas-all-36-v1",
            "technical_status":"PASS_PUBLISHER_WHOLE_FILE_AND_CHR10_ORIGINAL_ROWS",
            "scientific_status":"HOLD",
            "ancestry":args.ancestry,
            "source":integrity,"context":scanned,
            "actual_36_assay_groups_considered":len(group_rows),
            "unique_candidate_gene_ids":len({x["gene_id"] for x in group_rows}),
            "independent_GWAS_genomic_loci_replicated":0,
            "new_coloc_performed":0,
            "all_group_input_statuses":group_rows
        }
    except (OSError,ValueError,TypeError,KeyError,EOFError,
            urllib.error.URLError) as exc:
        result={
            "schema":"is-p3a-chr10-original-publisher-gwas-all-36-v1",
            "technical_status":"HOLD_ORIGINAL_SOURCE_OR_QC_ERROR",
            "scientific_status":"HOLD","ancestry":args.ancestry,
            "actual_36_assay_groups_considered":0,"new_coloc_performed":0,
            "error_type":type(exc).__name__,"error":str(exc)[:220]
        }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:result.get(k) for k in (
        "technical_status","scientific_status","ancestry",
        "actual_36_assay_groups_considered","new_coloc_performed","error")},ensure_ascii=False))
    return 0 if result["technical_status"].startswith("PASS_") else 1

if __name__=="__main__":
    raise SystemExit(main())
