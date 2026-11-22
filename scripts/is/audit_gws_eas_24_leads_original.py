#!/usr/bin/env python3
"""Original-source 24 NONPAL lead-variant reproduction across 30 stroke regions.

Historical source: pinned 2026-10-09 180-row region×GWAS table (GWS_EAS_AIS).
Original source: publisher GCST90104545 GRCh37 FULL file, SHA256+MD5 verified.

Read 24 real nonpal lead keys across multiple chromosomes; compare same-position,
allele pair, signed ALT beta, p, EAF to original effect allele. Four A/T or
C/G palindromes are withheld; two missing source leads are recorded missing.

Both old and new may derive from the same consortium GWAS! This demonstrates
file-to-archive reproduction, NOT independent external phenotype replication.
Never infer case fraction, per-SNP N or any P3A 36-group QTL coloc approval.
"""
from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import importlib
import io
import json
import math
from pathlib import Path
import re
import tempfile
import urllib.error
import urllib.request

GATE=importlib.import_module("scripts.is.p3a_publisher_md5_gate")
EXPECTED_LEGACY_BLOB="b34b4ae84ea4e5b0c72e3982d530a048cc9be5ca"
EXPECTED_LEGACY_COMMIT="0bf646fd6f29322d17da878253f37def24c01c74"
EXPECTED_FAMILY="GWS_EAS_AIS"
HEADER=("chromosome","base_pair_location","effect_allele","other_allele",
        "beta","standard_error","p_value","effect_allele_frequency")
PAIR=re.compile(r"^([1-9]|1[0-9]|2[0-2]):([1-9]\d*):([ACGT]):([ACGT])$")
PALIN={frozenset(("A","T")),frozenset(("C","G"))}
BETA_ABS_TOL=0.00006
BETA_REL_TOL=0.015
P_LOG_TOL=0.05
AF_ABS_TOL=0.025

def ensure(condition, reason):
    if not condition:
        raise ValueError(reason)

def read_legacy_180(path):
    raw=path.read_bytes()
    digest=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    ensure(digest==EXPECTED_LEGACY_BLOB,"original historical GWAS region table Git blob SHA drift")
    rows=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig")),delimiter="\t"))
    ensure(len(rows)==180 and all(None not in row for row in rows),
           "historical 180 source-region-table rows missing or malformed")
    ensure(len({(x["region_id"],x["source_dataset"]) for x in rows})==180,
           "duplicate original 30×6 region GWAS source key")
    eas=[x for x in rows if x["source_dataset"]==EXPECTED_FAMILY]
    ensure(len(eas)==30 and len({x["region_id"] for x in eas})==30,
           "original GWS_EAS_AIS 30 region source not complete")
    selected=[]; withheld=[]; missing=[]
    for x in eas:
        region=x["region_id"]
        ensure(re.fullmatch(r"IS_XDATA_G\d{4}",region) is not None,
               "unexpected genomic group identity")
        raw_key=x["group_source_lead_variant"].strip()
        if not raw_key:
            ensure(x["group_source_lead_allele_status"]=="NOT_PRESENT",
                   "missing lead status should be explicit")
            missing.append({"region_id":region,"status":"ARCHIVED_LEAD_NOT_PRESENT"})
            continue
        match=PAIR.fullmatch(raw_key)
        ensure(match is not None,"invalid GRCh37 original lead key: "+region)
        chrom,pos,ref,alt=match.groups()
        ensure(ref!=alt,"REF=ALT forbidden")
        ensure(x["chromosome_GRCh37"]==chrom
               and int(x["context_window_start"])<=int(pos)<=int(x["context_window_end"]),
               "source lead outside stated GRCh37 chromosome/window")
        beta=float(x["group_source_lead_beta_ALT"])
        p=float(x["group_source_lead_p"])
        eaf=float(x["group_source_lead_ALT_AF"])
        ensure(all(math.isfinite(y) for y in (beta,p,eaf)) and
               0<p<=1 and 0<=eaf<=1,"missing numeric effect/p/EAF: "+region)
        record={"region_id":region,"GRCh37_variant_key":raw_key,
                "chromosome":chrom,"position":int(pos),
                "historical_ref_allele":ref,"historical_alt_allele":alt,
                "historical_beta_ALT":beta,"historical_p":p,"historical_ALT_EAF":eaf,
                "historical_allele_flag":x["group_source_lead_allele_status"]}
        if frozenset((ref,alt)) in PALIN:
            withheld.append({"region_id":region,"GRCh37_variant_key":raw_key,
                             "status":"PALINDROMIC_STRAND_HOLD"})
        else:
            selected.append(record)
    ensure(len(selected)==24 and len(withheld)==4 and len(missing)==2,
           "original 24 nonpal / 4 palindromic / 2 absent source breakdown drift")
    ensure(len({(r["chromosome"],r["position"]) for r in selected})==24
           and len({r["chromosome"] for r in selected})>=10,
           "24 variants are not positioned across distinct chromosomes")
    return {"selected":selected,"withheld":withheld,"missing":missing,
            "historical_rows":180,"source_family":EXPECTED_FAMILY}

def inspect_original_row(row):
    try:
        chrom=str(row["chromosome"]).strip().upper().removeprefix("CHR")
        pos=int(row["base_pair_location"])
        beta=float(row["beta"])
        se=float(row["standard_error"])
        p=float(row["p_value"])
        af=float(row["effect_allele_frequency"])
        ea=row["effect_allele"].strip().upper()
        oa=row["other_allele"].strip().upper()
    except (KeyError,ValueError,TypeError):
        return None
    if not (chrom in [str(i) for i in range(1,23)] and pos>0
            and all(math.isfinite(v) for v in (beta,se,p,af))
            and se>0 and 0<=p<=1 and 0<=af<=1
            and ea in "ACGT" and oa in "ACGT" and ea!=oa):
        return None
    return {"chrom":chrom,"pos":pos,"beta":beta,"se":se,"p":p,
            "eaf":af,"effect_allele":ea,"other_allele":oa}

def compare(record,original_rows):
    if not original_rows:
        return {"region_id":record["region_id"],"GRCh37_variant_key":record["GRCh37_variant_key"],
                "status":"OFFICIAL_LEAD_POSITION_ABSENT",
                "position_covered":False,"nonpalindromic":True,
                "joint_effect_p_eaf_reproduced":False}
    compatible=[r for r in original_rows
                if {r["effect_allele"],r["other_allele"]}==
                   {record["historical_ref_allele"],record["historical_alt_allele"]}]
    if len(compatible)!=1 or len(original_rows)!=1:
        return {"region_id":record["region_id"],"GRCh37_variant_key":record["GRCh37_variant_key"],
                "status":"AMBIGUOUS_OR_INCOMPATIBLE_ORIGINAL_ALLELES_AT_POSITION",
                "position_covered":True,"original_rows_at_position":len(original_rows),
                "allele_pair_matches":len(compatible),"nonpalindromic":True,
                "joint_effect_p_eaf_reproduced":False}
    raw=compatible[0]
    eff_alt=(raw["effect_allele"]==record["historical_alt_allele"])
    oriented_beta=raw["beta"] if eff_alt else -raw["beta"]
    oriented_eaf=raw["eaf"] if eff_alt else 1.0-raw["eaf"]
    beta_pass=math.isclose(oriented_beta,record["historical_beta_ALT"],
                           rel_tol=BETA_REL_TOL,abs_tol=BETA_ABS_TOL)
    p_pass=raw["p"]>0 and abs(math.log10(raw["p"])-math.log10(
        record["historical_p"]))<=P_LOG_TOL
    eaf_diff=abs(oriented_eaf-record["historical_ALT_EAF"])
    af_pass=eaf_diff<=AF_ABS_TOL
    passed=bool(beta_pass and p_pass and af_pass)
    return {
        "region_id":record["region_id"],
        "GRCh37_variant_key":record["GRCh37_variant_key"],
        "chromosome":record["chromosome"],
        "nonpalindromic":True,"position_covered":True,
        "original_rows_at_position":len(original_rows),
        "allele_pair_matches":1,
        "effect_allele_was_historical_ALT":eff_alt,
        "ALT_oriented_beta_match":beta_pass,
        "p_value_log10_tolerance_match":p_pass,
        "ALT_oriented_eaf_match":af_pass,
        "absolute_eaf_difference":round(eaf_diff,7),
        "original_effect_se_finite_and_positive":raw["se"]>0,
        "joint_effect_p_eaf_reproduced":passed,
        "status":"SIGNED_ORIGINAL_EAS_LEAD_REPRODUCED_ARCHIVE" if passed
                 else "ORIGINAL_EAS_LEAD_NUMERIC_DISCORDANCE_HOLD",
        "original_study_case_fraction_from_regional_row_approved":False,
        "original_GWAS_SNP_effect_valid_for_P3A_coloc":False
    }

def scan_original_file(path,records,min_rows=1_000_000):
    requested={(x["chromosome"],x["position"]) for x in records["selected"]}
    collected=collections.defaultdict(list)
    count=0;qualified=0
    with gzip.open(path,"rt",encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f,delimiter="\t")
        ensure(set(HEADER).issubset(reader.fieldnames or []),
               "official GWAS original absent mandatory signed effect/SE/p/EAF")
        for row in reader:
            count+=1
            ensure(count<=9_000_000,"unexpected unbounded original SNP rows")
            chrom=str(row["chromosome"]).strip().upper().removeprefix("CHR")
            key=(chrom,row["base_pair_location"].strip())
            # Keep numeric keys; no per-SNP output of 6M unrequested variant effects.
            try:key=(chrom,int(key[1]))
            except ValueError:continue
            if key not in requested:continue
            parsed=inspect_original_row(row)
            if parsed is not None:
                collected[key].append(parsed)
                qualified+=1
    ensure(count>=min_rows and qualified>0,"full original EAS summary input not scanned")
    results=[compare(a,collected[(a["chromosome"],a["position"])])
             for a in records["selected"]]
    n_pass=sum(x["joint_effect_p_eaf_reproduced"] for x in results)
    n_pair=sum(x.get("allele_pair_matches")==1 for x in results)
    n_chrom=len({x["chromosome"] for x in results if
                 x.get("joint_effect_p_eaf_reproduced")})
    high_support=(n_pass>=23 and n_chrom>=10 and n_pair>=23)
    return {
        "schema":"is-gws-eas-30-provisional-region-24-nonpal-original-variant-audit-v1",
        "original_GWAS_full_variant_rows_scanned":count,
        "original_source_nominal_build":"GRCh37",
        "historical_GWAS_table_blob_sha1":EXPECTED_LEGACY_BLOB,
        "actual_original_header_fields_present":list(HEADER),
        "historical_region_rows":30,
        "independently_excluded_palindromic_source_leads":len(records["withheld"]),
        "source_missing_lead_regions":len(records["missing"]),
        "original_nonpalindromic_archived_leads_tested":len(results),
        "unique_chromosomes_with_concordant_leads":n_chrom,
        "original_same_position_allele_pair_match_count":n_pair,
        "original_signed_beta_p_eaf_all_three_match_count":n_pass,
        "original_signed_beta_p_eaf_mismatch_count":24-n_pass,
        "technical_source_archive_signed_reproduction_screen":
            "PASS_REVIEW_REQUIRED" if high_support else "HOLD_SOURCE_REPRODUCTION_INSUFFICIENT",
        "falsely_claimed_independent_stroke_replication":False,
        "source_provenance_technical_only":True,
        "published_original_study_Ncase_Ncontrol_independently_approved":False,
        "BBJ_stroke_source_provenance_verified":False,
        "matched_QTL_cis_effects_verified":False,
        "coloc_abf_executed":0,
        "coloc_susie_executed":0,
        "scientific_status":"HOLD",
        "tested_leads":results,
        "palindromic_leads_excluded":records["withheld"],
        "historical_missing_leads":records["missing"],
        "caveat":"All 24 markers are provisional different regional leads, not independent GWAS replications, genome-wide LD-pruned variants, clinical sample-size approvals or causal-gene evidence"
    }

def remote_original_file(destination):
    spec=GATE.SOURCES["EAS"]
    url,_=GATE.build_urls(spec)
    req=urllib.request.Request(url,headers={
        "User-Agent":"IS-Master-Research-24-Nonpal-OfficialGwasAudit/1.0"})
    sha=hashlib.sha256()
    md5=hashlib.md5(usedforsecurity=False)
    bytes_seen=0
    with urllib.request.urlopen(req,timeout=80) as response:
        ensure(response.status==200,"publisher source did not respond with HTTP 200")
        with destination.open("wb") as f:
            while True:
                chunk=response.read(4*1024*1024)
                if not chunk:break
                bytes_seen+=len(chunk)
                ensure(bytes_seen<=spec["bytes"],"original GWAS exceeds known file size")
                sha.update(chunk);md5.update(chunk);f.write(chunk)
    ensure(bytes_seen==spec["bytes"] and
           sha.hexdigest()==spec["full_sha256_pinned"] and
           md5.hexdigest()==spec["publisher_md5_pinned"],
           "official source compressed original bytes, SHA256 or publisher MD5 changed")
    return {"accession":spec["accession"],"compressed_original_file_bytes":bytes_seen,
            "official_publisher_md5_and_pinned_sha256":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-legacy",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--input-official-gwas",type=Path,
                   help="Only for authorized preverified local full original source")
    args=p.parse_args()
    try:
        prior=read_legacy_180(args.original_legacy)
        with tempfile.TemporaryDirectory(prefix="is_gwas_24_orig_") as temp:
            if args.input_official_gwas is None:
                path=Path(temp)/"official_GCST90104545_buildGRCh37.tsv.gz"
                source=remote_original_file(path)
            else:
                path=args.input_official_gwas
                # Local file paths without previously audited publisher digest
                # MUST undergo the same complete integrity check.
                spec=GATE.SOURCES["EAS"]
                sha=hashlib.sha256()
                md5=hashlib.md5(usedforsecurity=False)
                size=0
                with path.open("rb") as f:
                    while True:
                        block=f.read(4*1024*1024)
                        if not block:break
                        size+=len(block)
                        ensure(size<=spec["bytes"],"local source file too big")
                        sha.update(block);md5.update(block)
                ensure(size==spec["bytes"] and sha.hexdigest()==spec["full_sha256_pinned"] and
                       md5.hexdigest()==spec["publisher_md5_pinned"],
                       "provided local GWAS source not publisher SHA/MD5 approved")
                source={"accession":spec["accession"],"compressed_original_file_bytes":size,
                        "official_publisher_md5_and_pinned_sha256":True}
            result=scan_original_file(path,prior)
        result["verified_whole_file_original_source"]=source
        result["technical_status"]="PASS_FULL_PUBLISHER_FILE_AND_24_LEAD_SOURCE_SCREEN"
    except (OSError,KeyError,ValueError,TypeError,EOFError,urllib.error.URLError) as error:
        result={"schema":"is-gws-eas-30-provisional-region-24-nonpal-original-variant-audit-v1",
                "technical_status":"HOLD_SOURCE_MISSING_OR_INPUT_MISMATCH",
                "scientific_status":"HOLD","coloc_abf_executed":0,"coloc_susie_executed":0,
                "error_type":type(error).__name__,"error":str(error)[:300]}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:result.get(k) for k in (
        "technical_status","scientific_status",
        "original_nonpalindromic_archived_leads_tested",
        "original_signed_beta_p_eaf_all_three_match_count",
        "original_same_position_allele_pair_match_count",
        "unique_chromosomes_with_concordant_leads",
        "technical_source_archive_signed_reproduction_screen","error")},ensure_ascii=False))
    return 0 if result["technical_status"].startswith("PASS_") else 1

if __name__=="__main__":
    raise SystemExit(main())
