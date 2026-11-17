#!/usr/bin/env python3
"""Bounded original publisher prefix: check beta/log(OR) scaling for four EUR outcomes.

No full source download. No rsID, N, signed LDSC effects or statistical estimates.
Primary-source full-MD5 independent verification is in the existing full GWAS audit.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime,timezone
import json
import math
from pathlib import Path
from statistics import median
import urllib.request
from urllib.parse import urlsplit
import zlib

from scripts.ckm_stroke.verify_gigastroke_raw import source_urls,load_manifest
from scripts.ckm_stroke.audit_gigastroke_bmi_allele_harmonization import COLUMNS

READ_BYTES=131072
MAX_DECOMPRESSED=1_000_000
MIN_ROWS=500
MAX_ABS_DELTA=1e-10
EXPECTED={
    "AIS":"GCST90104540",
    "LAA":"GCST90104542",
    "CE":"GCST90104541",
    "SVD":"GCST90104543",
}


def validate_prefix(phenotype:str,accession:str,compressed:bytes):
    if phenotype not in EXPECTED or accession!=EXPECTED[phenotype]:
        raise ValueError("unexpected phenotype/accession or LAA/CE swap")
    if len(compressed)!=READ_BYTES or compressed[:2]!=b"\x1f\x8b":
        raise ValueError("publisher gzip prefix bytes invalid")
    decompressor=zlib.decompressobj(16+zlib.MAX_WBITS)
    uncompressed=decompressor.decompress(compressed,MAX_DECOMPRESSED)
    if len(uncompressed)>=MAX_DECOMPRESSED:
        raise ValueError("unexpectedly highly compressible source prefix")
    # A Range prefix ends midway through the last line: discard the untrusted tail.
    if uncompressed.count(b"\n")<MIN_ROWS:
        raise ValueError("not enough complete original rows in compressed prefix")
    completed=uncompressed.split(b"\n")[:-1]
    header=completed.pop(0).decode("utf-8-sig").strip("\r").split("\t")
    if tuple(header)!=COLUMNS:
        raise ValueError("original outcome EUR effect schema changed")
    if any(key in header for key in ("rsid","N","n_samples","samplesize","N_eff")):
        raise ValueError("unexpected per-SNP N or rsID in original source schema")
    indices={key:header.index(key) for key in (
        "beta","standard_error","odds_ratio","p_value","effect_allele_frequency",
        "ci_lower","ci_upper","effect_allele","other_allele")}
    observations=[]
    stats={"complete_rows_examined":len(completed),"parsed_effect_rows":0,
           "p_value_zero_count":0,"nonfinite_or_invalid":0,
           "CI_interval_not_containing_OR":0}
    for line in completed:
        cols=line.strip(b"\r").decode("ascii").split("\t")
        if len(cols)!=len(COLUMNS):
            raise ValueError("unexpected row width in complete source line")
        try:
            beta,se,odds,p,eaf,ci_low,ci_high=[
                float(cols[indices[k]]) for k in (
                    "beta","standard_error","odds_ratio","p_value",
                    "effect_allele_frequency","ci_lower","ci_upper")]
        except ValueError as exc:
            raise ValueError("non-numeric effect in first-party sample") from exc
        if not(all(map(math.isfinite,(beta,se,odds,p,eaf,ci_low,ci_high))) and
               se>0 and odds>0 and ci_low>0 and ci_high>0 and ci_low<=odds<=ci_high and
               0<=p<=1 and 0<=eaf<=1):
            stats["nonfinite_or_invalid"]+=1
            continue
        if any(len(cols[indices[k]])!=1 or cols[indices[k]] not in "ACGT"
               for k in ("effect_allele","other_allele")):
            raise ValueError("unexpected original nucleotide alleles")
        delta=abs(beta-math.log(odds))
        observations.append(delta)
        stats["parsed_effect_rows"]+=1
        if p==0:stats["p_value_zero_count"]+=1
        if delta>MAX_ABS_DELTA:
            raise ValueError("beta != log(OR) in original source prefix")
    if stats["parsed_effect_rows"]<MIN_ROWS or stats["nonfinite_or_invalid"]:
        raise ValueError("original prefix invalid source numeric rows")
    return {
        "phenotype":phenotype,"accession":accession,
        "source_header_columns":list(COLUMNS),
        "outcome_beta_scale":"LOG_ODDS_ratio_approx_beta",
        "effect_allele":"source_effect_allele",
        "source_per_variant_N_available":False,
        "source_rsid_available":False,
        "sample_prefix_compressed_bytes":len(compressed),
        "observed":stats,
        "maximum_abs_beta_minus_lnOR":max(observations),
        "median_abs_beta_minus_lnOR":median(observations),
        "bounded_sample_not_genomewide_proof":True,
        "derived_Z_is_beta_over_se_if_N_then_approved":False,
    }


def fetch_prefix(row,opener=urllib.request.urlopen):
    raw,_,_=source_urls(row)
    req=urllib.request.Request(raw,headers={
        "User-Agent":"CKM-Stroke-GIGASTROKE-logOR-bounded-scale-check/1.0",
        "Range":f"bytes=0-{READ_BYTES-1}",
    })
    with opener(req,timeout=60) as r:
        uri=urlsplit(r.geturl())
        if (r.status!=206 or uri.scheme!="https" or uri.hostname!="ftp.ebi.ac.uk" or
            r.headers.get("Content-Range","").split(" ",1)[0].lower()!="bytes"):
            raise ValueError("publisher HTTP byte-range not honored")
        data=r.read(READ_BYTES+1)
    if len(data)!=READ_BYTES:
        raise ValueError("publisher prefix truncated or excessive")
    return validate_prefix(row["phenotype"],row["accession"],data)


def audit(manifest:Path,opener=urllib.request.urlopen):
    result={
        "schema":"ckm-gigastroke-eur-bounded-original-logor-scale-v1",
        "checked_utc":datetime.now(timezone.utc).isoformat(),
        "technical_status":"HOLD_ORIGINAL_OUTCOME_EFFECT_SCALE_PREFIX",
        "scientific_status":"HOLD","ready_for_ldsc":False,
        "LDSC_h2_intercept_rg_executed":False,
        "full_source_verified_in_this_run":False,
        "prefixes":[]
    }
    try:
        records=load_manifest(manifest)
        if {x["phenotype"] for x in records}!=set(EXPECTED):
            raise ValueError("original EUR manifest phenotype schema missing")
        for row in records:
            result["prefixes"].append(fetch_prefix(row,opener))
        result["technical_status"]="PASS_FOUR_BOUNDED_LOGOR_EFFECT_SCALE_ONLY"
    except Exception as exc:
        result["error_type"]=type(exc).__name__
        result["error_detail"]=str(exc)[:250]
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest",type=Path,default=Path("config/ckm_stroke/gigastroke_eur_official.tsv"))
    parser.add_argument("--out",type=Path,required=True)
    a=parser.parse_args()
    result=audit(a.manifest)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"technical_status":result["technical_status"],
                      "sampled":{r["phenotype"]:r["observed"]["parsed_effect_rows"] for r in result["prefixes"]},
                      "scientific":result["scientific_status"],
                      "error":result.get("error_detail")},ensure_ascii=False))
    if result["technical_status"]!="PASS_FOUR_BOUNDED_LOGOR_EFFECT_SCALE_ONLY":
        raise SystemExit(1)


if __name__=="__main__":
    main()
