#!/usr/bin/env python3
"""Primary GIANT Locke 2015 EUR BMI full-source vs pinned 2018 GIANT+UKB BMI panel.

- Streams official primary original gzip once; complete SHA-256, exact bytes, gzip CRC.
- Compares all original rsIDs and ordered alleles to pinned 1,007,899-SNP BMI panel.
- Reports only aggregate N/allele/P QC and source provenance; no sumstats exported.
- Pre-UKB publication is NOT proof of GIANT vs stroke participant independence.
- Does not run LDSC, rg, MR, or infer a cross-study genetic effect.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime,timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
from statistics import median
from urllib.parse import urlsplit
import urllib.request

from scripts.ckm_stroke.verify_bmi_ldsc_candidate import inspect
from scripts.ckm_stroke.audit_gigastroke_bmi_allele_harmonization import orientation,GOOD

SOURCE_URL="https://giant-consortium.web.broadinstitute.org/images/1/15/SNP_gwas_mc_merge_nogc.tbl.uniq.gz"
DOC_URL="https://giant-consortium.web.broadinstitute.org/images/b/b9/Data_File_Description_2012-2015_GIANT_GWAS_Summary_Statistics.docx"
SOURCE_SIZE=35_869_809
# Pinned from complete 2026-11-17 first-party GIANT original stream; NOT publisher-provided MD5.
EXPECTED_SOURCE_SHA256="ddca3b30985e3428f4762f8e8aadd385504627b9b825cff1ecdc4a5e2e358b44"
EXPECTED_ORIGINAL_ROWS=2_554_637
EXPECTED_DUPLICATE_EXTRA_ROWS=14
EXPECTED_HARMONIZABLE=1_007_858
BMI_PANEL_SHA="ed8f713827cd3281bea5f88a9cb42a3eabfecb1ee67c912ca7e86f151e4399bf"
PANEL_COUNT=1_007_899
MAX_ORIGINAL_ROWS=4_000_000
MAX_LINE_BYTES=1024
HEAD=("SNP","A1","A2","Freq1.Hapmap","b","se","p","N")
MIN_N=1
MAX_N=1_000_000
MAX_BYTES=40_000_000


class DigestReader:
    """Apply SHA256 and bounded byte count to raw compressed stream."""
    def __init__(self,response):
        self.response=response
        self.sha=hashlib.sha256()
        self.n=0

    def read(self,count=-1):
        requested=1_048_576 if count is None or count<0 else min(count,1_048_576)
        raw=self.response.read(requested)
        self.n+=len(raw)
        if self.n>MAX_BYTES:raise ValueError("source compression size cap exceeded")
        self.sha.update(raw)
        return raw


def panel_from_candidate(candidate:Path,prov:Path):
    previous=json.loads(prov.read_text(encoding="utf-8"))
    checked=inspect(candidate,previous)
    if checked["candidate_rows_verified"]!=PANEL_COUNT or checked["candidate_plaintext_sha256"]!=BMI_PANEL_SHA:
        raise ValueError("pinned 2018 BMI candidate differs from prior verified 1M panel")
    panel={}
    with gzip.open(candidate,"rt",encoding="ascii") as gz:
        if gz.readline().rstrip("\r\n")!="SNP\tN\tZ\tA1\tA2":
            raise ValueError("candidate columns changed")
        for line in gz:
            parts=line.rstrip("\r\n").split("\t")
            if len(parts)!=5:raise ValueError("candidate malformed")
            snp,n,z,a1,a2=parts
            if snp in panel:raise ValueError("candidate duplicate")
            panel[snp]=(a1,a2)
    if len(panel)!=PANEL_COUNT:raise ValueError("BMI candidate row count changed")
    return panel


def compare_source(lines,panel,*,max_rows=MAX_ORIGINAL_ROWS):
    counts=Counter()
    original_seen=set()
    matched_seen=set()
    ambiguous_duplicate_rsids=set()
    candidate_record_by_snp={}
    match_n=[]
    n_min=None
    n_max=None
    first=next(lines,None)
    if first is None or tuple(first.rstrip("\r\n").split("\t"))!=HEAD:
        raise ValueError("official Locke 2015 original header changed")
    for num,line in enumerate(lines,start=1):
        if num>max_rows:raise ValueError("original BMI rows exceeded cap")
        counts["original_rows"]+=1
        parts=line.rstrip("\r\n").split("\t")
        if len(parts)!=8:raise ValueError("original row width invalid")
        snp,a1,a2,hapmap_freq,beta,se,p,n=parts
        if not snp:
            raise ValueError("original has missing SNP identifier")
        if snp in original_seen:
            counts["original_duplicate_rsid_rows"]+=1
            if snp not in ambiguous_duplicate_rsids:
                ambiguous_duplicate_rsids.add(snp)
                if snp in panel:
                    counts["candidate_duplicate_rsid_excluded"]+=1
            candidate_record_by_snp.pop(snp,None)
            continue
        original_seen.add(snp)
        counts["original_unique_rsids"]+=1
        if snp in panel:
            matched_seen.add(snp)
        try:
            nval=float(n)
        except ValueError:
            nval=float("nan")
        if not math.isfinite(nval) or not MIN_N<=nval<=MAX_N:
            counts["original_variant_N_missing_or_invalid"]+=1
            if snp in panel:counts["candidate_match_with_invalid_N"]+=1
            continue
        if not nval.is_integer():
            counts["original_fractional_N_rows"]+=1
            if snp in panel:counts["candidate_fractional_N_rows"]+=1
        n_min=nval if n_min is None else min(n_min,nval)
        n_max=nval if n_max is None else max(n_max,nval)
        if a1 not in "ACGT" or a2 not in "ACGT" or len(a1)!=1 or len(a2)!=1 or a1==a2:
            counts["original_invalid_alleles"]+=1
            if snp in panel:counts["candidate_match_with_invalid_original_alleles"]+=1
            continue
        try:
            nb=float(beta);ns=float(se);pv=float(p)
        except ValueError:
            counts["original_invalid_numeric"]+=1
            if snp in panel:counts["candidate_match_with_invalid_original_numeric"]+=1
            continue
        if not(all(map(math.isfinite,(nb,ns,pv))) and ns>0 and 0<=pv<=1):
            counts["original_invalid_numeric"]+=1
            if snp in panel:counts["candidate_match_with_invalid_original_numeric"]+=1
            continue
        # HapMap CEU frequency is *reference-only*, not exposure cohort EAF.
        # In the official primary data it can be 'NA'. Never use this field to
        # reject an otherwise well-defined BETA/SE/N association.
        try:
            hf=float(hapmap_freq)
        except ValueError:
            hf=float("nan")
        if not math.isfinite(hf) or not 0<=hf<=1:
            counts["original_hapmap_CEU_frequency_missing_or_invalid"]+=1
            if snp in panel:
                counts["candidate_hapmap_CEU_frequency_missing_or_invalid"]+=1
        counts["valid_original_variants"]+=1
        if snp not in panel:
            counts["original_not_in_2018_candidate"]+=1
            continue
        label=orientation(panel[snp][0],panel[snp][1],a1,a2)
        if snp not in ambiguous_duplicate_rsids:
            candidate_record_by_snp[snp]=(label,nval,pv)
    if not counts["original_rows"]:
        raise ValueError("empty original")
    if counts["original_unique_rsids"]+counts["original_duplicate_rsid_rows"]!=counts["original_rows"]:
        raise ValueError("full original duplicate accounting differs")
    counts["original_distinct_ambiguous_duplicate_rsids"]=len(ambiguous_duplicate_rsids)
    for snp,(label,nval,pv) in candidate_record_by_snp.items():
        counts["candidate_"+label]+=1
        if label in GOOD:
            counts["candidate_harmonizable"]+=1
            match_n.append(nval)
            if pv<=5e-8:counts["candidate_original_p_le_5e8"]+=1
            if pv==0:counts["candidate_original_p_zero"]+=1
    counts["candidate_absent_from_2015_original"]=len(panel)-len(matched_seen)
    counts["candidate_rsids_present_in_2015_original"]=len(matched_seen)
    counts["candidate_allele_or_numeric_rejected"]=len(matched_seen)-counts["candidate_harmonizable"]
    counts["candidate_harmonizable_missing_from_full_panel"]=len(panel)-counts["candidate_harmonizable"]
    if not match_n:raise ValueError("2015 original has no harmonizable candidate SNPs")
    match_n.sort()
    def q(frac):
        return match_n[int((len(match_n)-1)*frac)]
    return {
        "counts":dict(counts),
        "all_original_variant_N_min":n_min,
        "all_original_variant_N_max":n_max,
        "harmonizable_variant_N_summary":{
            "n":len(match_n),"min":match_n[0],"p10":q(.1),
            "median":median(match_n),"p90":q(.9),"max":match_n[-1],
        },
        "candidate_total":len(panel),
        "candidate_overlap_pct":round(100*counts["candidate_harmonizable"]/len(panel),4),
    }


def audit(candidate:Path,prov:Path,fetcher=urllib.request.urlopen):
    report={
        "schema":"ckm-locke-2015-eur-bmi-ukb-avoidance-candidate-source-v1",
        "checked_utc":datetime.now(timezone.utc).isoformat(),
        "source":SOURCE_URL,"publisher_document":DOC_URL,
        "original_pub_year":2015,
        "population":"European GIANT meta-analysis (Locke et al. 2015)",
        "ukb_status":"PRE_UKB_META_ANALYSIS; NOT PROOF OF ZERO OTHER SHARED COHORTS",
        "source_original_publisher_sha256_available":False,
        "scientific_status":"HOLD","ready_for_ldsc":False,
        "LDSC_h2_rg_MR_executed":False,"raw_source_persisted":False,
        "technical_status":"HOLD_SOURCE_OR_MATCH_QC"
    }
    try:
        panel=panel_from_candidate(candidate,prov)
        req=urllib.request.Request(SOURCE_URL,headers={"User-Agent":"CKM-Stroke-Locke2015-Source-Audit/1.0"})
        with fetcher(req,timeout=60) as resp:
            u=urlsplit(resp.geturl())
            if resp.status!=200 or u.scheme!="https" or u.hostname!="giant-consortium.web.broadinstitute.org":
                raise ValueError("official publisher source URL changed")
            length=resp.headers.get("Content-Length","")
            if not length.isdecimal() or int(length)!=SOURCE_SIZE:
                raise ValueError("publisher source compressed size changed")
            reader=DigestReader(resp)
            with gzip.GzipFile(fileobj=reader,mode="rb") as gz:
                def decode():
                    while True:
                        raw=gz.readline(MAX_LINE_BYTES+1)
                        if not raw:break
                        if len(raw)>MAX_LINE_BYTES:raise ValueError("oversize original line")
                        yield raw.decode("ascii")
                result=compare_source(decode(),panel)
            while reader.read(1_048_576):pass
            if reader.n!=SOURCE_SIZE:
                raise ValueError("truncated original compressed source")
            if reader.sha.hexdigest()!=EXPECTED_SOURCE_SHA256:
                raise ValueError("original publisher source differs from pinned complete SHA256 baseline")
            if (result["counts"]["original_rows"]!=EXPECTED_ORIGINAL_ROWS or
                result["counts"]["original_duplicate_rsid_rows"]!=EXPECTED_DUPLICATE_EXTRA_ROWS or
                result["counts"]["candidate_harmonizable"]!=EXPECTED_HARMONIZABLE):
                raise ValueError("original full-row / duplicate / candidate overlap pin changed")
        report.update({
            "technical_status":"PASS_ORIGINAL_SOURCE_SHA256_AND_BMI_PANEL_ALLELE_QC_ONLY",
            "original_compressed_bytes":reader.n,
            "original_compressed_sha256":reader.sha.hexdigest(),
            "original_gzip_integrity":"PASS",
            **result
        })
    except Exception as exc:
        report["technical_status"]="HOLD_SOURCE_OR_MATCH_QC"
        report["error_type"]=type(exc).__name__
        report["error_detail"]=str(exc)[:240]
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--candidate",type=Path,required=True)
    p.add_argument("--provenance",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    out=audit(args.candidate,args.provenance)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":out["technical_status"],
                     "source_sha256":out.get("original_compressed_sha256"),
                     "counts":out.get("counts"),
                     "matched_N":out.get("harmonizable_variant_N_summary"),
                     "error":out.get("error_detail"),
                     "scientific_status":out["scientific_status"]},ensure_ascii=False))
    if out["technical_status"]!="PASS_ORIGINAL_SOURCE_SHA256_AND_BMI_PANEL_ALLELE_QC_ONLY":
        raise SystemExit(1)


if __name__=="__main__":
    main()
