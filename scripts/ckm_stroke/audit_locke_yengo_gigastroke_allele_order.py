#!/usr/bin/env python3
"""Cross-check original EUR GIGASTROKE effect allele order against independently
selected Locke 2015 and Yengo 2018 BMI source A1/A2. Sampled GWAS only.
No LDSC, no approved case/control N_eff and no participant-independence claims.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
from urllib.parse import urlsplit
import urllib.request
import zlib

from scripts.ckm_stroke.audit_gigastroke_bmi_allele_harmonization import (
    build_bmi_coordinate_panel, orientation, COLUMNS, GOOD
)
from scripts.ckm_stroke.audit_gigastroke_logor_prefix import (
    READ_BYTES, EXPECTED as ACCESSIONS, validate_prefix
)
from scripts.ckm_stroke.audit_ldsc_preregression_readiness import (
    EXPECTED as LOCKE_PROFILE, reference_evidence_contract
)
from scripts.ckm_stroke.verify_gigastroke_raw import load_manifest, source_urls

EXPECTED_ROWS = {"AIS": 3413, "LAA": 3256, "CE": 3273, "SVD": 3259}
MIN_SHARED_BMI = 500000
MIN_PREFIX_MATCH = 100

def original_locke_2015_panel(path, provenance):
    record=json.loads(provenance.read_text("utf-8"))
    expected=LOCKE_PROFILE["locke_2015"]
    if (record.get("technical_status")!=expected["status"]
        or record.get("scientific_status")!="HOLD"
        or record.get("ready_for_ldsc") is not False
        or record.get("2018_exposure_candidates_used") is not False
        or record.get("counts",{}).get("selected_rows")!=expected["rows"]
        or record.get("2015_candidate_plaintext_sha256")!=expected["sha256"]):
        raise ValueError("independent Locke 2015 provenance not pinned")
    reference_evidence_contract("locke_2015",record)
    digest=hashlib.sha256()
    alleles={}
    fractional=0
    with gzip.open(path,"rb") as src:
        header=src.readline(4097)
        if header != b"SNP\tN_REPORTED\tZ_FROM_BETA_SE\tA1\tA2\n":
            raise ValueError("Locke 2015 candidate schema mismatch")
        digest.update(header)
        for raw in src:
            if len(raw)>4096:
                raise ValueError("Locke candidate row oversized")
            digest.update(raw)
            row=raw.rstrip(b"\r\n").decode("ascii","strict").split("\t")
            if len(row)!=5:
                raise ValueError("Locke candidate row malformed")
            snp,n,z,a1,a2=row
            if (not snp or snp in alleles or a1 not in ("A","C","G","T")
                or a2 not in ("A","C","G","T") or a1==a2
                or frozenset((a1,a2)) in (frozenset("AT"),frozenset("CG"))):
                raise ValueError("invalid Locke candidate ID or allele pair")
            nvalue=float(n);zvalue=float(z)
            if not math.isfinite(nvalue) or nvalue<=0 or not math.isfinite(zvalue):
                raise ValueError("nonfinite independent BMI N/Z")
            if not nvalue.is_integer():
                fractional+=1
            alleles[snp]=(a1,a2)
            if len(alleles)>1_500_000:
                raise ValueError("Locke candidate count cap exceeded")
    if (digest.hexdigest()!=expected["sha256"] or len(alleles)!=expected["rows"]
        or fractional!=11_506):
        raise ValueError("independent Locke original SHA, row total or fractional N drift")
    return alleles,{"candidate_rows":len(alleles),
                    "plain_sha256":digest.hexdigest(),
                    "fractional_N_rows":fractional,
                    "selected_without_Yengo_2018":True}

def compare_source_orders(locke, yengo_coordinates,*,strict=True):
    counts=Counter()
    for rsid,a1,a2 in yengo_coordinates.values():
        other=locke.get(rsid)
        if other is None:
            counts["not_in_independent_2015"]+=1
            continue
        counts[orientation(other[0],other[1],a1,a2)]+=1
    both=sum(counts[k] for k in GOOD)
    if strict and (sum(counts.values())!=1_007_899 or both<MIN_SHARED_BMI):
        raise ValueError("unverified cross-cohort variant intersection")
    return {"total_Yengo_candidate":sum(counts.values()),
            "both_original_BMI_order_comparable":both,
            "categories":dict(counts),
            "2015_order_reversed_vs_2018":counts["swapped"]+counts["complement_swapped"],
            "individual_participant_independence_confirmed":False}

def compare_rows(records, locke, yengo_positions):
    counts=Counter()
    examples=[]
    visited=set()
    for fields in records:
        if len(fields)!=len(COLUMNS):
            raise ValueError("original GWAS column count changed")
        counts["source_complete_rows"]+=1
        ch,bp=fields[0],fields[1]
        if not ch.isdecimal() or not bp.isdecimal() or not 1<=int(ch)<=22 or int(bp)<=0:
            counts["nonautosomal_or_invalid"]+=1
            continue
        key=(int(ch),int(bp))
        bmi=yengo_positions.get(key)
        if bmi is None:
            counts["no_2018_candidate"]+=1
            continue
        rsid,ya1,ya2=bmi
        old=locke.get(rsid)
        if old is None:
            counts["no_2015_independent_candidate"]+=1
            continue
        if key in visited:
            raise ValueError("duplicate GWAS coordinate within original bounded sample")
        visited.add(key)
        ea,oa=fields[9].upper(),fields[10].upper()
        label2018=orientation(ya1,ya2,ea,oa)
        label2015=orientation(old[0],old[1],ea,oa)
        counts["against_2018_"+label2018]+=1
        counts["against_2015_"+label2015]+=1
        counts["triple_intersection"]+=1
        if len(examples)<5:
            examples.append({"2018_class":label2018,"2015_independent_class":label2015})
    counts["triple_intersection"]+=0  # explicit zero when no independent source matches
    if counts["source_complete_rows"]!=sum(counts[k] for k in (
       "nonautosomal_or_invalid","no_2018_candidate",
       "no_2015_independent_candidate","triple_intersection")):
        raise ValueError("original bounded source partition mismatch")
    return {"counts":dict(counts),"category_examples_without_effects":examples,
            "GIGASTROKE_vs_Locke_sign_reversals":
              counts["against_2015_swapped"]+counts["against_2015_complement_swapped"],
            "GIGASTROKE_vs_Yengo_sign_reversals":
              counts["against_2018_swapped"]+counts["against_2018_complement_swapped"]}

def verify_prefix_rows(phenotype,accession,compressed):
    qc=validate_prefix(phenotype,accession,compressed)
    d=zlib.decompressobj(16+zlib.MAX_WBITS)
    data=d.decompress(compressed,1_000_000)
    if len(data)>=1_000_000:
        raise ValueError("original publisher prefix inflated beyond cap")
    lines=data.split(b"\n")[:-1]
    if tuple(lines.pop(0).decode("utf-8-sig").strip("\r").split("\t"))!=COLUMNS:
        raise ValueError("original outcome header changed")
    records=[line.rstrip(b"\r").decode("ascii","strict").split("\t") for line in lines]
    if any(len(line)>4096 for line in lines):
        raise ValueError("GWAS prefix row oversized")
    if (len(records)!=EXPECTED_ROWS[phenotype] or
        qc["observed"]["parsed_effect_rows"]!=EXPECTED_ROWS[phenotype]):
        raise ValueError("pinned original GWAS prefix count drift")
    return records

def fetch_official_range(row, opener=urllib.request.urlopen):
    raw,_,_=source_urls(row)
    req=urllib.request.Request(raw,headers={
        "User-Agent":"CKM-Stroke-Independent-BMI-GWAS-Allele-Audit",
        "Range":f"bytes=0-{READ_BYTES-1}"})
    with opener(req,timeout=75) as res:
        url=urlsplit(res.geturl())
        if (res.status!=206 or url.scheme!="https" or url.hostname!="ftp.ebi.ac.uk"
            or not res.headers.get("Content-Range","").startswith(f"bytes 0-{READ_BYTES-1}/")):
            raise ValueError("original EUR GWAS bounded HTTP range untrusted")
        data=res.read(READ_BYTES+1)
    if len(data)!=READ_BYTES:
        raise ValueError("original GWAS incomplete bounded prefix")
    return data

def audit(locke_path,locke_provenance,yengo_path,yengo_provenance,manifest,*,opener=urllib.request.urlopen):
    result={"schema":"ckm-stroke-2015-independent-gwas-allele-order-prefix-v1",
            "scientific_status":"HOLD","ready_for_ldsc":False,
            "approved_outcome_N_eff":False,"LDSC_h2_rg_MR_coloc_run":False,
            "original_four_entire_GWAS_redownloaded":False}
    try:
        old,old_metadata=original_locke_2015_panel(locke_path,locke_provenance)
        current,current_metadata=build_bmi_coordinate_panel(yengo_path,yengo_provenance)
        bmi=compare_source_orders(old,current)
        entries=load_manifest(manifest)
        if len(entries)!=4 or {e["phenotype"] for e in entries}!=set(ACCESSIONS):
            raise ValueError("missing four official EUR GWAS original accessions")
        samples=[]
        for entry in entries:
            ph,acc=entry["phenotype"],entry["accession"]
            if ACCESSIONS[ph]!=acc:
                raise ValueError("CE/LAA or other original accession mismatch")
            rows=verify_prefix_rows(ph,acc,fetch_official_range(entry,opener))
            summary=compare_rows(rows,old,current)
            if summary["counts"]["triple_intersection"]<MIN_PREFIX_MATCH:
                raise ValueError("too few independent three-study original SNP matches")
            samples.append({"phenotype":ph,"accession":acc,
                            "source_complete_prefix_rows":len(rows),
                            "bounded_original_source_MD5_verified_here":False,
                            **summary})
        result.update({
            "technical_status":"PASS_INDEPENDENT_LOCKE_ALLELE_ORDER_BOUNDED_GWAS_ONLY",
            "checked_at_utc":datetime.now(timezone.utc).isoformat(),
            "Locke2015_independent":old_metadata,
            "Yengo2018_pinned":current_metadata,
            "independent_BMI_order_comparison":bmi,
            "original_GIGASTROKE_bounded_prefixes":samples,
            "caveat":"A1/A2 source-order comparisons do NOT establish biological direction across all original GWAS SNPs or participant independence."
        })
    except Exception as exc:
        result["technical_status"]="HOLD_INDEPENDENT_SOURCE_DIRECTION_VALIDATION"
        result["error_type"]=type(exc).__name__
        result["error_detail"]=str(exc)[:260]
    return result

def main():
    ap=argparse.ArgumentParser()
    for name in ("locke","locke-provenance","yengo","yengo-provenance"):
        ap.add_argument("--"+name,required=True,type=Path)
    ap.add_argument("--manifest",type=Path,default=Path("config/ckm_stroke/gigastroke_eur_official.tsv"))
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    report=audit(a.locke,a.locke_provenance,a.yengo,a.yengo_provenance,a.manifest)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":report["technical_status"],
        "2015_2018_shared":report.get("independent_BMI_order_comparison",{}).get("both_original_BMI_order_comparable"),
        "GIGASTROKE_sampled":{x["phenotype"]:x["counts"]["triple_intersection"] for x in report.get("original_GIGASTROKE_bounded_prefixes",[])},
        "2015_sign_reversal":{x["phenotype"]:x["GIGASTROKE_vs_Locke_sign_reversals"] for x in report.get("original_GIGASTROKE_bounded_prefixes",[])},
        "scientific":report["scientific_status"],"error":report.get("error_detail")
    },ensure_ascii=False))
    if report["technical_status"]!="PASS_INDEPENDENT_LOCKE_ALLELE_ORDER_BOUNDED_GWAS_ONLY":
        raise SystemExit(1)

if __name__=="__main__":
    main()
