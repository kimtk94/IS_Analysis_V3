#!/usr/bin/env python3
"""Resolve chr11 HRS frequency outlier identity and chr4 P/Wald deviations.

Original EUR stroke gzip streams are read once, publisher full-byte MD5 and SHA
are verified, only chromosome 4/11 preselected Yengo2018/HM3/LDscore candidates
are evaluated, and previous full 4-GWAS aggregate chr QC must be reproduced.
Public identifiers and diagnostics for at most a bounded number of outliers
are exported; no genomewide effect/beta tables or LDSC analyses are produced.
"""
from __future__ import annotations
import argparse
from collections import Counter
from decimal import Decimal, InvalidOperation
import gzip
import hashlib
import json
import math
from pathlib import Path
from urllib.parse import urlsplit
import urllib.request

from scripts.ckm_stroke.audit_gigastroke_bmi_allele_harmonization import (
    COLUMNS, GOOD, HashCompressedReader, MAX_GZIP_LINE, EXPECTED_GZIP_SIZES,
    build_bmi_coordinate_panel, orientation,
)
from scripts.ckm_stroke.audit_gigastroke_eaf_gc import (
    diagnostic, summarize, DIAGNOSTIC_KEYS, KEEP_FREQ,
)
from scripts.ckm_stroke.audit_bmi_duplicate_zero_p import minus_log10_two_sided_p_from_z
from scripts.ckm_stroke.audit_gigastroke_eaf_gc import review_diagnostics
from scripts.ckm_stroke.verify_gigastroke_raw import (
    load_manifest, source_urls, check_sidecars, get_text, MAX_BYTES,
)

TARGET_CHRS=(4,11)
MAX_TARGET_HITS=800_000
MAX_OUTPUT_FREQUENCY=8
MAX_OUTPUT_WALD=40
MAX_SOURCE_ROWS=35_000_000

def rounded_decimal_interval(text):
    """Half-unit of last printed source decimal digit, not measurement SE."""
    try:
        d=Decimal(text)
        if not d.is_finite():
            raise ValueError("nonfinite decimal input")
        step=Decimal(5).scaleb(d.as_tuple().exponent-1)
        return float(d-step),float(d+step)
    except (InvalidOperation,OverflowError,TypeError) as exc:
        raise ValueError("source decimal rounding bound not computable") from exc

def p_compatible_with_source_rounding(fields):
    """Whether printed P could arise from rounded BETA and SE (not proof)."""
    p_text=fields[5]
    if float(p_text)==0:
        return None
    b_lo,b_hi=rounded_decimal_interval(fields[3])
    s_lo,s_hi=rounded_decimal_interval(fields[4])
    p_lo,p_hi=rounded_decimal_interval(p_text)
    if s_lo<=0 or p_hi<=0:
        return None
    vals=[minus_log10_two_sided_p_from_z(b/s)
          for b in (b_lo,b_hi) for s in (s_lo,s_hi)]
    low,high=min(vals),max(vals)
    report_min=-math.log10(min(1,p_hi))
    report_max=float("inf") if p_lo<=0 else -math.log10(p_lo)
    return not (high<report_min-1e-12 or low>report_max+1e-12)

def original_chr_target_lines(rows,panel,hrs_by_coord):
    """Target chr4/11, deduplicate coordinates and count original source lines."""
    records={}
    collided=set()
    original_rows=0
    raw_header=next(rows,None)
    if raw_header is None or tuple(raw_header.rstrip("\r\n").split("\t"))!=COLUMNS:
        raise ValueError("original GIGASTROKE EUR GWAS source header drift")
    for line in rows:
        original_rows+=1
        if original_rows>MAX_SOURCE_ROWS or len(line)>MAX_GZIP_LINE:
            raise ValueError("original GWAS source row safety cap exceeded")
        fields=line.rstrip("\r\n").split("\t")
        if len(fields)!=len(COLUMNS):
            raise ValueError("original source row width drift")
        if fields[0] not in ("4","11"):
            continue
        if not fields[1].isdecimal():
            continue
        ch=int(fields[0]); pos=int(fields[1])
        coord=(ch,pos)
        if coord not in panel:
            continue
        if coord in collided:
            continue
        if coord in records:
            records.pop(coord)
            collided.add(coord)
            continue
        bmi=panel[coord]
        order=orientation(bmi[1],bmi[2],fields[9].upper(),fields[10].upper())
        if order not in GOOD:
            records[coord]=(None,None,None)
        else:
            diag=diagnostic(fields,hrs_by_coord.get(coord),order)
            # Hold only two chromosomes' candidate diagnostic data.
            records[coord]=(diag,order,(bmi[0],fields,hrs_by_coord.get(coord)))
        if len(records)+len(collided)>MAX_TARGET_HITS:
            raise ValueError("target chromosome source candidate cap exceeded")
    retained={coord:row for coord,row in records.items() if row[0] is not None}
    return original_rows,retained,len(collided),len(records)

def check_against_original_chr_evidence(phenotype,source,retained):
    by_chr={ch:{coord:row[0] for coord,row in retained.items() if coord[0]==ch}
            for ch in TARGET_CHRS}
    for ch in TARGET_CHRS:
        observed=summarize(by_chr[ch])["by_chr"][ch-1]
        expected=source["eaf_gc_source_descriptive"]["by_chr"][ch-1]
        old_chr=source["chr_attrition_pre_munge"][ch-1]
        if observed["matched_diagnostic_SNPs"]!=expected["matched_diagnostic_SNPs"] or (
            observed["matched_diagnostic_SNPs"]!=old_chr["allele_aligned_numeric_valid_pre_munge"]
        ):
            raise ValueError("target chr original SNP count does not match full-source "+phenotype)
        for key in DIAGNOSTIC_KEYS:
            if observed["categories"].get(key,0)!=expected["categories"].get(key,0):
                raise ValueError("target chr full source diagnostic category drift: "+phenotype+" chr"+str(ch))

def targeted_records(retained,source):
    freqs=[]
    wald=[]
    for (ch,pos),(diag,order,(rsid,fields,hrs)) in retained.items():
        if ch==11 and diag[1]=="hrs_absdiff_ge_020":
            if hrs is None:
                raise ValueError("extreme frequency without HRS reference")
            eaf=float(fields[2]);aligned=eaf if order in KEEP_FREQ else 1-eaf
            freqs.append({
                "chromosome":ch,"position_GRCh37":pos,"rsid_from_pinned_GIANT":rsid,
                "source_reference_alleles_sorted":"".join(sorted(fields[9:11])),
                "allele_oriented_to_GIANT_tested":True,
                "GIANT_HRS_tested_allele_frequency":round(hrs,9),
                "GIGASTROKE_frequency_oriented_to_GIANT":round(aligned,9),
                "absolute_source_HRS_difference":round(abs(aligned-hrs),9),
            })
        if ch==4 and diag[2]=="abs_log10p_delta_010_to_1":
            b=float(fields[3]);s=float(fields[4]);p=float(fields[5])
            if p<=0:
                raise ValueError("P/Wald delta anomaly cannot be underflowed zero")
            observed=-math.log10(p)
            implied=minus_log10_two_sided_p_from_z(b/s)
            wald.append({
                "chromosome":ch,"position_GRCh37":pos,"rsid_from_pinned_GIANT":rsid,
                "reported_P_decimal_digits":len(fields[5]),
                "BETA_decimal_digits":len(fields[3]),
                "SE_decimal_digits":len(fields[4]),
                "abs_log10P_difference":round(abs(observed-implied),9),
                "reported_P_compatible_with_rounded_BETA_SE":p_compatible_with_source_rounding(fields),
                "source_association_effects_exported":False,
            })
    freqs.sort(key=lambda x:(x["chromosome"],x["position_GRCh37"],x["rsid_from_pinned_GIANT"]))
    wald.sort(key=lambda x:(x["chromosome"],x["position_GRCh37"],x["rsid_from_pinned_GIANT"]))
    expected_freq=source["eaf_gc_source_descriptive"]["by_chr"][10]["categories"].get("hrs_absdiff_ge_020",0)
    expected_wald=source["eaf_gc_source_descriptive"]["by_chr"][3]["categories"].get("abs_log10p_delta_010_to_1",0)
    if len(freqs)!=expected_freq or len(wald)!=expected_wald:
        raise ValueError("rare chr11 or chr4 original outlier count drift")
    if len(freqs)>MAX_OUTPUT_FREQUENCY or len(wald)>MAX_OUTPUT_WALD:
        raise ValueError("too many targeted variants for bounded identity export")
    return freqs,wald

def scan_original(row,panel,hrs,source,*,fetcher=urllib.request.urlopen):
    raw,meta,sums=source_urls(row)
    check_sidecars(row,get_text(meta,fetcher),get_text(sums,fetcher))
    req=urllib.request.Request(raw,headers={"User-Agent":"CKM-Stroke-Chr4-11-MD5-Triage/1.0"})
    with fetcher(req,timeout=90) as response:
        dest=urlsplit(response.geturl())
        if response.status!=200 or dest.scheme!="https" or dest.hostname!="ftp.ebi.ac.uk":
            raise ValueError("publisher original whole-file status/host changed")
        size=EXPECTED_GZIP_SIZES[row["phenotype"]]
        if response.headers.get("Content-Length")!=str(size):
            raise ValueError("publisher original length changed")
        reader=HashCompressedReader(response,MAX_BYTES)
        with gzip.GzipFile(fileobj=reader,mode="rb") as gz:
            def lines():
                while True:
                    chunk=gz.readline(MAX_GZIP_LINE+1)
                    if not chunk:
                        break
                    if len(chunk)>MAX_GZIP_LINE:
                        raise ValueError("source gzip line exceeds bound")
                    yield chunk.decode("utf-8-sig","strict")
            count,retained,collisions,valid_positions=original_chr_target_lines(lines(),panel,hrs)
        while reader.read(1048576):
            pass
        if (reader.read_bytes!=size or reader.md5.hexdigest()!=row["published_raw_md5"].lower()
            or reader.sha256.hexdigest()!=source["source_sha256"]):
            raise ValueError("publisher MD5/SHA256 full-source mismatch")
        if count!=source["counts"]["original_outcome_rows"]:
            raise ValueError("total original source rows changed")
        check_against_original_chr_evidence(row["phenotype"],source,retained)
        f,w=targeted_records(retained,source)
        return {
            "phenotype":row["phenotype"],"accession":row["accession"],
            "full_source_publisher_md5_and_prior_sha_verified":True,
            "total_original_source_rows":count,
            "original_source_bytes":reader.read_bytes,
            "chr4_11_duplicate_coordinates_quarantined":collisions,
            "chr4_11_distinct_noncandidate_excluded":True,
            "source_frequency_outliers_chr11":f,
            "source_p_wald_outliers_chr4":w,
            "source_effects_stored":False,
        }

def audit(candidate,provenance,manifest,original_evidence,*,fetcher=urllib.request.urlopen):
    previous=review_diagnostics(original_evidence)
    panel,meta,hrs=build_bmi_coordinate_panel(candidate,provenance,collect_hrs=True)
    if (meta["candidate_plaintext_sha256"]!=original_evidence["BMI_reference_panel"]["candidate_plaintext_sha256"]
        or len(panel)!=1007899):
        raise ValueError("candidate reference panel differs from pinned whole-source evidence")
    records=load_manifest(manifest)
    sources={item["phenotype"]:item for item in original_evidence["outcomes"]}
    results=[]
    for row in records:
        ph=row["phenotype"]
        if ph not in sources or row["published_raw_md5"]!=sources[ph]["published_md5"]:
            raise ValueError("original accession/MD5 is no longer pinned")
        results.append(scan_original(row,panel,hrs,sources[ph],fetcher=fetcher))
    observations=[x for r in results for x in r["source_frequency_outliers_chr11"]]
    locuskeys={(x["position_GRCh37"],x["rsid_from_pinned_GIANT"],x["source_reference_alleles_sorted"])
               for x in observations}
    same=(len(observations)==4 and len(locuskeys)==1)
    walds=[x for r in results for x in r["source_p_wald_outliers_chr4"]]
    wald_loci_by_trait={x["phenotype"]:{(v["position_GRCh37"],v["rsid_from_pinned_GIANT"])
                                        for v in x["source_p_wald_outliers_chr4"]}
                        for x in results}
    shared_AIS_CE=len(wald_loci_by_trait["AIS"] & wald_loci_by_trait["CE"])
    return {
        "schema":"ckm-GIGASTROKE-EUR-chr11-variant-identity-chr4-P-Wald-v1",
        "technical_status":"PASS_PINNED_FULL_MD5_TARGET_CHR4_11_VARIANT_IDENTITY",
        "based_on_previous_original_EAF_Wald_review":previous["technical_status"],
        "full_source_reverified_four_original_MD5_and_SHA256":True,
        "all_four_chr11_outlier_same_genomic_variant":same,
        "chr11_outlier_variant_identity":observations[0] if same else None,
        "chr11_source_outlier_frequency_cross_trait_details":observations,
        "chr4_P_Wald_moderate_outlier_shared_AIS_CE_variant_count":shared_AIS_CE,
        "chr4_source_precision_diagnostics_by_trait":results,
        "source_eaf_gc_factor_estimated":False,
        "original_allele_effects_or_complete_SNP_tables_saved":False,
        "outcome_per_variant_N_eff_approved":False,
        "LDSC_h2_rg_MR_coloc_executed":False,
        "scientific_status":"HOLD","ready_for_ldsc":False,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--candidate",required=True,type=Path)
    p.add_argument("--provenance",required=True,type=Path)
    p.add_argument("--manifest",type=Path,default=Path("config/ckm_stroke/gigastroke_eur_official.tsv"))
    p.add_argument("--previous-original-evidence",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    try:
        original=json.loads(args.previous_original_evidence.read_text("utf-8"))
        report=audit(args.candidate,args.provenance,args.manifest,original)
    except Exception as exc:
        report={
            "technical_status":"HOLD_CHR11_IDENTITY_OR_CHR4_WALD_PROVENANCE",
            "error_type":type(exc).__name__,"error_detail":str(exc)[:300],
            "scientific_status":"HOLD","ready_for_ldsc":False,
            "outcome_per_variant_N_eff_approved":False,
            "LDSC_h2_rg_MR_coloc_executed":False,
        }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n","utf-8")
    print(json.dumps({"status":report["technical_status"],
                      "same_chr11_variant":report.get("all_four_chr11_outlier_same_genomic_variant"),
                      "chr4_AIS_CE_shared_variant_count":report.get("chr4_P_Wald_moderate_outlier_shared_AIS_CE_variant_count"),
                      "science":report["scientific_status"],"error":report.get("error_detail")}))
    if report["technical_status"]!="PASS_PINNED_FULL_MD5_TARGET_CHR4_11_VARIANT_IDENTITY":
        raise SystemExit(1)

if __name__=="__main__":
    main()
