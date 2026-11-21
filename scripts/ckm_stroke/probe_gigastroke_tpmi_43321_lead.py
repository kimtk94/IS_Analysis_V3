#!/usr/bin/env python3
"""Original MD5-verified GIGASTROKE EUR each outcome at one TPMI lead.

Full original publisher gzip streamed once per outcome, 4 parallel matrix
shards with max-parallel=2, never re-host or retain full summary data. Exact
GRCh37 coordinate and allele direction from independently Ensembl-mapped
TPMI rs187759247. Missing outcome association is *NOT* a replication.
"""
from __future__ import annotations
import argparse
from collections import Counter
import gzip
import io
import hashlib
import json
import math
from pathlib import Path
import urllib.request
from urllib.parse import urlsplit

from scripts.ckm_stroke.verify_gigastroke_raw import (
    source_urls, load_manifest, check_sidecars, get_text, MAX_BYTES
)
from scripts.ckm_stroke.review_gigastroke_full_evidence import (
    PINNED, OBSERVED_COUNTS
)
from scripts.ckm_stroke.audit_gigastroke_bmi_allele_harmonization import (
    COLUMNS, HashCompressedReader, orientation, GOOD, aligned_effect
)
from scripts.ckm_stroke.map_tpmi_43321_lead_ensembl import SHA as TPMI_SHA

TARGET_RSID="rs187759247"
TARGET_GRCH37=("19",15573379)
TARGET_REF="T"
TARGET_EFFECT="C"
TARGET_TPMI_BETA=0.5215696989267228
TARGET_TPMI_SE=0.0662893
EXPECTED_TRAITS=("AIS","LAA","CE","SVD")

def gate_mapped(mapping):
    src=mapping.get("lead_GRCh38",{})
    to37=mapping.get("mapped_GRCh37",{})
    if (mapping.get("schema")!="ckm-TPMI43321-lead-variant-Ensembl-GRCh38-GRCh37-v1"
        or mapping.get("technical_status")!="PASS_BIDIRECTIONAL_ENSEMBL_BUILD_MAP_AND_REF_ALLELE"
        or mapping.get("original_source_SHA256_verified")!=TPMI_SHA
        or mapping.get("liftover_reverse_verified") is not True
        or src.get("reported_rsid")!=TARGET_RSID
        or src.get("chrom_GRCh38")!="19"
        or src.get("position_GRCh38")!=15462568
        or to37.get("chrom_GRCh37")!=TARGET_GRCH37[0]
        or to37.get("pos_GRCh37")!=TARGET_GRCH37[1]
        or to37.get("REF_GRCh37")!=TARGET_REF
        or to37.get("A1_effect_allele_GRCh37")!=TARGET_EFFECT
        or mapping.get("scientific_status")!="HOLD"
        or mapping.get("ready_for_ldsc") is not False):
        raise ValueError("original TPMI lead/Ensembl source variant provenance drift")
    return True

def effect_source(fields):
    if len(fields)!=len(COLUMNS):
        raise ValueError("source original GWAS row width drift at target")
    chrom,bp,eaf,beta,se,pvalue,odds,cil,cih,ea,oa=fields
    if (chrom,bp)!=(TARGET_GRCH37[0],str(TARGET_GRCH37[1])):
        raise ValueError("wrong original study target genomic locus")
    orient=orientation(TARGET_EFFECT,TARGET_REF,ea.upper(),oa.upper())
    row={
        "trait_target_GRCh37_chr":chrom,
        "trait_target_GRCh37_position":int(bp),
        "reported_GIGASTROKE_effect_allele":ea,
        "reported_GIGASTROKE_other_allele":oa,
        "allele_orientation_relative_to_TPMI_A1_C":orient,
        "original_effect_frequency":float(eaf) if eaf not in ("NA","",".") else None,
    }
    if orient not in GOOD:
        row["technical_harmonization_status"]="HOLD_ALLELE_ORIENTATION"
        row["allele_aligned_effects_verified"]=False
        return row
    try:
        beta=float(beta);se=float(se);p=float(pvalue)
        if not all(map(math.isfinite,(beta,se,p))) or se<=0 or not 0<=p<=1:
            raise ValueError("source effect, standard error or p not finite")
        adjusted,z=aligned_effect(beta,se,orient)
    except ValueError:
        row["technical_harmonization_status"]="HOLD_EFFECT_NUMERIC_INVALID"
        row["allele_aligned_effects_verified"]=False
        return row
    row.update({
        "technical_harmonization_status":"PASS_ORIGINAL_ALLELE_DIRECTION_AND_NUMERIC",
        "allele_aligned_effects_verified":True,
        "original_GIGASTROKE_beta_logOR_aligned_C":adjusted,
        "original_GIGASTROKE_SE":se,
        "original_GIGASTROKE_P":p,
        "original_GIGASTROKE_aligned_Z":z,
        "raw_GIGASTROKE_effect_allele_frequency":row["original_effect_frequency"],
        "TPMI_original_beta_logOR_aligned_C":TARGET_TPMI_BETA,
        "TPMI_original_SE":TARGET_TPMI_SE,
        "effect_sign_concordant":(adjusted>=0)==(TARGET_TPMI_BETA>=0),
        "GIGASTROKE_original_samples_case_control_per_variant_known":False,
        "TPMI_vs_GIGASTROKE_pooled_or_meta_analysis_NOT_PERFORMED":True,
    })
    return row

def stream_full_one(trait,mapping,fetcher=urllib.request.urlopen):
    gate_mapped(mapping)
    if trait not in EXPECTED_TRAITS:
        raise ValueError("invalid original GIGASTROKE EUR trait")
    rows=load_manifest(Path("config/ckm_stroke/gigastroke_eur_official.tsv"))
    manifest=next(x for x in rows if x["phenotype"]==trait)
    acc,published_md5,size=PINNED[trait]
    if (manifest["accession"]!=acc or manifest["published_raw_md5"]!=published_md5
        or manifest.get("build")!="GRCh37"):
        raise ValueError("GIGASTROKE publisher source manifest drift")
    url,meta,sums=source_urls(manifest)
    check_sidecars(manifest,get_text(meta,fetcher),get_text(sums,fetcher))
    region=(TARGET_GRCH37[0],str(TARGET_GRCH37[1]))
    matches=[]
    n=0
    with fetcher(urllib.request.Request(url,headers={
            "User-Agent":"CKM-IS-Targeted-Publisher-MD5-rs187759247/1.0"}),
            timeout=60) as resp:
        actual=urlsplit(resp.geturl())
        if (resp.status!=200 or actual.scheme!="https"
            or actual.hostname!="ftp.ebi.ac.uk"):
            raise ValueError("GIGASTROKE original publisher source redirect or HTTP drift")
        header_content_length=resp.headers.get("Content-Length")
        if header_content_length and int(header_content_length)!=size:
            raise ValueError("published compressed source length mismatch")
        reader=HashCompressedReader(resp,cap=MAX_BYTES)
        with gzip.GzipFile(fileobj=reader,mode="rb") as compressed:
            with io.TextIOWrapper(compressed,encoding="utf-8",newline="") as textfd:
                header=textfd.readline().rstrip("\r\n").split("\t")
                if tuple(header)!=COLUMNS:
                    raise ValueError("original full GWAS header changed")
                for line in textfd:
                    if len(line)>8192:
                        raise ValueError("unbounded original GIGASTROKE GWAS row")
                    fields=line.rstrip("\r\n").split("\t")
                    if len(fields)!=len(COLUMNS):
                        raise ValueError("original GIGASTROKE association row malformed")
                    n+=1
                    if (fields[0],fields[1])==region:
                        matches.append(effect_source(fields))
                        if len(matches)>3:
                            raise ValueError("multiple GIGASTROKE target-locus duplicates; refuse selection")
                    if n>OBSERVED_COUNTS[trait][0]:
                        raise ValueError("original GWAS row count exceeds pinned whole-genome count")
        # Drain all original compressed bytes to account for any trailing data.
        while reader.read(1_048_576):pass
    if (n!=OBSERVED_COUNTS[trait][0]
        or reader.read_bytes!=size
        or reader.md5.hexdigest()!=published_md5):
        raise ValueError("original publisher compressed full MD5, bytes or row count invalid")
    if len(matches)>1:
        targetstatus="HOLD_DUPLICATE_OR_MULTIALLELIC_TARGET"
    elif len(matches)==1:
        targetstatus=matches[0]["technical_harmonization_status"]
    else:
        targetstatus="NOT_IN_ORIGINAL_SOURCE_GWAS"
    return {
        "schema":"ckm-TPMI43321-GIGASTROKE-four-subtype-original-variant-probe-v1",
        "technical_status":"PASS_ORIGINAL_PUBLISHER_MD5_AND_EXACT_TARGET_LOCUS_PROBE",
        "trait":trait,
        "original_source_accession":acc,
        "original_source_genome_build":"GRCh37",
        "original_publisher_MD5":published_md5,
        "original_MD5_verified":True,
        "original_compressed_bytes":reader.read_bytes,
        "original_all_rows_streamed":n,
        "original_expected_all_rows":OBSERVED_COUNTS[trait][0],
        "searched_TPMI_variant_rsID":TARGET_RSID,
        "searched_GRCh37_chrom":TARGET_GRCH37[0],
        "searched_GRCh37_position":TARGET_GRCH37[1],
        "exact_locus_original_row_count":len(matches),
        "one_locus_original_source_rows":matches,
        "target_effect_aligned_and_valid":len(matches)==1 and targetstatus=="PASS_ORIGINAL_ALLELE_DIRECTION_AND_NUMERIC",
        "target_source_status":targetstatus,
        "absence_does_NOT_demonstrate_no_effect":True,
        "multi_cohort_STUDY_sample_overlap_checked":False,
        "scientific_independent_replication_asserted":False,
        "SNP_case_control_effective_N_available":False,
        "LDSC_MR_coloc_not_executed":True,
        "scientific_status":"HOLD","ready_for_ldsc":False
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--trait",choices=EXPECTED_TRAITS,required=True)
    ap.add_argument("--mapped",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    try:
        source=json.loads(a.mapped.read_text("utf-8"))
        result=stream_full_one(a.trait,source)
    except Exception as exc:
        result={"technical_status":"HOLD_GIGASTROKE_ORIGINAL_MD5_OR_TARGET_SOURCE",
            "trait":a.trait,"error_type":type(exc).__name__,
            "error_detail":str(exc)[:300],"scientific_status":"HOLD",
            "ready_for_ldsc":False,"LDSC_MR_coloc_not_executed":True}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n","utf-8")
    print(json.dumps({"trait":a.trait,
        "status":result["technical_status"],
        "locus":result.get("target_source_status"),
        "matches":result.get("exact_locus_original_row_count"),
        "effect":result.get("one_locus_original_source_rows"),
        "error":result.get("error_detail"),
        "science":result["scientific_status"]}))
    if result["technical_status"]!="PASS_ORIGINAL_PUBLISHER_MD5_AND_EXACT_TARGET_LOCUS_PROBE":
        raise SystemExit(1)

if __name__=="__main__":
    main()
