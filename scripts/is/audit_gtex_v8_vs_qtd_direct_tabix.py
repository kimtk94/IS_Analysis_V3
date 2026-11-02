#!/usr/bin/env python3
"""Bounded original GTEx V8 vs EBI reprocessed QTD000161 direct tabix comparison.

Two source archives, same 2020 GTEx v8 Brain Cerebellar Hemisphere tissue:
* imported original gene expression FastQTL results
* EBI eQTL Catalogue reprocessed gene expression QTD000161

Only narrow predeclared GRCh38 regions are queried; source URLs and original
MD5 manifest are pinned. No full 3.4GB downloads, paid GCS, SSH or coloc.
Never merge old GWAS or claim QTL model identity from SNP overlap alone.
"""
import argparse
import hashlib
import json
import math
import re
import statistics
import time
import urllib.request
from pathlib import Path

ROOT="https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/"
ORIGINAL=ROOT+"imported/GTEx_V8/ge/Brain_Cerebellar_Hemisphere.tsv.gz"
ORIGINAL_MD5_LIST=ROOT+"imported/GTEx_V8/ge/md5sum.txt"
ORIGINAL_MD5="c07a1731eb96aafbd7aed131cb97515a"
REPROCESSED=ROOT+"sumstats/QTS000015/QTD000161/QTD000161.all.tsv.gz"
REGIONS=[("4",79759900,79764900),("10",103079000,103084000)]
MAX_ROWS_EACH_REGION=25000
MAX_META_BYTES=20000
REQUIRED={"molecular_trait_id","chromosome","position","ref","alt","pvalue","beta","se"}

def verify_manifest(text):
    lines=text.splitlines()
    hits=[]
    for line in lines:
        match=re.fullmatch(r"([a-fA-F0-9]{32})\s+\*?([A-Za-z0-9._-]+)",line.strip())
        if match and match.group(2)=="Brain_Cerebellar_Hemisphere.tsv.gz":
            hits.append(match.group(1).lower())
    return len(hits)==1 and hits[0]==ORIGINAL_MD5

def read_md5_list(fetcher=None):
    if fetcher is not None:
        raw=fetcher(ORIGINAL_MD5_LIST)
    else:
        request=urllib.request.Request(ORIGINAL_MD5_LIST,
            headers={"User-Agent":"IS-Bounded-GTExv8-QTD-DirectComparison/1.0"})
        with urllib.request.urlopen(request,timeout=25) as response:
            if not response.url.startswith(ROOT+"imported/GTEx_V8/ge/"):
                raise ValueError("Unexpected MD5 manifest redirect")
            raw=response.read(MAX_META_BYTES+1)
    if len(raw)>MAX_META_BYTES:
        raise ValueError("MD5 list exceeds cap")
    return {"source_md5_manifest_sha256":hashlib.sha256(raw).hexdigest(),
            "original_tissue_file_md5_in_official_manifest":
                verify_manifest(raw.decode("utf-8","replace"))}

def normalize_header(lines):
    raw=[x.lstrip("#").strip() for x in lines]
    for line in raw:
        bits=line.split("\t")
        if {"beta","pvalue","se"}.issubset(set(bits)):
            return bits
    return None

def known_index_contig(contigs,chrom):
    if chrom in contigs:return chrom
    if "chr"+chrom in contigs:return "chr"+chrom
    return None

def select_qtl_rows(header, lines, *, line_limit=MAX_ROWS_EACH_REGION):
    index={name:i for i,name in enumerate(header)}
    if not REQUIRED.issubset(index):
        return {"status":"HOLD_MISSING_QTL_FIELDS",
                "columns":header,"rows":{},"observed":0}
    rows={}
    count=0
    for line in lines:
        count+=1
        if count>line_limit:
            return {"status":"HOLD_REGION_ROW_CAP","columns":header,
                    "rows":{},"observed":count}
        vals=line.rstrip("\n").split("\t")
        if len(vals)!=len(header):
            return {"status":"HOLD_MALFORMED_ROW","columns":header,
                    "rows":{},"observed":count}
        try:
            trait=vals[index["molecular_trait_id"]].split(".")[0]
            chrom=vals[index["chromosome"]].removeprefix("chr")
            pos=int(vals[index["position"]])
            ref=vals[index["ref"]].upper()
            alt=vals[index["alt"]].upper()
            p=float(vals[index["pvalue"]])
            beta=float(vals[index["beta"]])
            se=float(vals[index["se"]])
            if len(ref)!=1 or len(alt)!=1 or ref==alt or not all(x in "ACGT" for x in (ref,alt)):
                continue
            if not (math.isfinite(p) and 0<=p<=1 and math.isfinite(beta) and
                    math.isfinite(se) and se>0):
                continue
            key=(trait,chrom,pos,ref,alt)
            if key not in rows:
                rows[key]=(beta,se,p)
            elif rows[key]!=(beta,se,p):
                return {"status":"HOLD_CONFLICTING_DUPLICATE_KEY","columns":header,
                        "rows":{},"observed":count}
        except (TypeError,ValueError):
            continue
    return {"status":"SCREEN_VALID","columns":header,
            "rows":rows,"observed":count}

def compare_rows(original,reprocessed):
    first=original["rows"]
    second=reprocessed["rows"]
    common=set(first)&set(second)
    sign=0
    p_diff=[]
    beta_diff=[]
    sep=[]
    same_beta=0
    same_p=0
    for key in common:
        b1,s1,p1=first[key]
        b2,s2,p2=second[key]
        sign+=int((b1>=0)==(b2>=0))
        beta_diff.append(abs(b1-b2))
        sep.append(abs(s1-s2))
        if p1>0 and p2>0:
            p_diff.append(abs(math.log10(p1)-math.log10(p2)))
        same_beta+=int(math.isclose(b1,b2,abs_tol=1e-6,rel_tol=1e-5))
        same_p+=int(math.isclose(p1,p2,abs_tol=1e-9,rel_tol=1e-5))
    n=len(common)
    return {
        "original_region_rows":original["observed"],
        "reprocessed_region_rows":reprocessed["observed"],
        "original_distinct_valid_trait_variant_pairs":len(first),
        "reprocessed_distinct_valid_trait_variant_pairs":len(second),
        "matched_same_trait_and_exact_GRCh38_REF_ALT_pairs":n,
        "beta_sign_matches":sign,"beta_sign_fraction":sign/n if n else None,
        "beta_numerically_near_equal":same_beta,
        "p_numerically_near_equal":same_p,
        "median_abs_beta_difference":statistics.median(beta_diff) if beta_diff else None,
        "median_abs_SE_difference":statistics.median(sep) if sep else None,
        "median_abs_log10_p_difference":statistics.median(p_diff) if p_diff else None,
        "source_fields_complete":original["status"]=="SCREEN_VALID" and reprocessed["status"]=="SCREEN_VALID",
        "comparison_status":(
            "SOURCE_MODELS_DIFFERENT_OR_EQUALITY_REVIEW" if n>0
            else "HOLD_NO_EXACT_GENE_VARIANT_OVERLAP"),
        "sampled_region_limited":True,
        "causal_gene_inferred":False,
        "coloc_performed":False,
    }

def extract_region(tabix,url,region,*,sleep_seconds=2.25):
    chrom,start,end=region
    indexed=known_index_contig(set(tabix.contigs),chrom)
    if indexed is None:
        return {"status":"HOLD_CHROMOSOME_ABSENT","columns":[],
                "rows":{},"observed":0}
    # One bounded indexed fetch per original and reprocessed region.
    data=select_qtl_rows(normalize_header(tabix.header) or [], 
                         tabix.fetch(indexed,start-1,end))
    return data

def live_audit(tabix_module=None,md5_fetcher=None,wait_seconds=2.25):
    report={"schema":"IS_GTExv8_OriginalVsQTD000161_BoundedTabix_v1",
        "tissue":"Brain_Cerebellar_Hemisphere","GTEx_original_ftp":ORIGINAL,
        "reprocessed_QTL_ftp":REPROCESSED,"original_md5_manifest_expected":ORIGINAL_MD5,
        "reference_genome_context":"GTEx v8 GRCh38, must confirm reported REF/ALT in both",
        "scientific_status":"HOLD",
        "GTEx_original_complete_cis_universe_authenticated":False,
        "historical_inherited_source_proven_equal_to_GTEx_v8_original":False,
        "QTD000161_2023_release6_file_sha256_authenticated":False,
        "original_full_file_checksum_recomputed":False,
        "coloc_performed":False}
    try:
        report["official_original_md5_list"]=read_md5_list(md5_fetcher)
    except (OSError,ValueError,UnicodeError) as exc:
        report["official_original_md5_list"]={
            "original_tissue_file_md5_in_official_manifest":False,
            "failure_class":type(exc).__name__}
    if tabix_module is None:
        try:
            import pysam as tabix_module
        except ImportError:
            report["source_region_gate"]="HOLD_PYSAM_UNAVAILABLE"
            return report
    handles=[]
    try:
        for url in (ORIGINAL,REPROCESSED):
            try:
                h=tabix_module.TabixFile(url)
            except (OSError,ValueError) as exc:
                report["source_region_gate"]="HOLD_REMOTE_TABIX_SOURCE_UNAVAILABLE"
                report["failure_class"]=type(exc).__name__
                return report
            handles.append(h)
        report["source_headers"]={
            "original":normalize_header(handles[0].header),
            "reprocessed":normalize_header(handles[1].header)}
        report["source_contigs"]={
            "original_count":len(handles[0].contigs),
            "reprocessed_count":len(handles[1].contigs),
            "chr4_and_chr10_in_original":all(
                known_index_contig(handles[0].contigs,c) for c in ("4","10")),
            "chr4_and_chr10_in_reprocessed":all(
                known_index_contig(handles[1].contigs,c) for c in ("4","10")),
        }
        stats=[]
        for region in REGIONS:
            entries=[]
            for handle,url in zip(handles,(ORIGINAL,REPROCESSED)):
                time.sleep(wait_seconds)
                entries.append(extract_region(handle,url,region))
            record=compare_rows(*entries)
            record["region_GRCh38"]="%s:%d-%d"%region
            record["original_extract_status"]=entries[0]["status"]
            record["reprocessed_extract_status"]=entries[1]["status"]
            stats.append(record)
        report["sampled_regions"]=stats
        sum_matched=sum(x["matched_same_trait_and_exact_GRCh38_REF_ALT_pairs"] for x in stats)
        report["total_exact_pairs_two_small_regions"]=sum_matched
        report["source_region_gate"]=(
            "SMALL_REGION_GTEXv8_VS_QTD_QTL_NUMERIC_COMPARE_PASS_REVIEW_SOURCE_DIFFERENCES"
            if sum_matched>0 and all(x["source_fields_complete"] for x in stats)
            else "HOLD_INCOMPLETE_OR_NO_EXACT_QTL_OVERLAP")
    except (OSError,ValueError,KeyError,UnicodeError) as exc:
        report["source_region_gate"]="HOLD_REMOTE_REGION_FETCH_OR_SCHEMA"
        report["region_failure_class"]=type(exc).__name__
    finally:
        for h in handles:
            h.close()
    return report

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    report=live_audit()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("IS_GTEX_V8_QTD_TABIX="+json.dumps(report,sort_keys=True))

if __name__=="__main__":
    main()
