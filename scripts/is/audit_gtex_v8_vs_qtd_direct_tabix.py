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
import subprocess
import os
import time
import urllib.request
import zlib
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


def original_bgzf_first_line(url, fetcher=None):
    """Only read a bounded official BGZF prefix; Tabix headers may be unmarked."""
    if url not in (ORIGINAL,REPROCESSED):
        raise ValueError("QTL URL not allowlisted")
    if fetcher is None:
        req=urllib.request.Request(url, headers={
            "User-Agent":"IS-GTExv8-QTD-Header-Read/1.0",
            "Range":"bytes=0-131071"})
        with urllib.request.urlopen(req,timeout=30) as res:
            if res.url!=url:
                raise ValueError("Unexpected cross-origin or version redirect")
            prefix=res.read(131073)
    else:
        prefix=fetcher(url)
    if len(prefix)>131072 or not prefix.startswith(b"\x1f\x8b"):
        raise ValueError("Unexpected source BGZF prefix")
    first=zlib.decompressobj(wbits=16+zlib.MAX_WBITS).decompress(prefix,262144)
    if b"\n" not in first:
        raise ValueError("Incomplete BGZF header first line")
    first_line=first.split(b"\n",1)[0].decode("utf-8-sig","replace").strip()
    header=first_line.lstrip("#").split("\t")
    return {"header":header,
            "gzip_prefix_sha256_not_full_file":hashlib.sha256(prefix).hexdigest(),
            "header_from_first_bgzf_block":True}

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


class SystemTabixFile:
    """Use distribution htslib/libcurl with the system trust store, not bundled wheel."""
    def __init__(self,url):
        if url not in (ORIGINAL,REPROCESSED):
            raise ValueError("Tabix source not allowlisted")
        self.url=url
        # The EBI eQTL Catalogue data-access documentation explicitly supports
        # FTP Tabix. The official HTTPS root is retained for provenance.
        # Only switch transport, never hostname or exact published path.
        failures=[]
        for transport in (url,url.replace("https://","ftp://",1)):
            try:
                header=self._call(["tabix","-H",transport],cap=200_000).splitlines()
                contigs=self._call(["tabix","-l",transport],cap=20_000).splitlines()
                self.transport=transport
                self.header=header
                self.contigs=contigs
                break
            except (OSError,ValueError) as exc:
                failures.append(str(exc)[-160:])
        else:
            raise OSError("HTTPS and documented EBI FTP Tabix failed: "+" | ".join(failures))
    def _call(self,args,cap):
        env=os.environ.copy()
        env["CURL_CA_BUNDLE"]="/etc/ssl/certs/ca-certificates.crt"
        try:
            completed=subprocess.run(args,capture_output=True,check=False,
                                      timeout=110,env=env)
        except (OSError,subprocess.TimeoutExpired) as exc:
            raise OSError("system tabix invocation failed: "+type(exc).__name__) from exc
        if completed.returncode:
            tail=completed.stderr.decode("utf-8","replace")[-220:]
            raise OSError("system tabix command failed: "+str(completed.returncode)+": "+tail)
        if len(completed.stdout)>cap:
            raise ValueError("tabix output exceeded strict capture cap")
        return completed.stdout.decode("utf-8","replace")
    def fetch(self,chrom,start,end):
        if end-start>6000:
            raise ValueError("interval too broad")
        region=f"{chrom}:{start+1}-{end}"
        return iter(self._call(["tabix",self.transport,region],cap=8_000_000).splitlines())
    def close(self):
        return None

class SystemTabix:
    @staticmethod
    def TabixFile(url):
        return SystemTabixFile(url)

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
        # pysam wheel htslib has libcurl, but some Linux distributions need
        # an explicit CA bundle to avoid libcurl error 77 on HTTPS sources.
        os.environ["CURL_CA_BUNDLE"]="/etc/ssl/certs/ca-certificates.crt"
        os.environ["SSL_CERT_FILE"]="/etc/ssl/certs/ca-certificates.crt"
        try:
            import pysam as tabix_module
            report["source_tabix_backend"]="pysam_libcurl_explicit_system_CA"
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
                report["failure_reason_bounded"]=str(exc)[-250:]
                report["failure_url_type"]="original" if url==ORIGINAL else "reprocessed"
                return report
            handles.append(h)
        report["original_source_tabix_transport"]=handles[0].transport.split("://")[0] if hasattr(handles[0],"transport") else "mock"
        report["reprocessed_tabix_transport"]=handles[1].transport.split("://")[0] if hasattr(handles[1],"transport") else "mock"
        report["remote_transport_cryptographic_integrity_verified"]=False
        report["source_headers"]={}
        report["source_header_provenance"]={}
        for kind,handle,url in zip(("original","reprocessed"),handles,(ORIGINAL,REPROCESSED)):
            head=normalize_header(handle.header)
            if head is None:
                try:
                    source=original_bgzf_first_line(url)
                    head=source["header"]
                    report["source_header_provenance"][kind]={
                        "origin":"actual_original_gzip_Bgzf_prefix",
                        "prefix_sha256_not_full_file":source["gzip_prefix_sha256_not_full_file"],
                    }
                except (OSError,ValueError,zlib.error) as exc:
                    report["source_header_provenance"][kind]={
                        "origin":"HOLD_BGZF_HEADER_UNAVAILABLE",
                        "failure_type":type(exc).__name__}
            else:
                report["source_header_provenance"][kind]={"origin":"tabix_comment_header"}
            report["source_headers"][kind]=head
            handle.header=(["#"+"\t".join(head)] if head is not None else [])
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
