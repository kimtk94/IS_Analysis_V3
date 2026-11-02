#!/usr/bin/env python3
"""Bounded official GWAS Catalog FTP source-file audit for GCST90104545.

Uses public HTTPS only; no SSH, secrets, personal data, or full data artifacts.
In 'full' mode, compute compressed-file SHA256/MD5 and scan gzip records.
Source-file checksum/column presence is not variant harmonisation or coloc approval.
"""
import argparse
import csv
import gzip
import hashlib
from html.parser import HTMLParser
import json
import math
import re
from collections import Counter
import tempfile
import urllib.error
import urllib.parse
import urllib.request
import zlib
from pathlib import Path

ACCESSION = "GCST90104545"
BASE = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/" + ACCESSION + "/"
MAX_HTML = 2_000_000
MAX_YAML = 2_000_000
MAX_PREFIX = 262144
MAX_COMPRESSED = 1_500_000_000
MAX_DATA_ROWS = 20_000_000
REQUIRED_COLUMNS = {"chromosome", "base_pair_location", "effect_allele",
                    "other_allele", "p_value"}
SE_ALIASES = {"standard_error", "se", "beta_se", "sebeta"}
EFFECT_ALIASES = {"beta", "effect_size", "odds_ratio", "hazard_ratio"}
UA = "IS-Source-ReadOnly-Audit/20261102"
# Frozen from user's 2026-10-11 separate four-gene Drive crosswalk report.
# These are PUBLIC association-level p-value anchors, not original QTL lineage
# or a canonical-wide SNP-wise equivalence certification.
CANONICAL_PANEL = {
    ("4", 81182554): {"gene": "FGF5", "alleles": frozenset(("T","C")), "giga_ais_p": 2.888e-8},
    ("10", 105447782): {"gene": "CALHM2", "alleles": frozenset(("A","G")), "giga_ais_p": 1.332e-6},
    ("10", 105459834): {"gene": "CALHM2", "alleles": frozenset(("T","C")), "giga_ais_p": 7.792e-7},
}

class Hrefs(HTMLParser):
    def __init__(self):
        super().__init__()
        self.items = []
    def handle_starttag(self, tag, attrs):
        if tag == "a":
            url = dict(attrs).get("href")
            if url:
                self.items.append(url)

def safe_asset_url(filename):
    if not isinstance(filename, str):
        raise ValueError("Bad filename")
    if not re.fullmatch(r"GCST90104545[A-Za-z0-9._-]*", filename):
        raise ValueError("Asset outside exact accession namespace")
    if filename in (".","..") or ".." in filename:
        raise ValueError("Unsafe path")
    return BASE + urllib.parse.quote(filename, safe="._-")

def open_checked(url, *, method="GET", headers=None, opener=None):
    if not url.startswith(BASE):
        raise ValueError("Outside allowlisted accession")
    request = urllib.request.Request(url, headers={"User-Agent":UA, **(headers or {})},method=method)
    response = (opener or urllib.request.urlopen)(request, timeout=25)
    final = urllib.parse.urlsplit(response.url)
    base = urllib.parse.urlsplit(BASE)
    if final.scheme != "https" or final.hostname != base.hostname or not final.path.startswith(base.path):
        response.close()
        raise ValueError("Untrusted redirect")
    return response

def limited_get(url, limit, opener=None):
    with open_checked(url, opener=opener) as r:
        data = r.read(limit+1)
    if len(data) > limit:
        raise ValueError("Exceeded metadata byte budget")
    return data

def discover_assets(index_html):
    parser=Hrefs()
    parser.feed(index_html.decode("utf-8","replace"))
    files=[]
    for href in parser.items:
        base = urllib.parse.urlsplit(href).path.split("/")[-1]
        base = urllib.parse.unquote(base)
        if base.startswith(ACCESSION) and re.fullmatch(r"[A-Za-z0-9._-]+",base):
            if base not in files:
                files.append(base)
    return sorted(files)

def yaml_fields(raw):
    fields={}
    for line in raw.decode("utf-8","replace").splitlines():
        match=re.fullmatch(r"\s*([A-Za-z_][A-Za-z_0-9]*)\s*:\s*([^#]*?)\s*",line)
        if match:
            key,value=match.groups()
            fields[key]=value.strip().strip("'\"")
    return fields

def get_stats_headers(prefix):
    if prefix.startswith(b"\x1f\x8b"):
        decompressor=zlib.decompressobj(wbits=16+zlib.MAX_WBITS)
        text=decompressor.decompress(prefix, 1024*1024)
    else:
        return None
    lines=text.splitlines()
    if not lines or b"\t" not in lines[0]:
        return None
    return [x.strip().lower() for x in lines[0].decode("utf-8-sig","replace").split("\t")]

def probe(opener=None):
    result={
        "schema":"is-gigastroke-ais-original-ftp-gate-v1",
        "accession":ACCESSION,"directory_url":BASE,
        "metadata_checks":"HOLD",
        "original_file_byte_hash_verified":False,
        "original_gwas_N_SE_effect_allele_approved":False,
        "original_GRCh37_reference_mapping_approved":False,
        "source_to_local_EAS5_canonical_allele_alignment_verified":False,
        "coloc_performed":False,"scientific_status":"HOLD",
    }
    try:
        listing=limited_get(BASE,MAX_HTML,opener)
        result["directory_html_sha256"]=hashlib.sha256(listing).hexdigest()
        files=discover_assets(listing)
        result["discovered_asset_names"]=files
        metadata_files=[n for n in files if n.endswith("-meta.yaml")]
        if len(metadata_files)!=1:
            result["metadata_gate"]="HOLD_ZERO_OR_MULTIPLE_YAML_FILES"
            return result
        metadata_file=metadata_files[0]
        yaml_raw=limited_get(safe_asset_url(metadata_file),MAX_YAML,opener)
        y=yaml_fields(yaml_raw)
        result["yaml_filename"]=metadata_file
        result["yaml_sha256"]=hashlib.sha256(yaml_raw).hexdigest()
        exposed_keys=("gwas_id","data_file_name","data_file_md5sum","file_type",
                      "trait_description","case_count","control_count",
                      "genome_assembly","is_harmonised","harmonisation_reference")
        result["yaml_selected_fields"]={key:y[key] for key in exposed_keys if key in y}
        filename=y.get("data_file_name","")
        if filename not in files:
            result["metadata_gate"]="HOLD_YAML_DATA_FILE_NOT_IN_DIRECTORY"
            return result
        if not filename.endswith(".tsv.gz"):
            result["metadata_gate"]="HOLD_UNSUPPORTED_RAW_FILE_FORMAT"
            return result
        raw_md5=y.get("data_file_md5sum","")
        if not re.fullmatch(r"[0-9a-fA-F]{32}",raw_md5):
            result["metadata_gate"]="HOLD_MISSING_OFFICIAL_RAW_CHECKSUM"
            return result
        if y.get("gwas_id") != ACCESSION:
            result["metadata_gate"]="HOLD_MISMATCHED_ACCESSION"
            return result
        result["raw_filename"]=filename
        try:
            with open_checked(safe_asset_url(filename),method="HEAD",opener=opener) as head:
                length=head.headers.get("Content-Length")
                result["compressed_content_length"]=int(length) if length and length.isdecimal() else None
        except (OSError, ValueError):
            result["compressed_content_length"]=None
        with open_checked(safe_asset_url(filename),headers={"Range":f"bytes=0-{MAX_PREFIX-1}"},opener=opener) as response:
            prefix=response.read(MAX_PREFIX)
            result["prefix_response_status"]=getattr(response,"status",None)
        result["prefix_bytes_read"]=len(prefix)
        headers=get_stats_headers(prefix)
        result["original_gwas_header"]=headers
        if headers is None:
            result["metadata_gate"]="HOLD_UNREADABLE_HEADER"
        else:
            names=set(headers)
            result["column_presence"]={
                "mandatory_position_and_allele_and_p":REQUIRED_COLUMNS.issubset(names),
                "any_effect_estimate":bool(EFFECT_ALIASES & names),
                "standard_error":bool(SE_ALIASES & names),
            }
            result["metadata_gate"]=(
                "METADATA_HEADER_SCREEN_PASS" if all(result["column_presence"].values())
                else "HOLD_REQUIRED_OR_SE_COLUMN_NOT_FOUND")
        result["metadata_checks"]="COMPLETE_PROBE"
    except (OSError, ValueError, UnicodeError, zlib.error, EOFError) as e:
        result["metadata_gate"]="HOLD_NETWORK_OR_METADATA_PROBE"
        result["failure_type"]=type(e).__name__
    return result

def verify_full(probed, opener=None, *, max_bytes=MAX_COMPRESSED, max_rows=MAX_DATA_ROWS):
    result=dict(probed)
    if result.get("metadata_gate") != "METADATA_HEADER_SCREEN_PASS":
        result["full_file_gate"]="HOLD_HEADER_GATE_NOT_PASSED"
        return result
    expected_md5=result["yaml_selected_fields"]["data_file_md5sum"].lower()
    filename=result["raw_filename"]
    content_length=result.get("compressed_content_length")
    if content_length is not None and content_length > max_bytes:
        result["full_file_gate"]="HOLD_EXCEEDS_BOUNDED_BYTE_BUDGET"
        return result
    try:
        with tempfile.TemporaryDirectory(prefix="is-giga-public-") as tmp:
            dest=Path(tmp)/"source.tsv.gz"
            md5=hashlib.md5()
            sha=hashlib.sha256()
            total=0
            with open_checked(safe_asset_url(filename),opener=opener) as stream, dest.open("wb") as out:
                while True:
                    part=stream.read(1024*1024)
                    if not part:
                        break
                    total+=len(part)
                    if total>max_bytes:
                        result["full_file_gate"]="HOLD_EXCEEDS_BOUNDED_BYTE_BUDGET"
                        return result
                    md5.update(part);sha.update(part);out.write(part)
            result["full_compressed_bytes"]=total
            result["full_compressed_md5"]=md5.hexdigest()
            result["full_compressed_sha256"]=sha.hexdigest()
            if content_length is not None and total != content_length:
                result["full_file_gate"]="HOLD_TRUNCATED_TRANSFER"
                return result
            if md5.hexdigest()!=expected_md5:
                result["full_file_gate"]="HOLD_YAML_MD5_MISMATCH"
                return result
            rows=0
            qc=Counter()
            panel_hits={key:[] for key in CANONICAL_PANEL}
            with gzip.open(dest,"rt",encoding="utf-8-sig",newline="") as f:
                reader=csv.reader(f,delimiter="\t")
                first=next(reader)
                cols={name.lower():i for i,name in enumerate(first)}
                needed=("beta","standard_error","p_value","effect_allele",
                        "other_allele","effect_allele_frequency")
                if not set(needed).issubset(cols):
                    result["full_file_gate"]="HOLD_NUMERIC_QC_COLUMNS_MISSING"
                    return result
                for row in reader:
                    rows+=1
                    if rows>max_rows:
                        result["full_file_gate"]="HOLD_EXCEEDS_ROW_SCAN_BOUND"
                        return result
                    if len(row)!=len(first):
                        result["full_file_gate"]="HOLD_MALFORMED_GWAS_ROW"
                        return result
                    ea=row[cols["effect_allele"]].strip().upper()
                    oa=row[cols["other_allele"]].strip().upper()
                    if not (ea in ("A","C","G","T") and oa in ("A","C","G","T") and ea!=oa):
                        qc["non_acgt_or_equal_allele_rows"]+=1
                    elif {ea,oa} in ({"A","T"},{"C","G"}):
                        qc["palindromic_at_or_cg_rows"]+=1
                    try:
                        beta=float(row[cols["beta"]])
                        se=float(row[cols["standard_error"]])
                        p=float(row[cols["p_value"]])
                        eaf=float(row[cols["effect_allele_frequency"]])
                        if not (math.isfinite(beta) and math.isfinite(se) and se>0
                                and math.isfinite(p) and 0<=p<=1):
                            raise ValueError("invalid beta/se/p")
                        qc["finite_beta_se_p_positive_se_rows"]+=1
                        chrom=row[cols["chromosome"]].strip().lower().removeprefix("chr")
                        try:
                            locus=(chrom,int(row[cols["base_pair_location"]]))
                        except ValueError:
                            locus=None
                        if locus in CANONICAL_PANEL:
                            panel_hits[locus].append({
                                "allele_pair_match":frozenset((ea,oa))==CANONICAL_PANEL[locus]["alleles"],
                                "source_p":p,
                            })
                        if not math.isfinite(eaf) or not 0<=eaf<=1:
                            qc["invalid_eaf_rows"]+=1
                        else:
                            qc["valid_eaf_rows"]+=1
                        if rows%101==0:
                            z=abs(beta/se)
                            p_wald=math.erfc(z/math.sqrt(2))
                            if p>0 and p_wald>0:
                                qc["wald_two_sided_comparable_rows"]+=1
                                if abs(math.log10(p)-math.log10(p_wald))>0.2:
                                    qc["wald_p_log10_discordant_sample_rows"]+=1
                            else:
                                qc["wald_p_underflow_or_zero_sample_rows"]+=1
                    except (ValueError,OverflowError):
                        qc["missing_or_invalid_beta_se_p_or_eaf_rows"]+=1
            result["full_gwas_data_rows"]=rows
            result["gwas_numeric_qc"]=dict(qc)
            panel=[]
            for key, anchor in CANONICAL_PANEL.items():
                hits=panel_hits[key]
                matching=(len(hits)==1 and hits[0]["allele_pair_match"] and
                          hits[0]["source_p"]>0 and
                          abs(math.log10(hits[0]["source_p"])-math.log10(anchor["giga_ais_p"]))<=0.05)
                panel.append({
                    "gene":anchor["gene"], "chromosome_GRCh37":key[0],
                    "base_pair_location_GRCh37":key[1],
                    "candidate_alleles_unordered":"".join(sorted(anchor["alleles"])),
                    "published_crosswalk_expected_p":anchor["giga_ais_p"],
                    "source_rows_at_exact_position":len(hits),
                    "source_original_p":hits[0]["source_p"] if len(hits)==1 else None,
                    "unordered_alleles_and_p_match":bool(matching),
                })
            result["original_to_drive_three_anchor_panel"]=panel
            result["original_to_drive_panel_gate"]=(
                "THREE_POSITION_ALLELE_AND_P_SCREEN_PASS_NOT_FULL_CANONICAL_VERIFICATION"
                if all(x["unordered_alleles_and_p_match"] for x in panel)
                else "HOLD_ORIGINAL_TO_DRIVE_3_POSITION_ANCHOR_MISMATCH")
            result["numeric_qc_status"]=(
                "BASIC_NUMERIC_FIELDS_SCREEN_PASS" if
                qc["finite_beta_se_p_positive_se_rows"]==rows
                and qc["non_acgt_or_equal_allele_rows"]==0
                and qc["invalid_eaf_rows"]==0
                else "HOLD_NUMERIC_OR_ALLELE_QC_REVIEW")
            result["original_file_byte_hash_verified"]=True
            result["full_file_gate"]="ORIGINAL_FILE_MD5_AND_GZIP_ROWS_PASS_REVIEW_REQUIRED"
    except (OSError,ValueError,EOFError,csv.Error,StopIteration) as e:
        result["full_file_gate"]="HOLD_FULL_DOWNLOAD_OR_GZIP_ERROR"
        result["full_failure_type"]=type(e).__name__
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--mode",choices=["probe","full"],default="probe")
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    result=probe()
    if args.mode=="full":
        result=verify_full(result)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("IS_GIGA_AIS_FTP="+json.dumps(result,sort_keys=True))
if __name__=="__main__":
    main()
