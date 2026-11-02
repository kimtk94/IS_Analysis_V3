#!/usr/bin/env python3
"""Independent official source gates: GIGASTROKE EAS AIS N and GTEx QTL release lineage.

Read-only bounded HTTPS source metadata. Does not authorize coloc/MR, infer
per-SNP N, claim original GTEx v8 equivalence, or download full QTL/GWAS files.
"""
import argparse
import hashlib
from html.parser import HTMLParser
import json
import math
import re
import urllib.parse
import urllib.request
from pathlib import Path

ACCESSION="GCST90104545"
ORIGINAL_GWAS_MD5="904e12e33ec3c2eeefe4dc93bec4995c"
YAML_URL=("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
          "GCST90104001-GCST90105000/GCST90104545/"
          "GCST90104545_buildGRCh37.tsv.gz-meta.yaml")
QTL_URL=("https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/"
         "sumstats/QTS000015/QTD000161/QTD000161.all.tsv.gz")
EQTL_STUDIES="https://www.ebi.ac.uk/eqtl/Studies/"
EQTL_RELEASE="https://www.ebi.ac.uk/eqtl/Release_notes/"
DOCUMENT_URLS=(YAML_URL,QTL_URL,EQTL_STUDIES,EQTL_RELEASE)
PUBLISHED_REFERENCE={"cases":19032,"controls":237242,"label":"GIGASTROKE East Asian Any Ischemic Stroke"}
MAX_METADATA=1_000_000
MAX_QTL_PREFIX=131072

class TextParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.text=[]
    def handle_data(self,data):
        self.text.append(data)
    def value(self):
        return " ".join(" ".join(self.text).split())

def metadata_fetch(url, *, fetcher=None):
    if url not in DOCUMENT_URLS:
        raise ValueError("Unapproved research source URL")
    if fetcher:
        return fetcher(url)
    prefix=(url==QTL_URL)
    req=urllib.request.Request(url,headers={
        "User-Agent":"IS-Research-QTL-and-N-Lineage/1.0",
        **({"Range":f"bytes=0-{MAX_QTL_PREFIX-1}"} if prefix else {})
    })
    with urllib.request.urlopen(req,timeout=22) as r:
        dest=urllib.parse.urlsplit(r.url)
        orig=urllib.parse.urlsplit(url)
        if dest.scheme!="https" or dest.hostname!=orig.hostname:
            raise ValueError("Source redirect crosses trusted hostname")
        maximum=MAX_QTL_PREFIX if prefix else MAX_METADATA
        data=r.read(maximum+1)
        if len(data)>maximum:
            raise ValueError("Official response exceeds bounded read")
        content_type=r.headers.get_content_type()
        if prefix:
            if not data.startswith(b"\x1f\x8b"):
                raise ValueError("QTL source prefix is not gzip")
        elif content_type not in ("text/plain","text/html","application/x-yaml",
                                  "text/x-yaml","application/octet-stream"):
            raise ValueError("Unexpected source metadata MIME")
        return data

def scalar_yaml(raw):
    # Strictly read only simple scalar fields, never infer unobserved nested values.
    obj={}
    for line in raw.decode("utf-8","replace").splitlines():
        match=re.fullmatch(r"([A-Za-z_][A-Za-z0-9_]*)\s*:\s*(.*?)\s*",line)
        if not match: continue
        k,v=match.groups()
        v=v.split(" #",1)[0].strip().strip('"').strip("'")
        if k in obj and obj[k]!=v:
            raise ValueError("Conflicting duplicate YAML metadata key")
        obj[k]=v
    return obj

def integer_metadata(val):
    if val is None:
        return None
    norm=str(val).strip().replace(",","")
    if not re.fullmatch(r"[0-9]+",norm):return None
    return int(norm)

def inspect_gwas_yaml(raw):
    metadata=scalar_yaml(raw)
    fields={key:metadata.get(key) for key in (
        "gwas_id","data_file_name","data_file_md5sum","genome_assembly",
        "file_type","is_harmonised","case_control_study","sample_size",
        "case_count","control_count","sample_ancestry","sample_ancestry_category",
        "trait_description") if key in metadata}
    same_source=(
        metadata.get("gwas_id")==ACCESSION
        and metadata.get("data_file_name")=="GCST90104545_buildGRCh37.tsv.gz"
        and metadata.get("data_file_md5sum","").lower()==ORIGINAL_GWAS_MD5
        and metadata.get("genome_assembly")=="GRCh37")
    case=integer_metadata(metadata.get("case_count"))
    control=integer_metadata(metadata.get("control_count"))
    total=integer_metadata(metadata.get("sample_size"))
    n_consistent=(
        case is not None and control is not None and case>0 and control>0
        and (total is None or total==case+control))
    # In a GWAS meta-analysis pooled N can vary by SNP. Matching public N
    # cannot by itself make a per-SNP N or case fraction usable for coloc.
    fraction=(case/(case+control) if n_consistent else None)
    eff=(4*case*control/(case+control) if n_consistent else None)
    return {
        "source_file_identity_and_yaml_digest":"PASS" if same_source else "HOLD",
        "source_yaml_sha256":hashlib.sha256(raw).hexdigest(),
        "metadata_fields":fields,
        "file_specific_case_control_fields_available":n_consistent,
        "yaml_case_count":case,"yaml_control_count":control,"yaml_sample_size":total,
        "yaml_internal_N_consistency":n_consistent,
        "published_reference_case_count_only":PUBLISHED_REFERENCE["cases"],
        "published_reference_control_count_only":PUBLISHED_REFERENCE["controls"],
        "published_reference_matches_yaml":(
            case==PUBLISHED_REFERENCE["cases"] and control==PUBLISHED_REFERENCE["controls"]
            if n_consistent else None),
        "illustrative_s_from_yaml_not_adopted":fraction,
        "illustrative_Neff_from_yaml_not_adopted":eff,
        "per_variant_N_available_in_original_11_col_header":False,
        "association_level_N_or_case_fraction_approved":False,
        "coloc_approved":False,
        "N_gate":("YAML_N_PRESENT_REQUIRES_VARIANT_AND_META_ANALYSIS_REVIEW" if
                  same_source and n_consistent else
                  "HOLD_MISSING_OR_INCONSISTENT_FILE_SPECIFIC_CASE_CONTROL_N"),
    }

def inspect_qtl_metadata(studies_raw,releases_raw,qtl_prefix=None):
    studies=TextParser()
    studies.feed(studies_raw.decode("utf-8","replace"))
    release=TextParser()
    release.feed(releases_raw.decode("utf-8","replace"))
    s,r=studies.value(),release.value()
    # This is documentation-based semantics, NOT a QTD raw-stat identity check.
    original_imported=("Imported studies" in s and "GTEx_v8" in s)
    reprocessed=("QTS000015" in s and "GTEx" in s)
    r6=("release 6" in r.lower() and ".all.tsv.gz" in r)
    r8=("r8 pre-release" in r.lower())
    v10_reuse=("GTEx_v10" in r and "ids" in r.lower())
    return {
        "official_studies_page_sha256":hashlib.sha256(studies_raw).hexdigest(),
        "official_release_page_sha256":hashlib.sha256(releases_raw).hexdigest(),
        "qtl_dataset_id":"QTD000161",
        "qtl_study_id":"QTS000015",
        "official_reprocessed_GTEx_listing":reprocessed,
        "official_separate_imported_GTEx_v8_listing":original_imported,
        "official_release6_complete_cis_documented":r6,
        "official_r8_pre_release_documented":r8,
        "GTEx_v10_ID_reuse_warning_observed":v10_reuse,
        "qtl_prefix_gzip_accessible":(
            qtl_prefix is not None and qtl_prefix.startswith(b"\x1f\x8b")),
        "qtl_prefix_sha256":hashlib.sha256(qtl_prefix).hexdigest() if qtl_prefix is not None else None,
        "historical_2023_eqtl_release6_file_bytes_sha256_authenticated":False,
        "current_QTD000161_is_original_GTExv8_FastQTL":False,
        "same_GTEx_v8_donor_model_and_covariates_verified":False,
        "original_GTExv8_cis_results_verified":False,
        "QTL_gate":("REPROCESSED_VS_IMPORTED_SOURCE_DISTINGUISHED_ORIGINAL_HOLD" if
                    original_imported and reprocessed and r6
                    else "HOLD_OFFICIAL_QTL_LINEAGE_INCOMPLETE"),
    }

def audit(fetcher=metadata_fetch):
    out={"schema":"is-eas-ais-N-and-GTEx-eQTL-lineage-v1",
         "GWAS_accession":ACCESSION,
         "QTL_study":"QTS000015","QTL_dataset":"QTD000161",
         "scientific_status":"HOLD","MR_performed":False,"coloc_performed":False}
    try:
        out["gwas_yaml"]=inspect_gwas_yaml(fetcher(YAML_URL))
    except (OSError,ValueError,UnicodeError) as exc:
        out["gwas_yaml"]={"N_gate":"HOLD_GWAS_YAML_UNAVAILABLE",
                          "source_file_identity_and_yaml_digest":"HOLD",
                          "failure_class":type(exc).__name__,
                          "association_level_N_or_case_fraction_approved":False}
    try:
        out["qtl_provenance"]=inspect_qtl_metadata(
            fetcher(EQTL_STUDIES),fetcher(EQTL_RELEASE))
        try:
            prefix=fetcher(QTL_URL)
            out["qtl_provenance"]["qtl_prefix_gzip_accessible"]=prefix.startswith(b"\x1f\x8b")
            out["qtl_provenance"]["qtl_prefix_sha256"]=hashlib.sha256(prefix).hexdigest()
        except (OSError,ValueError,UnicodeError) as exc:
            out["qtl_provenance"]["qtl_prefix_gate"]="HOLD_QTL_PREFIX_UNAVAILABLE"
            out["qtl_provenance"]["qtl_prefix_failure_class"]=type(exc).__name__
    except (OSError,ValueError,UnicodeError) as exc:
        out["qtl_provenance"]={"QTL_gate":"HOLD_OFFICIAL_QTL_METADATA_UNAVAILABLE",
                               "failure_class":type(exc).__name__,
                               "current_QTD000161_is_original_GTExv8_FastQTL":False}
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    result=audit()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("IS_EAS_QTL_N_LINEAGE="+json.dumps(result,sort_keys=True,separators=(",",":")))

if __name__=="__main__":
    main()
