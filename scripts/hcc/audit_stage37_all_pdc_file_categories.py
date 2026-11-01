#!/usr/bin/env python3
"""Stage37: enumerate ALL official PDC study file categories, not only CDAP.

Read-only, strictly bounded PDC GraphQL metadata scan for two Gao HCC studies.
No raw file names, PDC identifiers, sample labels, URLs, or protein values are
logged. Discovery != complete 330-aliquot assay provenance.
"""
import argparse
import collections
import datetime as dt
import hashlib
import json
import pathlib
import re
import time

from probe_stage30_gao_catalog import PDC, pdc_query

PAGE_SIZE=75
MAX_PAGES=8
MAX_RECORDS=PAGE_SIZE*MAX_PAGES
PREFERRED_CATEGORIES=("Alternate Processing Pipeline","Other Metadata","Protein Assembly")
ALIASES={
    "gct":re.compile(r"\.gct(?:\.gz)?$",re.I),
    "cdap_sample":re.compile(r"\.sample\.(?:txt|tsv)(?:\.gz)?$",re.I),
    "author_design":re.compile(r"(sample.?map|sample.?sheet|experiment.?design|plex.?map|tmt.?design|aliquot.?map)",re.I),
    "author_quant":re.compile(r"(protein|phospho|quant|abundan|matrix|expression|spectral)",re.I),
    "author_supplement":re.compile(r"(supplement|supp\.?|clinical|metadata|patient)",re.I),
}

def filename_flags(s):
    n=str(s or "")
    return {k:bool(v.search(n)) for k,v in ALIASES.items()}

def query_page(pid,offset,fetch=pdc_query):
    if pid not in PDC or offset not in range(0,MAX_RECORDS,PAGE_SIZE):
        raise ValueError("unapproved study or page")
    # Deliberately NO data_category filter. Official filesPerStudy permits
    # pdc_study_id and pagination, per PDC public API reference.
    q=('{filesPerStudy(pdc_study_id:'+json.dumps(pid)+
       ',offset:'+str(offset)+',limit:'+str(PAGE_SIZE)+')'
       '{file_id file_name data_category file_format file_type md5sum file_size}}')
    data,digest=fetch(q)
    files=data.get("filesPerStudy")
    if not isinstance(files,list) or len(files)>PAGE_SIZE:
        raise ValueError("source pagination schema drift")
    if any(not isinstance(x,dict) for x in files):
        raise ValueError("malformed official PDC file record")
    return files,digest

def source(pid,fetch=pdc_query):
    o={"pdc_study_id":pid,"status":"ACCESS_HOLD",
       "patient_aliquot_or_filename_exported":False,
       "signed_url_or_raw_quant_exported":False,
       "quant_330_aliquot_join_certified":False,
       "RNA_protein_pairing_certified":False}
    pages=[]
    records=[]
    try:
        for offset in range(0,MAX_RECORDS,PAGE_SIZE):
            fetched=None
            for attempt in range(2):
                try:
                    fetched,digest=query_page(pid,offset,fetch)
                    break
                except (TimeoutError,ConnectionError,OSError):
                    if attempt==0:time.sleep(1)
                    else:raise
            pages.append({"offset":offset,"items_n":len(fetched),
                          "response_sha256":digest})
            records.extend(fetched)
            if len(fetched)<PAGE_SIZE:
                break
        ids=[str(f.get("file_id") or "") for f in records]
        if not records or any(not i for i in ids) or len(ids)!=len(set(ids)):
            raise ValueError("official PDC missing or duplicate file identity")
        exhaustive=bool(pages and pages[-1]["items_n"]<PAGE_SIZE)
        if len(pages)==MAX_PAGES and not exhaustive:
            o["status"]="TRUNCATED_MAX_PAGES_REVIEW"
        else:
            o["status"]="PDC_ALL_FILE_CATEGORIES_ENUMERATED"
        cats=collections.defaultdict(list)
        for f in records:
            cat=str(f.get("data_category") or "Unspecified")
            cats[cat].append(f)
        agg={}
        for cat,fs in sorted(cats.items()):
            flags=collections.Counter()
            formats=collections.Counter()
            for f in fs:
                flags.update(filename_flags(f.get("file_name")))
                formats[str(f.get("file_format") or "Unspecified")[:50]]+=1
            agg[cat]={
                "files_n":len(fs),
                "format_counts":dict(sorted(formats.items())),
                "candidate_suffix_or_keyword_counts":dict(sorted(flags.items())),
                "source_file_md5_metadata_present_n":
                    sum(bool(re.fullmatch(r"[0-9a-fA-F]{32}",str(f.get("md5sum") or ""))) for f in fs),
                "GCT_original_values_and_headers_verified":False,
                "CDAP_alias_to_official_aliquot_verified":False,
            }
        o.update({
            "pages":pages,
            "total_distinct_official_files_n":len(records),
            "all_file_categories":agg,
            "data_category_filter_not_used":True,
            "all_pages_exhausted":exhaustive,
            "alternate_processing_pipeline_files_n":len(cats.get("Alternate Processing Pipeline",[])),
            "other_metadata_files_n":len(cats.get("Other Metadata",[])),
            "protein_assembly_files_n":len(cats.get("Protein Assembly",[])),
            "any_new_official_source_category_found":
                any(k not in PREFERRED_CATEGORIES for k in cats),
            "original_PDC_GCT_file_catalog_candidate_n":
                sum(filename_flags(x.get("file_name"))["gct"] for x in records),
            "published_quant_matrix_can_now_join_by_position":False,
        })
    except Exception as e:
        o["status"]="ACCESS_HOLD"
        o["error_class"]=type(e).__name__
        # No ambiguous partial file inventory is declared exhaustive.
        o["all_pages_exhausted"]=False
    return o

def verify(doc):
    if doc.get("schema")!="hcc-stage37-all-PDC-source-file-categories-v1":
        raise ValueError("source schema")
    s=doc.get("studies")
    if not isinstance(s,list) or {x.get("pdc_study_id") for x in s}!=set(PDC):
        raise ValueError("study coverage")
    for a in s:
        for k in ("patient_aliquot_or_filename_exported","signed_url_or_raw_quant_exported",
                  "quant_330_aliquot_join_certified","RNA_protein_pairing_certified"):
            if a.get(k) is not False:raise ValueError("unsafe source inference")
        if a["status"]=="PDC_ALL_FILE_CATEGORIES_ENUMERATED":
            if a.get("all_pages_exhausted") is not True or not a.get("data_category_filter_not_used"):
                raise ValueError("unverified category coverage")
            if a["total_distinct_official_files_n"]!=sum(
                    f["files_n"] for f in a["all_file_categories"].values()):
                raise ValueError("source count drift")
            if a["total_distinct_official_files_n"]>MAX_RECORDS:
                raise ValueError("source bound exceeded")
        elif a["status"] not in ("ACCESS_HOLD","TRUNCATED_MAX_PAGES_REVIEW"):
            raise ValueError("invalid status")
        if a.get("published_quant_matrix_can_now_join_by_position") is True:
            raise ValueError("no position-based specimen identity")
    return {
        "schema":"hcc-stage37-source-category-postflight-v1",
        "complete_source_inventory_studies_n":
            sum(a["status"]=="PDC_ALL_FILE_CATEGORIES_ENUMERATED" for a in s),
        "studies_with_alternate_processing_category_n":
            sum(a.get("alternate_processing_pipeline_files_n",0)>0 for a in s),
        "G30_B":"HOLD_UNTIL_330_SOURCE_QUANT_COLUMN_TO_ALIQUOT_ID_PROVEN",
        "G30_C":"HOLD_NODE_OEP000321_ORIGINAL_RNA_MATRIX",
        "G30_D":"STOP_NO_INDEPENDENT_SAME_BIOPSY_RNA_PROTEIN",
        "biological_inference_allowed":False,
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True)
    ap.add_argument("--review-out",required=True)
    a=ap.parse_args()
    doc={"schema":"hcc-stage37-all-PDC-source-file-categories-v1",
         "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
         "studies":[source(pid) for pid in PDC]}
    v=verify(doc)
    for fn,x in ((a.out,doc),(a.review_out,v)):
        p=pathlib.Path(fn);p.parent.mkdir(parents=True,exist_ok=True)
        p.write_text(json.dumps(x,indent=2)+"\n",encoding="utf-8")
    print("HCC_STAGE37_CATALOG",json.dumps(doc))
    print("HCC_STAGE37_REVIEW",json.dumps(v))

if __name__=="__main__":main()
