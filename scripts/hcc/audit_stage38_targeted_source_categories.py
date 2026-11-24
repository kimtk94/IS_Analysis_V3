#!/usr/bin/env python3
"""Stage38: targeted PDC per-category author source discovery.

Global filesPerStudy pagination is heavily occupied by raw/processed spectra
(>=600 records per Gao study). Query source-relevant categories directly and
confirm the category filter is honored. Never interpret an empty result as
proof of absence without independent official filesCount inventory.
"""
import argparse
import collections
import datetime as dt
import json
import pathlib
import re
from probe_stage30_gao_catalog import PDC,pdc_query
from audit_stage37_all_pdc_file_categories import filename_flags

CATEGORIES=("Alternate Processing Pipeline","Other Metadata","Protein Assembly")
PAGE_LIMIT=75
MAX_PAGES=3

def declared_category_counts(pid,fetch=pdc_query):
    query=('{study(pdc_study_id:'+json.dumps(pid)+
           '){pdc_study_id filesCount {data_category file_type files_count}}}')
    result,digest=fetch(query)
    entries=result.get("study")
    if isinstance(entries,dict):entries=[entries]
    if not isinstance(entries,list) or len(entries)!=1 or entries[0].get("pdc_study_id")!=pid:
        raise ValueError("official study identity ambiguous or unavailable")
    counts=entries[0].get("filesCount")
    if not isinstance(counts,list):raise ValueError("missing official category count")
    per=collections.Counter()
    for v in counts:
        if not isinstance(v,dict) or not isinstance(v.get("files_count"),int):
            raise ValueError("invalid official file category count")
        if v["files_count"]<0:raise ValueError("negative category count")
        per[str(v.get("data_category") or "Unspecified")]+=v["files_count"]
    return dict(per),digest

def query_category(pid,cat,fetch=pdc_query):
    if pid not in PDC or cat not in CATEGORIES:
        raise ValueError("unapproved source study/category")
    unique=set()
    counts=collections.Counter()
    fmts=collections.Counter()
    pages=[]
    out={"category":cat,"status":"ACCESS_HOLD",
         "original_file_names_urls_or_ids_exported":False,
         "same_biopsy_RNA_protein_join_verified":False}
    try:
        for i in range(MAX_PAGES):
            q=('{filesPerStudy(pdc_study_id:'+json.dumps(pid)+
                ',data_category:'+json.dumps(cat)+
                ',offset:'+str(i*PAGE_LIMIT)+',limit:'+str(PAGE_LIMIT)+')'
                '{file_id file_name data_category file_format md5sum}}')
            obj,digest=fetch(q)
            rows=obj.get("filesPerStudy")
            if not isinstance(rows,list) or len(rows)>PAGE_LIMIT:
                raise ValueError("PDC category page invalid")
            for f in rows:
                if not isinstance(f,dict) or f.get("data_category")!=cat:
                    raise ValueError("official PDC data category filter ignored")
                fid=f.get("file_id")
                if not isinstance(fid,str) or not fid or fid in unique:
                    raise ValueError("duplicate or missing official PDC file ID")
                unique.add(fid)
                fmts[str(f.get("file_format") or "Unspecified")[:60]]+=1
                counts.update(filename_flags(f.get("file_name")))
            pages.append({"items_n":len(rows),"source_response_sha256":digest})
            if len(rows)<PAGE_LIMIT:break
        end=bool(pages and pages[-1]["items_n"]<PAGE_LIMIT)
        out.update({"status":"CATEGORY_SOURCE_LIST_CHECKED" if end else "TRUNCATED_CATEGORY_SOURCE_LIST",
                    "category_page_exhausted":end,
                    "candidate_source_files_n":len(unique),
                    "category_file_format_counts":dict(sorted(fmts.items())),
                    "candidate_filename_indicators_n":dict(sorted(counts.items())),
                    "category_pages":pages,
                    "source_alias_to_aliquot_330_verified":False})
    except Exception as e:
        out["error_class"]=type(e).__name__
    return out

def study(pid,fetch=pdc_query):
    out={"pdc_study_id":pid,"status":"TARGETED_CATEGORY_AUDIT",
         "patient_or_raw_source_names_exported":False,
         "data_join_complete":False,
         "official_study_counts_status":"ACCESS_HOLD"}
    try:
        declared,digest=declared_category_counts(pid,fetch)
        out.update({"official_study_counts_status":"SOURCE_COUNTS_VERIFIED",
                    "official_study_file_counts_source_SHA256":digest,
                    "official_declared_source_categories":declared})
    except Exception as e:
        out["official_study_file_counts_error_class"]=type(e).__name__
    out["categories"]=[query_category(pid,c,fetch) for c in CATEGORIES]
    if out["official_study_counts_status"]=="SOURCE_COUNTS_VERIFIED":
        out["category_count_reconciliation"]={
            v["category"]:{
                "official_declared_count_n":out["official_declared_source_categories"].get(v["category"],0),
                "api_source_list_n":v.get("candidate_source_files_n"),
                "count_matches_official_inventory":v.get("status")=="CATEGORY_SOURCE_LIST_CHECKED"
                    and v.get("candidate_source_files_n")==out["official_declared_source_categories"].get(v["category"],0)
            } for v in out["categories"]}
    return out

def review(d):
    if d.get("schema")!="hcc-stage38-targeted-PDC-author-category-v1":
        raise ValueError("schema mismatch")
    arr=d.get("studies")
    if not isinstance(arr,list) or {r.get("pdc_study_id") for r in arr}!=set(PDC):
        raise ValueError("unsupported studies")
    for v in arr:
        if v.get("patient_or_raw_source_names_exported") is not False or v.get("data_join_complete") is not False:
            raise ValueError("unsafe patient claim")
        cats=v.get("categories")
        if not isinstance(cats,list) or {x.get("category") for x in cats}!=set(CATEGORIES):
            raise ValueError("target categories omitted")
        for c in cats:
            if c.get("original_file_names_urls_or_ids_exported") is not False or c.get("same_biopsy_RNA_protein_join_verified") is not False:
                raise ValueError("unsafe source disclosure")
            if c["status"] not in ("CATEGORY_SOURCE_LIST_CHECKED","TRUNCATED_CATEGORY_SOURCE_LIST","ACCESS_HOLD"):
                raise ValueError("invalid source status")
            if c["status"]=="CATEGORY_SOURCE_LIST_CHECKED" and c.get("category_page_exhausted") is not True:
                raise ValueError("unproven complete category")
    return {"schema":"hcc-stage38-official-category-audit-postflight-v1",
        "studies_with_reconciled_PDC_file_category_counts_n":
            sum(v["official_study_counts_status"]=="SOURCE_COUNTS_VERIFIED"
                and all(t["count_matches_official_inventory"] for t in v.get("category_count_reconciliation",{}).values())
                for v in arr),
        "G30_B":"HOLD_UNTIL_ORIGINAL_330_QUANT_COLUMNS_LINK_OFFICIAL_ALIQUOTS",
        "G30_C":"HOLD_ORIGINAL_OEP000321_RNA_QUANTITATION",
        "G30_D":"STOP_NO_SAME_BIOPSY_TRANSCRIPT_PROTEIN_JOIN",
        "new_biological_significance_verified":False}

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True);p.add_argument("--review-out",required=True)
    a=p.parse_args()
    d={"schema":"hcc-stage38-targeted-PDC-author-category-v1",
       "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
       "studies":[study(p) for p in PDC]}
    v=review(d)
    for target,obj in ((a.out,d),(a.review_out,v)):
        file=pathlib.Path(target);file.parent.mkdir(parents=True,exist_ok=True)
        file.write_text(json.dumps(obj,indent=2)+"\n")
    print("HCC_STAGE38_CATEGORY_SOURCES",json.dumps(d))
    print("HCC_STAGE38_REVIEW",json.dumps(v))

if __name__=="__main__":main()
