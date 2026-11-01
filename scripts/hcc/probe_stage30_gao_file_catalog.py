#!/usr/bin/env python3
"""Stage30 read-only PDC study file catalog + candidate quantitation inventory.

This never downloads matrix files or exports sample/file names, file IDs,
patient IDs, signed URLs, or any per-case data. Filenames are categorized
in-memory. Catalog candidates DO NOT prove quantitative matrix availability.
"""
import argparse
import collections
import datetime as dt
import json
import pathlib
import re
from probe_stage30_gao_catalog import PDC,STUDY_FIELDS,pdc_query

def probe(pid):
    out={"pdc_study_id":pid,"status":"ACCESS_HOLD",
         "PDC_study_uuid_not_exported":True,
         "source_file_names_and_identifiers_not_exported":True,
         "actual_gene_by_patient_quantitation_matrix_downloaded":False}
    try:
        data,_=pdc_query('{ study(pdc_study_id:"'+pid+'"){ '+STUDY_FIELDS+' } }')
        item=data.get("study")
        if isinstance(item,list):
            if len(item)!=1:raise ValueError("PDC ambiguous study")
            item=item[0]
        if not isinstance(item,dict) or item.get("pdc_study_id")!=pid:
            raise ValueError("no PDC study")
        study_id=item.get("study_id")
        if not isinstance(study_id,str) or not re.fullmatch(r"[0-9a-fA-F-]{32,36}",study_id):
            raise ValueError("study UUID unavailable")
        q='{ getPaginatedFiles(study_id:"'+study_id+'",offset:0,limit:40) { total files { file_name file_type file_format data_category md5sum } pagination {count} } }'
        d,digest=pdc_query(q)
        v=d.get("getPaginatedFiles")
        if not isinstance(v,dict) or not isinstance(v.get("files"),list):
            raise ValueError("PDC paginated files API not accessible")
        total=v.get("total")
        if not isinstance(total,int) or total<0:
            raise ValueError("missing PDC file count")
        rows=v["files"]
        if len(rows)>40:raise ValueError("file preview limit exceeded")
        categories=collections.Counter()
        formats=collections.Counter()
        filetypes=collections.Counter()
        probable=0
        for x in rows:
            if not isinstance(x,dict):continue
            cat=str(x.get("data_category") or "")[:100]
            fmt=str(x.get("file_format") or "")[:80]
            typ=str(x.get("file_type") or "")[:100]
            categories[cat]+=1
            formats[fmt]+=1
            filetypes[typ]+=1
            # Only a lexical signal on first page of public file metadata.
            # This is not a quantitative matrix read/validation.
            name=str(x.get("file_name") or "")
            if re.search(r"protein|quant|abundan|matrix|ratio|normalized|tmt|report|txt|tsv",name,re.I):
                probable+=1
        out.update({"status":"VERIFIED_PDC_FILE_PREVIEW",
            "source_response_sha256":digest,
            "study_case_count":int(item.get("cases_count") or 0),
            "total_file_records_reported":total,
            "preview_file_n":len(rows),
            "preview_data_categories":dict(sorted(categories.items())),
            "preview_formats":dict(sorted(formats.items())),
            "preview_file_types":dict(sorted(filetypes.items())),
            "preview_filename_lexical_quant_candidate_n":probable,
            "patient_aliquot_biospecimen_manifest_join_verified":False,
            "matrix_numerical_cells_verified":False})
    except Exception as e:
        out["error_class"]=type(e).__name__
        # Never print error bodies as they might contain public filenames, URLs,
        # submitted query values or retrievable clinical identifiers.
    return out

def verify(obj):
    if obj.get("schema")!="hcc-stage30-pdc-file-catalog-preview-v1":
        raise ValueError("invalid preview schema")
    p=obj.get("pdc_studies")
    if not isinstance(p,list) or {x.get("pdc_study_id") for x in p}!=set(PDC):
        raise ValueError("unexpected PDC target studies")
    for x in p:
        for k in ("PDC_study_uuid_not_exported",
                  "source_file_names_and_identifiers_not_exported"):
            if x.get(k) is not True:raise ValueError("public identifier leak risk")
        for k in ("actual_gene_by_patient_quantitation_matrix_downloaded",
                  "patient_aliquot_biospecimen_manifest_join_verified",
                  "matrix_numerical_cells_verified"):
            if x.get(k) is True:raise ValueError("unexpected quantitative claim")
        if x.get("status")=="VERIFIED_PDC_FILE_PREVIEW":
            if not 0<=x["preview_file_n"]<=40 or x["preview_file_n"]>x["total_file_records_reported"]:
                raise ValueError("invalid file preview counts")
        elif x.get("status")!="ACCESS_HOLD":raise ValueError("invalid PDC access status")
    return {"schema":"hcc-stage30-pdc-file-catalog-review-v1",
        "verified_study_file_catalog_n":sum(x["status"]=="VERIFIED_PDC_FILE_PREVIEW" for x in p),
        "G30_B":"HOLD_QUANT_MATRIX_NOT_DOWNLOADED",
        "G30_D":"STOP_NO_CASE_SAMPLE_ALIQUOT_JOIN",
        "independent_biological_replication_completed":False,
        "numeric_matrix_success_is_not_implied_by_a_filename":True}

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--out",required=True);a.add_argument("--review-out",required=True)
    args=a.parse_args()
    v={"schema":"hcc-stage30-pdc-file-catalog-preview-v1",
       "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
       "pdc_studies":[probe(s) for s in PDC],
       "file_page_offset":0,"file_page_max_limit":40,
       "no_quantitative_data_downloaded":True}
    out=pathlib.Path(args.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(v,indent=2)+"\n")
    reviewed=verify(v)
    p=pathlib.Path(args.review_out);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(reviewed,indent=2)+"\n")
    print(json.dumps(reviewed,indent=2))
    print("PDC_FILE_PREVIEW",json.dumps({x["pdc_study_id"]:{
        k:x.get(k) for k in ("status","total_file_records_reported","preview_file_n",
                              "preview_data_categories","preview_formats","preview_file_types",
                              "preview_filename_lexical_quant_candidate_n","error_class")
        } for x in v["pdc_studies"]},separators=(",",":")))
if __name__=="__main__":main()
