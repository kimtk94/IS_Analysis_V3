#!/usr/bin/env python3
"""Check alternate PDC file categories for source GCT / sample mapping files.

Only anonymous category counts are output; no sample ID or filename export.
The first 25 files in each available category are a preview, not exhaustive.
"""
import argparse,datetime as dt,json,pathlib,re
from probe_stage30_gao_catalog import PDC,pdc_query,STUDY_FIELDS

def lookup(pid):
    out={"pdc_study_id":pid,"status":"SOURCE_ACCESS_HOLD",
         "file_names_exported":False,"sample_ids_exported":False}
    try:
        d,_=pdc_query('{study(pdc_study_id:"'+pid+'"){'+STUDY_FIELDS+'}}')
        st=d.get("study")
        if isinstance(st,list):st=st[0] if len(st)==1 else None
        if not isinstance(st,dict) or st.get("pdc_study_id")!=pid:
            raise ValueError("official PDC study unavailable")
        uuid=st.get("study_id")
        if not isinstance(uuid,str) or not re.fullmatch(r"[0-9a-f-]{36}",uuid,re.I):
            raise ValueError("PDC version invalid")
        cats={str(row["data_category"]):int(row["files_count"])
              for row in st.get("filesCount",[]) if row.get("data_category")}
        results={}
        target=[c for c in cats if any(q in c.lower() for q in
               ("supplement","alternate","metadata","assembly"))]
        for cat in target[:5]:
            try:
                query='{filesPerStudy(study_id:"'+uuid+'",data_category:"'+cat+'",offset:0,limit:25){file_name data_category}}'
                r,sha=pdc_query(query)
                v=r.get("filesPerStudy")
                if not isinstance(v,list):raise ValueError("PDC files missing")
                names=[str(x.get("file_name") or "") for x in v]
                results[cat]={
                  "status":"SOURCE_FILE_PREVIEW_VERIFIED",
                  "official_category_total_n":cats[cat],
                  "first_page_file_n":len(v),
                  "file_listing_may_be_truncated":cats[cat]>len(v),
                  "gct_candidate_name_n":sum(bool(re.search(r"\.gct(\.gz)?$",n,re.I)) for n in names),
                  "sample_txt_candidate_name_n":sum(bool(re.search(r"\.sample\.txt(\.gz)?$",n,re.I)) for n in names),
                  "sample_design_name_n":sum(bool(re.search(r"sample.?map|sample.?sheet|experimental.?design",n,re.I)) for n in names),
                  "file_metadata_sha256":sha}
            except Exception as e:
                results[cat]={"status":"SOURCE_FILE_ACCESS_HOLD",
                              "official_category_total_n":cats[cat],
                              "error_kind":type(e).__name__}
        out.update({"status":"FILE_CATEGORY_PREVIEW_COMPLETED",
                    "source_categories_considered_n":len(target),
                    "source_category_previews":results,
                    "source_GCT_content_or_patient_link_verified":False})
    except Exception as e:
        out["error_kind"]=type(e).__name__
    return out

def review(data):
    if data.get("schema")!="hcc-stage31-PDC-other-file-category-v1":
        raise ValueError("schema changed")
    studies=data.get("studies",[])
    if {x.get("pdc_study_id") for x in studies}!=set(PDC):
        raise ValueError("study IDs drift")
    if any(x.get("file_names_exported") is not False or x.get("sample_ids_exported") is not False for x in studies):
        raise ValueError("identifier exposure")
    return {"schema":"hcc-stage31-PDC-other-file-postflight-v1",
      "PDC_source_category_previews_completed_n":sum(x["status"]=="FILE_CATEGORY_PREVIEW_COMPLETED" for x in studies),
      "G30_B":"HOLD_GCT_OR_SAMPLE_DESIGN_CONTENT_AND_EXACT_ALIQUOT_JOIN",
      "G30_C":"HOLD_RNA_NODE",
      "G30_D":"STOP_NO_VALIDATED_EXTERNAL_RNA_PROTEIN_CASES"}

if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--out",required=True);a.add_argument("--review-out",required=True)
    x=a.parse_args()
    o={"schema":"hcc-stage31-PDC-other-file-category-v1",
       "time_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
       "studies":[lookup(p) for p in PDC]}
    f=pathlib.Path(x.out);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(o,indent=2)+"\n")
    v=review(o);g=pathlib.Path(x.review_out);g.parent.mkdir(parents=True,exist_ok=True)
    g.write_text(json.dumps(v,indent=2)+"\n")
    print(json.dumps(v,indent=2))
    print("STAGE31_FILE_CATEGORY_AGGREGATES",json.dumps(o["studies"]))
