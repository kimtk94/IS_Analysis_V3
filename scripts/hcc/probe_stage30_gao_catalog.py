#!/usr/bin/env python3
"""Stage30 read-only Gao 2019 PDC catalog and NODE landing-access probe.

No PDC/NODE per-case, patient or sample identifiers, signed file URLs,
downloadable molecular matrices or raw response bodies are saved or printed.
A successful PDC STUDY catalog check is NOT matched-RNA/proteome access.
"""
import argparse
import datetime as dt
import hashlib
import json
import pathlib
import re
import urllib.request
import urllib.error

GRAPHQL="https://pdc.cancer.gov/graphql"
NODE="https://www.biosino.org/node/project/detail/OEP000321"
PDC=("PDC000198","PDC000199")
LIMIT=1024*1024
AGENT="HCC-Stage30-Gao-readonly-metadata-QC/1.0"
STUDY_FIELDS="study_id pdc_study_id study_name disease_type experiment_type cases_count aliquots_count filesCount { data_category file_type files_count }"
CATALOG_FIELDS="pdc_study_id versions { study_id study_version is_latest_version }"

def http(url,payload=None,limit=LIMIT,timeout=28):
    headers={"User-Agent":AGENT,"Accept":"application/json"}
    if payload is not None:
        headers["Content-Type"]="application/json"
        payload=json.dumps({"query":payload},separators=(",",":")).encode()
    req=urllib.request.Request(url,data=payload,headers=headers)
    with urllib.request.urlopen(req,timeout=timeout) as rsp:
        raw=rsp.read(limit+1)
        if len(raw)>limit:raise ValueError("remote response exceeds hard limit")
        return {"bytes":raw,"http_status":rsp.status,
                "response_sha256":hashlib.sha256(raw).hexdigest(),
                "media_type":rsp.headers.get("Content-Type","")[:70]}

def pdc_query(q):
    result=http(GRAPHQL,q)
    obj=json.loads(result["bytes"])
    if not isinstance(obj,dict) or obj.get("errors") or not isinstance(obj.get("data"),dict):
        raise ValueError("PDC GraphQL returned errors or non-data response")
    return obj["data"],result["response_sha256"]

def probe_study(pid):
    if pid not in PDC:raise ValueError("unapproved study")
    try:
        data,digest=pdc_query('{ study (pdc_study_id: "'+pid+'") { '+STUDY_FIELDS+' } }')
        val=data.get("study")
        if isinstance(val,list):
            if len(val)!=1:raise ValueError("ambiguous PDC study")
            val=val[0]
        if not isinstance(val,dict) or val.get("pdc_study_id")!=pid:
            raise ValueError("returned wrong PDC study")
        for key in ("cases_count","aliquots_count"):
            n=val.get(key)
            if type(n) not in (int,float) or n<0:raise ValueError("missing case or aliquot count")
        cats=val.get("filesCount") or []
        if not isinstance(cats,list):raise ValueError("unexpected file counts")
        groups=[]
        for row in cats:
            n=row.get("files_count")
            if type(n) not in (int,float) or n<0:raise ValueError("bad file group")
            groups.append({"category":str(row.get("data_category",""))[:100],
                "type":str(row.get("file_type",""))[:120],"n":int(n)})
        return {
            "status":"VERIFIED_STUDY_METADATA","pdc_study_id":pid,
            "study_version_uuid_present":bool(val.get("study_id")),
            "study_name":str(val.get("study_name") or "")[:180],
            "disease_type":str(val.get("disease_type") or "")[:130],
            "experiment_type":str(val.get("experiment_type") or "")[:120],
            "published_study_case_count":int(val["cases_count"]),
            "published_study_aliquot_count":int(val["aliquots_count"]),
            "study_metadata_response_sha256":digest,
            "file_category_counts":groups,
            "full_quantitation_file_downloaded":False,
            "matched_RNA_protein_patient_count_verified":False}
    except (urllib.error.URLError,TimeoutError,ValueError,KeyError,json.JSONDecodeError,OSError) as exc:
        kind=type(exc).__name__
        if isinstance(exc,urllib.error.HTTPError):
            kind="HTTP_"+str(exc.code)
        return {
            "status":"ACCESS_HOLD","pdc_study_id":pid,
            "failure_category":kind,
            "full_quantitation_file_downloaded":False,
            "matched_RNA_protein_patient_count_verified":False
        }

def probe_catalog():
    try:
        data,digest=pdc_query("{ studyCatalog { "+CATALOG_FIELDS+" } }")
        catalog=data.get("studyCatalog")
        if not isinstance(catalog,list):raise ValueError("catalog schema changed")
        summary={}
        for item in catalog:
            if item.get("pdc_study_id") in PDC:
                ver=item.get("versions") or []
                summary[item["pdc_study_id"]]={
                    "version_count":len(ver),
                    "latest_version_count":sum(x.get("is_latest_version") is True for x in ver)}
        return {"status":"VERIFIED_STUDY_CATALOG","sha256":digest,
                "requested_study_ids_found":summary,
                "raw_study_catalog_not_exported":True}
    except (urllib.error.URLError,TimeoutError,ValueError,KeyError,json.JSONDecodeError,OSError) as exc:
        return {"status":"ACCESS_HOLD","failure_category":type(exc).__name__,
                "raw_study_catalog_not_exported":True}

def probe_node():
    try:
        rsp=http(NODE,limit=220000,timeout=24)
        body=rsp["bytes"].decode("utf-8","replace")
        # Public landing page may be a JS shell, not a data-record response.
        mentions_id="OEP000321" in body
        has_hcc=bool(re.search(r"hepatocellular|HCC|liver.cancer",body,re.I))
        return {"status":"LANDING_PAGE_RESPONDED",
            "http_status":rsp["http_status"],"response_sha256":rsp["response_sha256"],
            "page_mentions_exact_project_id":mentions_id,
            "page_mentions_liver_cancer":has_hcc,
            "actual_RNA_quantitation_file_verified":False,
            "source_patient_identifier_bridge_verified":False}
    except (urllib.error.URLError,TimeoutError,ValueError,KeyError,json.JSONDecodeError,OSError) as exc:
        return {"status":"ACCESS_HOLD","failure_category":type(exc).__name__,
                "actual_RNA_quantitation_file_verified":False,
                "source_patient_identifier_bridge_verified":False}

def review(v):
    if v.get("schema")!="hcc-stage30-gao-pdc-node-source-catalog-v1" or len(v.get("pdc_studies",[]))!=2:
        raise ValueError("Stage30 immutable target contract invalid")
    if {x.get("pdc_study_id") for x in v["pdc_studies"]}!=set(PDC):
        raise ValueError("unexpected or missing PDC studies")
    if v.get("node_project_accession")!="OEP000321":
        raise ValueError("NODE identifier drift")
    for s in v["pdc_studies"]:
        for k in ("full_quantitation_file_downloaded","matched_RNA_protein_patient_count_verified"):
            if s.get(k) is not False:raise ValueError("unsupported quantitative evidence")
        if s.get("status") not in ("VERIFIED_STUDY_METADATA","ACCESS_HOLD"):
            raise ValueError("unsupported source catalog claim")
    if v["node"].get("actual_RNA_quantitation_file_verified") is not False or v["node"].get("source_patient_identifier_bridge_verified") is not False:
        raise ValueError("unsupported NODE source molecular access")
    if v.get("case_sample_aliquot_cross_resource_join_completed") is not False or v.get("independent_cohort_gene_effect_replication_completed") is not False:
        raise ValueError("premature biological replication claim")
    return {"schema":"hcc-stage30-gao-catalog-postflight-v1",
        "PDC_study_metadata_verified_n":sum(x["status"]=="VERIFIED_STUDY_METADATA" for x in v["pdc_studies"]),
        "PDC_study_metadata_access_hold_n":sum(x["status"]=="ACCESS_HOLD" for x in v["pdc_studies"]),
        "NODE_landing_page_responded":v["node"].get("status")=="LANDING_PAGE_RESPONDED",
        "source_198_and_199_case_counts_are_not_cohort_join_counts":True,
        "G30_A":"SOURCE_CATALOG_PARTIAL" if all(x["status"]=="VERIFIED_STUDY_METADATA" for x in v["pdc_studies"]) else "ACCESS_HOLD",
        "G30_B":"HOLD_NO_QUANT_MATRIX",
        "G30_C":"HOLD_NO_RNA_MATRIX",
        "G30_D":"STOP_NO_VALIDATED_MATCHED_PATIENT_JOIN",
        "G30_E":"NOT_STARTED",
        "independent_HCC_biological_replication_claim_allowed":False}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True)
    ap.add_argument("--review-out",required=True)
    args=ap.parse_args()
    d={"schema":"hcc-stage30-gao-pdc-node-source-catalog-v1",
       "checked_at_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
       "source_type":"read_only_official_API_and_landing_page",
       "pdc_studies":[probe_study(x) for x in PDC],
       "study_catalog":probe_catalog(),
       "node_project_accession":"OEP000321","node":probe_node(),
       "case_sample_aliquot_cross_resource_join_completed":False,
       "independent_cohort_gene_effect_replication_completed":False,
       "patient_aliquot_identifiers_or_signed_urls_exported":False}
    p=pathlib.Path(args.out);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(d,indent=2)+"\n",encoding="utf-8")
    q=review(d)
    out=pathlib.Path(args.review_out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(q,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(q,indent=2))
    print("PDC_STUDY_STATUS_COUNTS",json.dumps({x["pdc_study_id"]:{"status":x["status"],
      "cases":x.get("published_study_case_count"),"aliquots":x.get("published_study_aliquot_count"),
      "categories":x.get("file_category_counts",[])} for x in d["pdc_studies"]}))
    print("NODE_LANDING",json.dumps(d["node"]))
if __name__=="__main__":main()
