#!/usr/bin/env python3
"""Stage54: strict source-verified NODE OEP000321 public data-id discovery.

Endpoints and JSON structures are from NIAID Data Ecosystem's open NODE crawler
(commit 8b0e9804afbf01af380d574c83b46a535d80a53f). Query exactly one
project; never enumerate the entire archive, guess an OED accession, invoke a
download URL or request private data. NIAID indexed metadata is a second,
independent evidence channel and is NOT substituted for NODE source authority.
No patient identifiers, file names except expected triplet, or clinical values
are emitted. Only explicitly public OED data accessions may appear in output.
"""
import argparse
import datetime as dt
import hashlib
import json
import pathlib
import re
import urllib.error
import urllib.parse
import urllib.request

PROJECT="OEP000321"
NODE_ORIGIN="https://www.biosino.org"
NIAID_ORIGIN="https://api.data.niaid.nih.gov"
SEARCH=f"{NODE_ORIGIN}/node/api/app/browse/search"
DATA_LIST=f"{NODE_ORIGIN}/node/api/app/browseDetail/getDataList"
NIAID_QUERY=NIAID_ORIGIN+"/v1/query"
TRIPLET=("HCCDB25_mRNA_level3.txt","HCCDB25.sample.txt","HCCDB25.patient.txt")
MAX_BYTES=2*1024*1024
MAX_SEARCH=15
PAGE_LIMIT=3
PAGESIZE=100
MAX_PUBLIC_REPORT=25
OED=re.compile(r"^OED[0-9]{5,}$",re.I)
INTERNAL=re.compile(r"^[A-Za-z0-9_-]{1,45}$")

class SourceError(ValueError):
    pass

def search_payload():
    return {"queryWord":PROJECT,"advQueryWord":"","sortKey":"modifiedDate","sortType":"desc",
            "pageNum":1,"pageSize":20,
            "leftStatQueries":[{"type":"Level","name":"project","treeRefId":"0_lsTree__Level"}]}

def source_request(url,method="GET",payload=None,opener=urllib.request.urlopen):
    p=urllib.parse.urlsplit(url)
    allowed=(p.scheme=="https" and p.hostname in ("www.biosino.org","api.data.niaid.nih.gov"))
    if not allowed or p.username or p.password or p.fragment:
        raise SourceError("SOURCE_ORIGIN_REJECTED")
    if p.hostname=="www.biosino.org":
        paths=( "/node/api/app/browse/search", "/node/api/app/browseDetail/getDataList" )
        general="/node/api/app/project/getGeneralInfo/"
        if not (p.path in paths or (p.path.startswith(general)
               and INTERNAL.fullmatch(p.path[len(general):]))):
            raise SourceError("NODE_ENDPOINT_NOT_CRAWLER_AUTHORIZED")
    else:
        if p.path!="/v1/query":
            raise SourceError("NIAID_ENDPOINT_NOT_DOCUMENTED")
    data=json.dumps(payload).encode("utf-8") if payload is not None else None
    if (method=="POST")!=(data is not None):
        raise SourceError("HTTP_METHOD_BODY_MISMATCH")
    req=urllib.request.Request(url,data=data,method=method,headers={
        "User-Agent":"HCC-Stage54-public-NODE-source-inventory/1.0",
        "Accept":"application/json,text/plain,*/*",
        "Content-Type":"application/json"})
    with opener(req,timeout=25) as r:
        if r.geturl()!=url or r.status!=200:
            raise SourceError("SOURCE_REDIRECT_OR_STATUS_UNEXPECTED")
        status=r.status
        raw=r.read(MAX_BYTES+1)
    if len(raw)>MAX_BYTES:
        raise SourceError("SOURCE_JSON_TOO_LARGE")
    try:
        data=json.loads(raw)
    except (ValueError,UnicodeError) as exc:
        raise SourceError("SOURCE_INVALID_JSON") from exc
    return data,{"http_status":status,"json_bytes":len(raw),
                 "response_sha256":hashlib.sha256(raw).hexdigest()}

def data_envelope(response):
    if not isinstance(response,dict) or "data" not in response or response["data"] is None:
        raise SourceError("NODE_RESPONSE_DATA_MISSING")
    return response["data"]

def search_candidates(payload):
    obj=data_envelope(payload)
    if not isinstance(obj,dict):
        raise SourceError("NODE_SEARCH_BAD_DATA")
    info=obj.get("pageInfo")
    if not isinstance(info,dict) or not isinstance(info.get("content"),list):
        raise SourceError("NODE_SEARCH_PAGE_INFO_MISSING")
    entries=info["content"]
    if len(entries)>100 or not all(isinstance(x,dict) for x in entries):
        raise SourceError("NODE_SEARCH_RESULT_CAP")
    ids=[]
    exact_hit=0
    for e in entries:
        if str(e.get("projectNo","")).upper()==PROJECT:
            exact_hit+=1
        candidate=e.get("id")
        if (isinstance(candidate,(int,str)) and not isinstance(candidate,bool)
            and INTERNAL.fullmatch(str(candidate))):
            candidate=str(candidate)
            if candidate not in ids:ids.append(candidate)
    return {"candidates":ids[:MAX_SEARCH],
            "listed_search_projects_n":len(entries),
            "explicit_accession_matches_in_search_n":exact_hit,
            "search_candidates_requiring_general_info_n":min(len(ids),MAX_SEARCH),
            "source_query_pages_reported":info.get("totalPages")
                if type(info.get("totalPages")) is int else None}

def project_general_url(internal):
    if not INTERNAL.fullmatch(str(internal)):
        raise SourceError("INVALID_INTERNAL_NODE_ID")
    return f"{NODE_ORIGIN}/node/api/app/project/getGeneralInfo/{internal}"

def datalist_payload(internal,page):
    if not INTERNAL.fullmatch(str(internal)) or page not in (1,2,3):
        raise SourceError("NON_SOURCE_PROJECT_OR_PAGE")
    return {"type":"project","typeNo":internal,"sortKey":"","sortType":"",
            "pageNum":page,"pageSize":PAGESIZE}

def parse_file_page(raw):
    obj=data_envelope(raw)
    if not isinstance(obj,dict) or type(obj.get("total")) is not int or obj["total"]<0:
        raise SourceError("NODE_DATALIST_TOTAL_UNVERIFIED")
    records=obj.get("list")
    if not isinstance(records,list) or len(records)>PAGESIZE:
        raise SourceError("NODE_DATALIST_RECORD_SCHEMA_UNVERIFIED")
    if not all(isinstance(r,dict) for r in records):
        raise SourceError("NODE_DATALIST_RECORDS_NOT_OBJECTS")
    return obj["total"],records

def summarize_records(pages,total):
    records=[r for page in pages for r in page]
    seen=set()
    public_ids=set()
    invalid_public=0
    private_n=0
    triplet={name:0 for name in TRIPLET}
    public_triplet={name:0 for name in TRIPLET}
    for item in records:
        no=item.get("datNo")
        security=item.get("security")
        name=item.get("name")
        safe=str(no or "")
        if safe and safe in seen:
            raise SourceError("DUPLICATED_NODE_DATA_ACCESSION")
        if safe:seen.add(safe)
        allowed=isinstance(security,str) and security.casefold()=="public"
        if not allowed:
            private_n+=1
        elif not OED.fullmatch(safe):
            invalid_public+=1
        else:
            public_ids.add(safe.upper())
        if isinstance(name,str) and name in triplet:
            triplet[name]+=1
            if allowed and OED.fullmatch(safe):
                public_triplet[name]+=1
    public_ids=sorted(public_ids)
    return {"official_file_records_seen_n":len(records),
            "official_file_records_reported_total_n":total,
            "complete_file_listing":len(records)==total,
            "public_OED_ids_verified_n":len(public_ids),
            "public_OED_ids_first_25":public_ids[:MAX_PUBLIC_REPORT],
            "public_OED_ids_output_truncated":len(public_ids)>MAX_PUBLIC_REPORT,
            "invalid_public_OED_id_records_n":invalid_public,
            "nonpublic_or_unknown_security_records_n":private_n,
            "exact_triplet_filename_count":triplet,
            "public_exact_triplet_filename_count":public_triplet,
            "patient_or_nontriplet_filenames_exported":False,
            "data_files_downloaded":False}

def inspect_node(opener=urllib.request.urlopen):
    out={"status":"NODE_SOURCE_HOLD","project_accession":PROJECT,
         "source_authorized_project_identity":False,
         "source_authorized_public_OED_list":False,
         "source_files_downloaded":False,
         "PDC_patient_biopsy_join_verified":False}
    try:
        result,meta=source_request(SEARCH,"POST",search_payload(),opener)
        out["search_http"]=meta
        hits=search_candidates(result)
        candidates=hits.pop("candidates")
        out["search"]=hits
        valid=[]
        out["general_info_schema_audits"]=[]
        def check_internal_id(internal,kind):
            entry={"candidate_origin":kind,"status":"UNVERIFIED"}
            try:
                general,gmeta=source_request(project_general_url(internal),opener=opener)
                d=data_envelope(general)
                if not isinstance(d,dict):
                    entry["status"]="GENERAL_INFO_NOT_OBJECT"
                else:
                    # Schema keys are not patient IDs and help debug source changes.
                    entry["top_level_field_names"]=sorted(str(x) for x in d.keys())[:35]
                    exact=str(d.get("projectNo","")).upper()==PROJECT
                    entry["projectNo_matches_OEP000321"]=exact
                    entry["status"]="EXACT_PROJECT_MATCH" if exact else "PROJECT_ACCESSION_MISMATCH"
                    if exact:
                        valid.append({"internal_id":internal,
                                      "response_sha256":gmeta["response_sha256"]})
            except urllib.error.HTTPError as exc:
                entry["status"]="GENERAL_INFO_HTTP_ERROR"
                entry["http_status"]=exc.code
            except Exception as exc:
                entry["status"]="GENERAL_INFO_ACCESS_HOLD"
                entry["error_class"]=type(exc).__name__
                if isinstance(exc,SourceError):entry["error_code"]=str(exc)
            out["general_info_schema_audits"].append(entry)
        for internal in candidates:
            check_internal_id(internal,"official_search_hit")
        # The accession is publicly documented, never an invented internal ID.
        # Probe that exact accession only if the search-derived IDs did not resolve.
        if not valid and PROJECT not in candidates:
            check_internal_id(PROJECT,"published_project_accession")
        out["exact_general_info_matches_n"]=len(valid)
        if len(valid)!=1:
            out["status"]="NODE_PROJECT_IDENTITY_UNCONFIRMED"
            return out
        internal=valid[0]["internal_id"]
        out["source_authorized_project_identity"]=True
        out["general_info_sha256"]=valid[0]["response_sha256"]
        total=None
        records=[]
        pages=[]
        for page in range(1,PAGE_LIMIT+1):
            payload=datalist_payload(internal,page)
            response,m=source_request(DATA_LIST,"POST",payload,opener)
            ptotal,chunk=parse_file_page(response)
            if total is not None and total!=ptotal:
                raise SourceError("NODE_DATALIST_TOTAL_CHANGED")
            total=ptotal
            records.append(chunk)
            pages.append({"page":page,"response_sha256":m["response_sha256"],
                          "records_n":len(chunk)})
            if sum(map(len,records))>=total or not chunk:
                break
        summary=summarize_records(records,total)
        out["official_datalist"]=summary
        out["page_metadata"]=pages
        out["source_authorized_public_OED_list"]=summary["public_OED_ids_verified_n"]>0
        out["status"]=("NODE_SOURCE_AUTHORIZED_DATALIST_COMPLETE" if
            summary["complete_file_listing"] else "NODE_DATALIST_PARTIAL_HOLD")
    except urllib.error.HTTPError as exc:
        out["status"]="NODE_API_HTTP_HOLD"
        out["http_status"]=exc.code
    except Exception as exc:
        out["status"]="NODE_API_SOURCE_HOLD"
        out["error_class"]=type(exc).__name__
        out["error_code"]=str(exc) if isinstance(exc,SourceError) else "UNVERIFIED_SOURCE_ERROR"
    return out

def niaid_query_url():
    return NIAID_QUERY+"?"+urllib.parse.urlencode({"q":PROJECT,"size":20,"fields":"all"})

def niaid_summary(query):
    if not isinstance(query,dict) or not isinstance(query.get("hits"),list):
        raise SourceError("NIAID_QUERY_SCHEMA_HOLD")
    hits=query["hits"][:50]
    exact=[]
    for record in hits:
        if not isinstance(record,dict):continue
        ids=record.get("identifier",[])
        values=ids if isinstance(ids,list) else [ids]
        ident=any(str(v).upper()==PROJECT for v in values)
        source_id=str(record.get("_id","")).lower()=="node_oep000321"
        if ident and source_id:exact.append(record)
    observed=set()
    for x in exact:
        downloads=x.get("distribution",[])
        if isinstance(downloads,dict):downloads=[downloads]
        if not isinstance(downloads,list):continue
        for f in downloads[:1000]:
            if not isinstance(f,dict):continue
            n=f.get("identifier","")
            if isinstance(n,str) and OED.fullmatch(n):
                observed.add(n.upper())
    return {"indexed_first_page_records_n":len(hits),
            "exact_NODE_dataset_identity_records_n":len(exact),
            "indexed_public_OED_like_identifiers_n":len(observed),
            "indexed_first_25_OED_identifiers":sorted(observed)[:MAX_PUBLIC_REPORT],
            "source_confirmed_project_membership":False,
            "index_is_not_original_file_download":True}

def inspect_niaid(opener=urllib.request.urlopen):
    out={"status":"NIAID_INDEX_HOLD","project_accession":PROJECT,
         "source_authorized_public_OED_list":False}
    try:
        obj,meta=source_request(niaid_query_url(),opener=opener)
        out["http"]=meta
        out["index"]=niaid_summary(obj)
        out["status"]="NIAID_INDEX_SEARCH_COMPLETE"
    except urllib.error.HTTPError as exc:
        out["status"]="NIAID_INDEX_HTTP_HOLD"
        out["http_status"]=exc.code
    except Exception as exc:
        out["status"]="NIAID_INDEX_SOURCE_HOLD"
        out["error_class"]=type(exc).__name__
    return out

def review(source):
    if source.get("source_files_downloaded") is not False or (
        source.get("PDC_patient_biopsy_join_verified") is not False):
        raise SourceError("INVALID_EVIDENCE_PROMOTION")
    authorized=source.get("source_authorized_project_identity",False)
    data=source.get("official_datalist",{})
    return {"schema":"hcc-stage54-node-source-postflight-v1",
            "official_project_identity_verified":authorized,
            "public_OED_ids_verified_n":data.get("public_OED_ids_verified_n",0),
            "source_RNA_expression_file_download_verified":False,
            "G30_B":"HOLD_PDC_QUANT_TO_OFFICIAL_ALIQUOT",
            "G30_C":"HOLD_NODE_RNA_MATRIX_OR_SAMPLE_LINK",
            "G30_D":"STOP_NO_INDEPENDENT_SAME_BIOPSY_RNA_PROTEIN",
            "biological_claims_allowed":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--mode",choices=("node","niaid"),required=True)
    p.add_argument("--out",required=True)
    args=p.parse_args()
    out={"schema":"hcc-stage54-authorized-source-file-metadata-v1",
         "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
         "mode":args.mode}
    if args.mode=="node":
        source=inspect_node()
        out["source"]=source
        out["review"]=review(source)
    else:
        out["source"]=inspect_niaid()
        out["review"]={"G30_B":"HOLD","G30_C":"HOLD","G30_D":"STOP",
                       "source_authoritative_node_access":False}
    path=pathlib.Path(args.out)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print("HCC_STAGE54_SOURCE",json.dumps(out))

if __name__=="__main__":main()
