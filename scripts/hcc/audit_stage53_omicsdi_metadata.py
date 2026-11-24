#!/usr/bin/env python3
"""Stage53: documented OmicsDI accession and file-metadata discovery.

Use only the published /ws/dataset/search, /ws/dataset/{source}/{id} and
/ws/dataset/{source}/{id}/files interfaces, with accession IDs taken from an
actual exact-match search response. No guessed OED IDs, no private NODE APIs,
no arbitrary download links followed, and no raw clinical or patient data.
Outputs aggregates, source hashes and fail-closed biological review gates.
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

BASE="https://www.omicsdi.org"
TARGETS={"OEP000321","HCCDB25"}
MAX_RESPONSE=1024*1024
MAX_ENTRIES=50
MAX_FILE_ITEMS=500
OED=re.compile(r"^OED[0-9]{5,}$",re.I)
ID=re.compile(r"^[A-Za-z0-9_.-]{1,80}$")
SOURCE=re.compile(r"^[A-Za-z][A-Za-z0-9_-]{0,40}$")

class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self,request,fp,code,msg,headers,newurl):
        return None

def documented_search_url(query):
    if query not in TARGETS:
        raise ValueError("unapproved exact-accession lookup")
    return BASE+"/ws/dataset/search?"+urllib.parse.urlencode(
        {"query":query,"size":20,"start":0})

def exact_dataset_url(source,accession,with_files=False):
    if not SOURCE.fullmatch(str(source)) or not ID.fullmatch(str(accession)):
        raise ValueError("untrusted OmicsDI source or accession")
    url=BASE+"/ws/dataset/"+urllib.parse.quote(str(source),safe="")+"/"+urllib.parse.quote(str(accession),safe="")
    return url+("/files" if with_files else "")

def bounded_json(url,opener=None):
    if not url.startswith(BASE+"/ws/dataset/"):
        raise ValueError("not an official documented OmicsDI URL")
    parsed=urllib.parse.urlsplit(url)
    if parsed.hostname!="www.omicsdi.org" or parsed.scheme!="https":
        raise ValueError("untrusted OmicsDI host")
    handler=opener or urllib.request.build_opener(NoRedirect()).open
    request=urllib.request.Request(url,headers={
        "User-Agent":"HCC-Stage53-public-OmicsDI-metadata/1.0",
        "Accept":"application/json",
        "Range":f"bytes=0-{MAX_RESPONSE-1}"})
    with handler(request,timeout=18) as resp:
        if resp.geturl()!=url or resp.status not in (200,206):
            raise ValueError("unexpected OmicsDI redirect/status")
        mime=resp.headers.get("Content-Type","").split(";")[0].strip().lower()
        if mime not in ("application/json","application/hal+json","text/json"):
            raise ValueError("unexpected OmicsDI content type")
        data=resp.read(MAX_RESPONSE+1)
    if len(data)>MAX_RESPONSE:
        raise ValueError("OmicsDI metadata exceeds bounded limit")
    return json.loads(data),{
        "http_status":resp.status,
        "sha256_payload":hashlib.sha256(data).hexdigest(),
        "json_bytes":len(data),
    }

def dataset_items(payload):
    if not isinstance(payload,dict):
        raise ValueError("search response is not an object")
    entries=payload.get("datasets")
    if entries is None:
        entries=payload.get("dataset",[])
    if not isinstance(entries,list) or len(entries)>MAX_ENTRIES:
        raise ValueError("unexpected search result cardinality")
    if not all(isinstance(x,dict) for x in entries):
        raise ValueError("nonobject dataset item")
    raw_count=payload.get("count")
    count=raw_count if type(raw_count) is int and raw_count>=0 else None
    return entries,count

def identifier_match(entry,query):
    """Only source-defined primary dataset accessions grant a detail lookup."""
    return any(str(entry.get(k,"")).upper()==query for k in ("id","accession"))

def search_hit_summary(payload,query):
    entries,count=dataset_items(payload)
    exact=[x for x in entries if identifier_match(x,query)]
    mention=[x for x in entries if not identifier_match(x,query)
             and query.lower() in " ".join(str(x.get(k,"")) for k in
                                          ("title","description","keywords")).lower()]
    safe_ids=[]
    for x in exact:
        source=x.get("source")
        key=x.get("id") or x.get("accession")
        if isinstance(source,str) and isinstance(key,str) and SOURCE.fullmatch(source) and ID.fullmatch(key):
            if key.upper()==query:
                safe_ids.append((source,key))
    return {"listed_count_first_page_n":len(entries),
            "total_index_search_matches_n":count,
            "exact_primary_accession_hit_n":len(exact),
            "text_mentions_without_exact_accession_n":len(mention),
            "exact_records_with_valid_source_metadata_n":len(safe_ids),
            "exact_source_identifiers":list(dict.fromkeys(safe_ids))[:3]}

def flatten_file_items(payload):
    if isinstance(payload,list):
        records=payload
    elif isinstance(payload,dict):
        records=next((payload[k] for k in ("files","dataFiles","file") if isinstance(payload.get(k),list)),[])
    else:
        records=[]
    if len(records)>MAX_FILE_ITEMS or not all(isinstance(x,dict) for x in records):
        raise ValueError("unexpected OmicsDI file record structure")
    return records

def file_metadata_summary(payload):
    records=flatten_file_items(payload)
    accessions=[]
    triplet={"HCCDB25_mRNA_level3.txt":0,
             "HCCDB25.sample.txt":0,
             "HCCDB25.patient.txt":0}
    for record in records:
        for k in ("id","accession","fileId","file_id"):
            value=record.get(k)
            if isinstance(value,str) and OED.fullmatch(value):
                accessions.append(value.upper())
        for k in ("name","filename","fileName"):
            value=record.get(k)
            if isinstance(value,str) and value in triplet:
                triplet[value]+=1
    return {"file_records_n":len(records),
            "OED_accession_fields_n":len(accessions),
            "unique_OED_accession_fields_n":len(set(accessions)),
            "HCCDB25_exact_file_names":triplet,
            "file_names_raw_ids_and_urls_exported":False,
            "NODE_project_membership_source_certified":False,
            "expression_file_bytes_downloaded":False}

def inspect_query(query,opener=None):
    out={"queried_accession":query,"status":"OMICSDI_SOURCE_METADATA_HOLD",
         "NODE_official_OED_file_id_verified":False,
         "actual_HCCDB25_processed_RNA_file_downloaded":False,
         "official_PDC_RNA_case_join_verified":False}
    try:
        data,meta=bounded_json(documented_search_url(query),opener)
        scan=search_hit_summary(data,query)
        safe_sources=scan.pop("exact_source_identifiers")
        out.update({"search_transport":meta,"search_response":scan})
        out["status"]="INDEX_QUERY_COMPLETE_NO_EXACT_MATCH"
        if not safe_sources:
            return out
        out["status"]="INDEX_EXACT_ACCESSION_MATCHED"
        out["source_detail_checks"]=[]
        for source,identifier in safe_sources:
            detail={"source_label":source.lower(),"status":"DETAIL_SOURCE_HOLD"}
            try:
                d,dm=bounded_json(exact_dataset_url(source,identifier),opener)
                detail["detail_sha256"]=dm["sha256_payload"]
                ids={str(d.get(k,"")).upper() for k in ("id","accession")}
                detail["detail_id_matches_target"]=query in ids
                if query not in ids:
                    detail["status"]="DETAIL_ID_MISMATCH_HOLD"
                else:
                    detail["status"]="OFFICIAL_OMDI_INDEX_DETAILS_FOUND"
                    try:
                        files,fm=bounded_json(exact_dataset_url(source,identifier,True),opener)
                        detail["files_response_sha256"]=fm["sha256_payload"]
                        detail["file_metadata"]=file_metadata_summary(files)
                        detail["status"]="OMDI_FILE_LIST_METADATA_READ"
                    except Exception as exc:
                        detail["file_error_class"]=type(exc).__name__
            except Exception as exc:
                detail["error_class"]=type(exc).__name__
            out["source_detail_checks"].append(detail)
    except urllib.error.HTTPError as exc:
        out["search_http_status"]=exc.code
        out["search_error_class"]="HTTPError"
    except Exception as exc:
        out["search_error_class"]=type(exc).__name__
    return out

def review(out):
    if out.get("NODE_official_OED_file_id_verified") is not False or (
        out.get("actual_HCCDB25_processed_RNA_file_downloaded") is not False):
        raise ValueError("unverified original source")
    return {"schema":"hcc-stage53-omicsdi-metadata-postflight-v1",
            "queried_accession":out["queried_accession"],
            "index_search_status":out["status"],
            "official_NODE_project_OED_data_membership_verified":False,
            "original_RNA_triplet_verified":False,
            "G30_B":"HOLD_PDC_OFFICIAL_ALIQUOT_MAPPING",
            "G30_C":"HOLD_NODE_HCCDB25_RNA_SOURCE_IDS_AND_FILE_BYTES",
            "G30_D":"STOP_NO_SAME_BIOPSY_RNA_PROTEIN_JOIN",
            "clinical_or_biological_claims_allowed":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--query",required=True,choices=sorted(TARGETS))
    p.add_argument("--out",required=True)
    args=p.parse_args()
    item={"schema":"hcc-stage53-omicsdi-documented-source-inventory-v1",
          "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
          "source":inspect_query(args.query)}
    item["review"]=review(item["source"])
    target=pathlib.Path(args.out)
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(item,indent=2)+"\n",encoding="utf-8")
    print("HCC_STAGE53_OMDI",json.dumps(item))

if __name__=="__main__":main()
