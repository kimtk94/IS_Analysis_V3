#!/usr/bin/env python3
"""Bounded read-only PDC GraphQL quantDataMatrix access shape audit.

The direct API may expose gene-by-aliquot log2 ratios without downloading
mass spectrometry. Never retain or emit matrix cells, names, case identifiers,
aliquot IDs, raw HTTP bodies, signed URLs or per-gene measurements. Report
ONLY source response status/schema/counts. Limit response bytes to 2MiB and
network timeout. Full matrix analysis remains HOLD by construction.
"""
import argparse,datetime as dt,hashlib,json,pathlib,urllib.error,urllib.request
from probe_stage30_gao_catalog import PDC

ENDPOINT="https://pdc.cancer.gov/graphql"
MAX_BYTES=2*1024*1024

def direct(pid):
    out={"pdc_study_id":pid,"requested_data_type":"log2_ratio",
         "status":"ACCESS_HOLD","full_matrix_fetched_and_verified":False,
         "sample_or_gene_names_exported":False,
         "cell_values_exported":False}
    if pid not in PDC:raise ValueError("unexpected study")
    query='{quantDataMatrix(pdc_study_id:"'+pid+'",data_type:"log2_ratio")}'
    req=urllib.request.Request(ENDPOINT,data=json.dumps({"query":query}).encode(),
       headers={"Content-Type":"application/json","Accept":"application/json",
                "User-Agent":"HCC-Stage30-QuantMatrix-BoundedSchema/1"},method="POST")
    try:
        with urllib.request.urlopen(req,timeout=35) as rsp:
            buffer=rsp.read(MAX_BYTES+1)
            out["http_status"]=rsp.status
            out["response_declared_content_length_bytes"]=(int(rsp.headers["Content-Length"])
                 if rsp.headers.get("Content-Length","").isdigit() else None)
        if len(buffer)>MAX_BYTES:
            out.update({"status":"RESPONSE_OVER_BYTE_LIMIT",
                        "download_prefix_bytes_capped":MAX_BYTES,
                        "source_response_prefix_sha256":hashlib.sha256(buffer[:MAX_BYTES]).hexdigest()})
            return out
        result=json.loads(buffer)
        out["whole_response_sha256"]=hashlib.sha256(buffer).hexdigest()
        if not isinstance(result,dict):raise ValueError("JSON root is not object")
        if result.get("errors"):
            messages=[str(x.get("message","")) for x in result["errors"] if isinstance(x,dict)]
            out["status"]="GRAPHQL_SCHEMA_OR_API_HOLD"
            out["reported_error_n"]=len(messages)
            out["error_category"]="FIELD_SELECTION_REQUIRED" if any("selection" in m.lower() for m in messages) else (
                "INVALID_ARGUMENT_OR_UNSUPPORTED_TYPE" if any("argument" in m.lower() for m in messages) else "OTHER_GRAPHQL")
            return out
        item=(result.get("data") or {}).get("quantDataMatrix")
        if item is None:
            out["status"]="MATRIX_NOT_RETURNED"
        elif isinstance(item,list):
            out["status"]="BOUNDED_QUANT_RESPONSE_STRUCTURE_VISIBLE"
            out["response_rows_under_byte_cap_n"]=len(item)
            # In memory only, do not export any actual keys:
            if item and isinstance(item[0],dict):
                out["first_row_field_count_n"]=len(item[0])
                out["rows_with_dictionary_shape_n"]=sum(isinstance(row,dict) for row in item)
        elif isinstance(item,str):
            out["status"]="STRING_PAYLOAD_RETURNED_UNINTERPRETED"
            out["string_payload_size_bytes"]=len(item.encode())
        else:
            out["status"]="UNRECOGNIZED_RESPONSE_KIND"
        return out
    except urllib.error.HTTPError as e:
        # Never log HTTP error body or endpoint query (can contain source details).
        out["status"]="HTTP_ACCESS_HOLD";out["http_status"]=e.code
    except Exception as e:
        out["status"]="ACCESS_HOLD";out["failure_category"]=type(e).__name__
    return out

def review(v):
    if v.get("schema")!="hcc-stage30-gao-direct-quant-api-bounded-v1" or len(v.get("studies",[]))!=2 or {x.get("pdc_study_id") for x in v["studies"]}!=set(PDC):
        raise ValueError("source quantitative preview contract invalid")
    for x in v["studies"]:
        if x.get("full_matrix_fetched_and_verified") is not False or x.get("sample_or_gene_names_exported") is not False or x.get("cell_values_exported") is not False:
            raise ValueError("potential biological/identifier privacy overclaim")
        if x.get("status") not in ("ACCESS_HOLD","HTTP_ACCESS_HOLD","RESPONSE_OVER_BYTE_LIMIT",
          "GRAPHQL_SCHEMA_OR_API_HOLD","MATRIX_NOT_RETURNED","BOUNDED_QUANT_RESPONSE_STRUCTURE_VISIBLE",
          "STRING_PAYLOAD_RETURNED_UNINTERPRETED","UNRECOGNIZED_RESPONSE_KIND"):
            raise ValueError("unexpected source quant status")
    return {"schema":"hcc-stage30-direct-PDC-quant-matrix-endpoint-postflight-v1",
       "studies_with_bounded_quant_API_structure_visible_n":
          sum(x["status"]=="BOUNDED_QUANT_RESPONSE_STRUCTURE_VISIBLE" for x in v["studies"]),
       "studies_with_over_byte_budget_n":
          sum(x["status"]=="RESPONSE_OVER_BYTE_LIMIT" for x in v["studies"]),
       "G30_B":"HOLD_FULL_QUANTITATIVE_GENE_ALIQUOT_MATRIX_AND_UNITS_NOT_VERIFIED",
       "G30_C":"HOLD_EXTERNAL_NODE_RNA_MATRIX",
       "G30_D":"STOP_NO_MATCHED_RNA_PROTEIN_PATIENT_COHORT",
       "independent_biological_validation_performed":False}

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True);p.add_argument("--review-out",required=True)
    opts=p.parse_args()
    res={"schema":"hcc-stage30-gao-direct-quant-api-bounded-v1",
         "checked_at_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
         "bounded_max_response_bytes":MAX_BYTES,
         "studies":[direct(pid) for pid in PDC],
         "case_and_sample_matching_performed":False}
    tgt=pathlib.Path(opts.out);tgt.parent.mkdir(parents=True,exist_ok=True)
    tgt.write_text(json.dumps(res,indent=2)+"\n")
    review_data=review(res);p2=pathlib.Path(opts.review_out);p2.parent.mkdir(parents=True,exist_ok=True)
    p2.write_text(json.dumps(review_data,indent=2)+"\n")
    print(json.dumps(review_data,indent=2))
    print("QUANT_API_BOUNDED_STATUS",json.dumps([{k:x.get(k) for k in
      ("pdc_study_id","status","http_status","failure_category","reported_error_n",
       "error_category","response_rows_under_byte_cap_n","first_row_field_count_n",
       "download_prefix_bytes_capped")} for x in res["studies"]]))
if __name__=="__main__":main()
