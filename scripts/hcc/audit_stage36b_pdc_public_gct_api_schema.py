#!/usr/bin/env python3
"""PDC public GraphQL schema-only inspection for authentic GCT/heatmap API route.

Does not make up REST endpoints. Uses introspection of exposed Query field names
only. No study sample identifiers, GCT download URLs or quantification values.
"""
import argparse
import datetime as dt
import json
import pathlib
import re
from probe_stage30_gao_catalog import pdc_query

SCHEMA_QUERY='{__schema{queryType{name}}}'
FIELDS_QUERY_TEMPLATE='{__type(name:%s){name fields{name args{name type{name kind ofType{name kind}}}}}}'
PATTERN=re.compile(r"gct|heat.?map|quant|morpheus|export",re.I)

def introspect(fetch=pdc_query):
    out={"status":"ACCESS_HOLD","source_schema_is_original_study_GCT_file":False,
         "GCT_330_biological_mapping_verified":False,"patient_ids_or_data_exported":False}
    try:
        schema,sha=fetch(SCHEMA_QUERY)
        qt=((schema.get("__schema") or {}).get("queryType") or {}).get("name")
        if not isinstance(qt,str) or not re.fullmatch(r"[_A-Za-z][_A-Za-z0-9]{0,127}",qt):
            raise ValueError("public GraphQL Query schema unavailable")
        query=FIELDS_QUERY_TEMPLATE%json.dumps(qt)
        obj,digest=fetch(query)
        fields=(obj.get("__type") or {}).get("fields")
        if not isinstance(fields,list) or len(fields)>3000 or not fields:
            raise ValueError("GraphQL Query type fields unavailable")
        matched=[]
        for f in fields:
            n=f.get("name","")
            if PATTERN.search(n):
                args=f.get("args") or []
                matched.append({"public_field_name":n,
                                "public_argument_names":[a["name"] for a in args]})
        out.update({"status":"PUBLIC_GRAPHQL_SCHEMA_INTROSPECTED",
                    "source_Query_type_fields_n":len(fields),
                    "potential_GCT_or_quant_routes_n":len(matched),
                    "public_api_relevant_fields":matched,
                    "GCT_specific_source_fields_n":sum(bool(re.search(r"gct|heat.?map|morpheus",d["public_field_name"],re.I)) for d in matched),
                    "schema_response_sha256":sha,"query_fields_response_sha256":digest})
    except Exception as e:
        out["failure_class"]=type(e).__name__
    return out

def verify(o):
    if o.get("schema")!="hcc-stage36b-PDC-public-quant-GCT-GraphQL-route-v1":
        raise ValueError("bad schema")
    ev=o.get("evidence")
    if not isinstance(ev,dict) or ev.get("status") not in ("ACCESS_HOLD","PUBLIC_GRAPHQL_SCHEMA_INTROSPECTED"):
        raise ValueError("unexpected evidence")
    for field in ("source_schema_is_original_study_GCT_file","GCT_330_biological_mapping_verified","patient_ids_or_data_exported"):
        if ev.get(field) is not False:
            raise ValueError("unverified sample proof or source value export")
    return {"schema":"hcc-stage36b-PDC-GraphQL-source-route-postflight-v1",
        "public_quant_source_schema_checked":ev["status"]=="PUBLIC_GRAPHQL_SCHEMA_INTROSPECTED",
        "GCT_download_verified":False,"PDC_quant_matrix_330_unique_aliquot_join_verified":False,
        "G30_B":"HOLD_GCT_SOURCE_NOT_DOWNLOADED_OR_MATCHED",
        "G30_C":"HOLD_NODE_RNA",
        "G30_D":"STOP_NO_MATCHED_EXTERNAL_RNA_PROTEIN",
        "biology_inference_allowed":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True);p.add_argument("--review-out",required=True)
    a=p.parse_args()
    data={"schema":"hcc-stage36b-PDC-public-quant-GCT-GraphQL-route-v1",
          "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
          "evidence":introspect()}
    decision=verify(data)
    for fn,x in ((a.out,data),(a.review_out,decision)):
        target=pathlib.Path(fn);target.parent.mkdir(parents=True,exist_ok=True)
        target.write_text(json.dumps(x,indent=2)+"\n",encoding="utf-8")
    print("HCC_STAGE36B_GRAPHQL_ROUTE",json.dumps(data))
    print("HCC_STAGE36B_POSTFLIGHT",json.dumps(decision))

if __name__=="__main__":main()
