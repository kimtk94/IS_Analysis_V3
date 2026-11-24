#!/usr/bin/env python3
"""Bounded read-only PDC GraphQL quantDataMatrix access shape audit.

The direct API may expose gene-by-aliquot log2 ratios without downloading
mass spectrometry. Never retain or emit matrix cells, names, case identifiers,
aliquot IDs, raw HTTP bodies, signed URLs or per-gene measurements. Report
ONLY source response status/schema/counts. Limit response bytes to 2MiB and
network timeout. Full matrix analysis remains HOLD by construction.
"""
import argparse,datetime as dt,hashlib,json,pathlib,urllib.error,urllib.request,re,math
from probe_stage30_gao_catalog import PDC

ENDPOINT="https://pdc.cancer.gov/graphql"
MAX_BYTES=2*1024*1024

def classify_nested_array_prefix(text, start, limit=12):
    """Bounded anonymous shape of first completed values inside nested arrays.

    Never export source strings/keys or numerical values. No raw sample or
    protein identifier may leave this function.
    """
    decoder=json.JSONDecoder()
    pos=start
    depth=0
    # Enter exactly ONE wrapper array: each completed inner list/object is
    # an independent record whose shape can be checked on a bounded prefix.
    # Descending repeatedly would strip the first record's own array wrapper
    # and misclassify its individual fields as separate records.
    while pos<len(text) and text[pos] in " \t\r\n,":
        pos+=1
    if pos<len(text) and text[pos]=="[":
        pos+=1
        depth=1
    counts={"NUMBER":0,"STRING":0,"OBJECT":0,"LIST":0,"BOOLEAN":0,"NULL":0}
    numeric_field_n=0
    anonymous_field_n=0
    object_field_n=[]
    list_item_n=[]
    completed=0
    for _ in range(limit):
        while pos<len(text) and text[pos] in " \t\r\n,":
            pos+=1
        if pos>=len(text) or text[pos] in "]}":break
        try:
            obj,end=decoder.raw_decode(text,pos)
        except json.JSONDecodeError:
            break
        completed+=1
        if isinstance(obj,bool):
            counts["BOOLEAN"]+=1
        elif obj is None:
            counts["NULL"]+=1
        elif isinstance(obj,(int,float)):
            counts["NUMBER"]+=1
            if math.isfinite(obj):numeric_field_n+=1
        elif isinstance(obj,str):
            counts["STRING"]+=1
        elif isinstance(obj,dict):
            counts["OBJECT"]+=1
            object_field_n.append(len(obj))
            for value in obj.values():
                anonymous_field_n+=1
                if isinstance(value,bool) or value is None:continue
                try:
                    if math.isfinite(float(value)):numeric_field_n+=1
                except (ValueError,TypeError,OverflowError):pass
        elif isinstance(obj,list):
            counts["LIST"]+=1
            list_item_n.append(len(obj))
            for value in obj[:1000]:
                anonymous_field_n+=1
                if isinstance(value,bool) or value is None:continue
                try:
                    if math.isfinite(float(value)):numeric_field_n+=1
                except (ValueError,TypeError,OverflowError):pass
        pos=end
    return {
        "nested_array_depth_from_first_element_n":depth,
        "completed_first_nested_values_n":completed,
        "completed_nested_value_type_counts":counts,
        "max_first_object_field_count_n":max(object_field_n) if object_field_n else None,
        "max_first_list_item_count_n":max(list_item_n) if list_item_n else None,
        "finite_numeric_values_examined_n":numeric_field_n,
        "anonymized_nested_fields_examined_n":anonymous_field_n,
        "gene_sample_identifiers_exported":False,
    }

def inspect_first_complete_rows(text, start, max_rows=5):
    """Read at most 5 whole first-level rows from capped JSON array prefix.

    An early string-valued first row can be a header. No strings (including
    gene/aliquot identifiers), individual values, or any raw row leave memory.
    Subsequent rows' numeric cell fractions are evidence of tabular data,
    NOT proof of correct case mapping or expression units.
    """
    decoder=json.JSONDecoder()
    pos=start
    rows=[]
    for idx in range(max_rows):
        while pos<len(text) and text[pos] in " \t\r\n,":
            pos+=1
        if pos>=len(text) or text[pos] in "]}":break
        try:
            item,end=decoder.raw_decode(text,pos)
        except json.JSONDecodeError:
            break
        if not isinstance(item,list):break
        n=len(item)
        numeric=0;string_n=0;missing_n=0
        for v in item:
            if v is None or v=="":missing_n+=1;continue
            if isinstance(v,bool):continue
            if isinstance(v,str):string_n+=1
            try:
                if math.isfinite(float(v)):numeric+=1
            except (TypeError,ValueError,OverflowError):pass
        rows.append({
            "column_n":n,
            "string_cells_n":string_n,
            "finite_numeric_cells_n":numeric,
            "missing_cells_n":missing_n,
            "numeric_cell_fraction":round(numeric/n,6) if n else 0})
        pos=end
    first_n=rows[0]["column_n"] if rows else None
    return {
      "complete_outer_array_rows_seen_n":len(rows),
      "first_row_column_count_n":first_n,
      "all_inspected_rows_same_width":all(r["column_n"]==first_n for r in rows) if rows else False,
      "first_row_is_non_numeric_header_candidate":bool(
          rows and rows[0]["string_cells_n"]>=max(1,first_n-1) and
          rows[0]["finite_numeric_cells_n"]==0),
      "first_data_rows_inspected_n":max(0,len(rows)-1),
      "first_data_rows_finite_numeric_cells_n":sum(r["finite_numeric_cells_n"] for r in rows[1:]),
      "first_data_rows_total_cells_n":sum(r["column_n"] for r in rows[1:]),
      "first_data_rows_numeric_fraction":round(
          sum(r["finite_numeric_cells_n"] for r in rows[1:])/
          sum(r["column_n"] for r in rows[1:]),6)
          if len(rows)>1 and sum(r["column_n"] for r in rows[1:]) else None,
      "per_row_anonymous_counts":rows,
      "per_cell_names_or_values_exported":False,
    }

def inspect_oversized_matrix_prefix(raw):
    """Parse only first <=3 complete per-gene JSON objects from 2MiB prefix.

    Never return field keys, expression values, gene names, aliquot labels,
    or original text. A sample numeric count does NOT establish full matrix
    units, sample identity, gene annotations or longitudinal reproducibility.
    """
    text=raw.decode("utf-8","replace")
    match=re.search(r'"quantDataMatrix"\s*:\s*\[',text[:100000])
    if not match:
        return {"array_prefix_marker_observed":False,
                "first_complete_gene_records_seen_n":0,
                "numeric_cells_in_complete_first_records_n":0}
    pos=match.end()
    while pos<len(text) and text[pos] in " \\t\\r\\n,":
        pos+=1
    first = text[pos:pos+1]
    kinds = {'{':'OBJECT', '"':'JSON_STRING', '[':'NESTED_ARRAY',
             'n':'NULL', 't':'TRUE', 'f':'FALSE'}
    payload_kind = kinds.get(first,'OTHER_OR_INCOMPLETE')
    # Only anonymous structural counts, never property names, patient
    # identifiers, signed URLs, annotations or expression magnitudes.
    structural = {
        "first_array_element_JSON_kind":payload_kind,
        "first_two_MiB_object_open_marker_n":text.count('{'),
        "first_two_MiB_object_close_marker_n":text.count('}'),
        "first_two_MiB_array_open_marker_n":text.count('['),
        "first_two_MiB_array_close_marker_n":text.count(']'),
        "first_two_MiB_json_keyvalue_colon_marker_n":text.count(':'),
        "first_two_MiB_newline_n":text.count('\\n'),
        "first_two_MiB_literal_tab_n":text.count('\\t'),
    }
    decoder=json.JSONDecoder()
    shapes=[]
    numeric=0;values_total=0
    gene_header_seen=False
    for _ in range(3):
        while pos<len(text) and text[pos] in " \t\r\n,":
            pos+=1
        if pos>=len(text) or text[pos]!="{" :break
        try:
            record,consumed=decoder.raw_decode(text[pos:])
        except json.JSONDecodeError:
            break
        if not isinstance(record,dict):break
        shapes.append(len(record))
        gene_header_seen|=("Gene/Aliquot" in record)
        for key,value in record.items():
            if key=="Gene/Aliquot":continue
            values_total+=1
            if value in (None,""):continue
            try:
                number=float(value)
                if math.isfinite(number):numeric+=1
            except (ValueError,TypeError,OverflowError):pass
        pos+=consumed
    nested=classify_nested_array_prefix(text,pos) if payload_kind=="NESTED_ARRAY" else {}
    outer_rows=inspect_first_complete_rows(text,pos) if payload_kind=="NESTED_ARRAY" else {}
    return {"array_prefix_marker_observed":True,
            **structural,
            "nested_matrix_prefix_shapes":nested,
            "first_complete_matrix_rows_summary":outer_rows,
            "first_complete_gene_records_seen_n":len(shapes),
            "sample_gene_field_present_in_partial_rows":gene_header_seen,
            "first_record_field_count_n":shapes[0] if shapes else None,
            "min_first_three_record_field_count_n":min(shapes) if shapes else None,
            "max_first_three_record_field_count_n":max(shapes) if shapes else None,
            "numeric_cells_in_complete_first_records_n":numeric,
            "non_gene_cells_inspected_in_first_records_n":values_total,
            "numeric_cell_fraction_in_first_records":round(numeric/values_total,6) if values_total else None}

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
                        "source_response_prefix_sha256":hashlib.sha256(buffer[:MAX_BYTES]).hexdigest(),
                        "first_gene_record_numeric_structure":inspect_oversized_matrix_prefix(buffer[:MAX_BYTES])})
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
       "studies_with_structural_gene_and_numeric_value_prefix_n":sum(
          x.get("first_gene_record_numeric_structure",{}).get("first_complete_gene_records_seen_n",0)>0 and
          (x.get("first_gene_record_numeric_structure",{}).get("numeric_cells_in_complete_first_records_n",0)>0 or
           x.get("first_gene_record_numeric_structure",{}).get("first_complete_matrix_rows_summary",{}).get("first_data_rows_finite_numeric_cells_n",0)>0)
          for x in v["studies"]),
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
       "download_prefix_bytes_capped","first_gene_record_numeric_structure")} for x in res["studies"]]))
if __name__=="__main__":main()
