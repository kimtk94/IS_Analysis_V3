#!/usr/bin/env python3
"""Bounded full public PDC000198/199 direct quantDataMatrix numerical QC.

Live GitHub CI only, never persist unredacted molecular data, raw JSON,
per-feature or per-patient labels, case IDs, signed URLs, or per-sample values.
The artifact is aggregate table shape/completeness, SHA256 of public source,
and safe eligibility gates. No RNA OEP000321, no case-level join, no biology.
"""
import argparse,datetime as dt,hashlib,json,math,pathlib,statistics,urllib.request
from probe_stage30_gao_catalog import PDC
from probe_stage30_gao_direct_quant_endpoint import ENDPOINT
LIMIT=48*1024*1024

def qtile(x,p):
    if not x:return None
    s=sorted(x)
    k=(len(s)-1)*p
    low=int(k);high=min(len(s)-1,low+1)
    return round(s[low]+(k-low)*(s[high]-s[low]),7)

def audit_matrix(obj):
    if not isinstance(obj,list) or len(obj)<3 or len(obj)>200000:
        raise ValueError("unusual PDC 2D quantitation array")
    header=obj[0]
    if not isinstance(header,list) or len(header)!=331:
        raise ValueError("PDC Gao quant header source columns not exactly 331")
    if any(not isinstance(v,str) for v in header):
        raise ValueError("source first header contains non-string")
    samples=header[1:]
    if len(samples)!=330 or any(not x.strip() for x in samples) or len(set(samples))!=330:
        raise ValueError("source private 330 quant columns duplicated/empty")
    n=330
    per_sample=[0]*n
    per_row=[]
    valid=0;missing=0
    zero=positive=negative=0
    gene_empty=0
    complete=0
    # avoid exposing original gene/protein label; only count distinct row label
    seen=set();duplicates=0
    for row in obj[1:]:
        if not isinstance(row,list) or len(row)!=331:raise ValueError("source row width not 331")
        label=row[0]
        if not isinstance(label,str) or not label.strip():gene_empty+=1
        else:
            if label in seen:duplicates+=1
            seen.add(label)
        number=0
        for j,v in enumerate(row[1:]):
            if v is None or (isinstance(v,str) and not v.strip()):
                missing+=1;continue
            try:
                if isinstance(v,bool):raise ValueError("boolean not a quantitative ratio")
                z=float(v)
                if not math.isfinite(z):raise ValueError("nonfinite ratio")
                number+=1;valid+=1;per_sample[j]+=1
                if z==0:zero+=1
                elif z>0:positive+=1
                else:negative+=1
            except (TypeError,ValueError,OverflowError):
                missing+=1
        per_row.append(number)
        if number==n:complete+=1
    if not per_row:raise ValueError("source has no quantified feature rows")
    total=len(per_row)*n
    if valid+missing!=total:raise ValueError("invalid source tally")
    return {"quant_source_matrix_row_n":len(per_row),
      "quant_header_feature_label_n":1,
      "quant_source_anonymous_sample_columns_n":n,
      "all_rows_width_331":True,
      "unique_nonempty_feature_label_n":len(seen),
      "empty_feature_label_row_n":gene_empty,
      "duplicate_feature_label_row_n":duplicates,
      "total_numeric_candidate_cells_n":total,
      "finite_numeric_cell_n":valid,
      "nonfinite_missing_or_nonnumeric_cell_n":missing,
      "finite_numeric_fraction":round(valid/total,9),
      "positive_value_cell_n":positive,
      "negative_value_cell_n":negative,
      "zero_value_cell_n":zero,
      "rows_with_all_330_numeric_cells_n":complete,
      "rows_with_at_least_one_numeric_cell_n":sum(v>0 for v in per_row),
      "per_feature_numeric_count_q10":qtile(per_row,.1),
      "per_feature_numeric_count_median":qtile(per_row,.5),
      "per_feature_numeric_count_q90":qtile(per_row,.9),
      "per_sample_numeric_feature_count_q10":qtile(per_sample,.1),
      "per_sample_numeric_feature_count_median":qtile(per_sample,.5),
      "per_sample_numeric_feature_count_q90":qtile(per_sample,.9),
      "gene_or_patient_names_or_individual_values_exported":False,
      "complete_numeric_biological_sample_ID_join_verified":False}

def source(pid):
    out={"pdc_study_id":pid,"status":"ACCESS_HOLD",
         "unredacted_source_saved_to_public_artifacts":False,
         "same_patient_RNA_protein_join_verified":False}
    req=urllib.request.Request(ENDPOINT,data=json.dumps({"query":
        '{quantDataMatrix(pdc_study_id:"'+pid+'",data_type:"log2_ratio")}'}).encode(),
        headers={"Content-Type":"application/json","Accept":"application/json",
                 "User-Agent":"HCC-Stage30-Full-Quant-Numeric-QC-ReadOnly/1"},method="POST")
    try:
        with urllib.request.urlopen(req,timeout=105) as rsp:
            raw=rsp.read(LIMIT+1)
            status=rsp.status
        if status!=200:raise ValueError("PDC source HTTP status not 200")
        if len(raw)>LIMIT:
            out.update({"status":"SOURCE_EXCEEDS_48MiB_BUDGET",
                        "source_bytes_read_cap_n":LIMIT})
            return out
        digest=hashlib.sha256(raw).hexdigest()
        payload=json.loads(raw)
        if not isinstance(payload,dict) or payload.get("errors") or not isinstance(payload.get("data"),dict):
            raise ValueError("PDC source GraphQL returned error")
        matrix=payload["data"].get("quantDataMatrix")
        metrics=audit_matrix(matrix)
        out.update({"status":"FULL_PDC_DIRECT_QUANT_TABLE_NUMERIC_QC_SUCCEEDED",
                    "source_HTTP_status":status,
                    "complete_source_response_bytes_n":len(raw),
                    "original_full_source_SHA256":digest,
                    "requested_PDC_quant_data_type":"log2_ratio",
                    "no_per_patient_or_gene_values_persisted":True,
                    **metrics})
    except Exception as e:
        out["source_failure_class"]=type(e).__name__
    return out

def verify(v):
    if v.get("schema")!="hcc-stage30-PDC-full-quant-table-numeric-QC-v1" or {x.get("pdc_study_id") for x in v.get("studies",[])}!=set(PDC):
        raise ValueError("unexpected PDC source")
    for x in v["studies"]:
        if x.get("unredacted_source_saved_to_public_artifacts") is not False or x.get("same_patient_RNA_protein_join_verified") is not False:
            raise ValueError("unverified biology or source privacy")
        if x["status"] not in ("ACCESS_HOLD","SOURCE_EXCEEDS_48MiB_BUDGET",
             "FULL_PDC_DIRECT_QUANT_TABLE_NUMERIC_QC_SUCCEEDED"):
            raise ValueError("unsupported PDC source state")
        if x["status"]=="FULL_PDC_DIRECT_QUANT_TABLE_NUMERIC_QC_SUCCEEDED":
            if x["quant_source_anonymous_sample_columns_n"]!=330 or x["quant_source_matrix_row_n"]<10:
                raise ValueError("source matrix missing")
            if x["finite_numeric_cell_n"]+x["nonfinite_missing_or_nonnumeric_cell_n"]!=x["total_numeric_candidate_cells_n"]:
                raise ValueError("numeric mass conservation")
            if x["gene_or_patient_names_or_individual_values_exported"] is not False or x["complete_numeric_biological_sample_ID_join_verified"] is not False:
                raise ValueError("source identifiers leaked or unverified join")
    n=sum(x["status"]=="FULL_PDC_DIRECT_QUANT_TABLE_NUMERIC_QC_SUCCEEDED" for x in v["studies"])
    return {"schema":"hcc-stage30-PDC-full-matrix-QC-postflight-v1",
       "full_original_PDC_quantitative_matrix_numeric_QC_study_n":n,
       "G30_B":"GO_SOURCE_NUMERIC_QC_HOLD_CASE_ALIQUOT_CHANNEL_LINK" if n==2 else "HOLD_PDC_FULL_NUMERIC_QC",
       "G30_C":"HOLD_NODE_OEP000321_RNA_FILE_AND_PATIENT_ID",
       "G30_D":"STOP_NO_SAME_PERSON_INDEPENDENT_RNA_PROTEIN_JOIN",
       "external_gene_by_gene_biological_replication_done":False}

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--out",required=True);ap.add_argument("--review-out",required=True)
    opts=ap.parse_args()
    d={"schema":"hcc-stage30-PDC-full-quant-table-numeric-QC-v1",
       "checked_at_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
       "download_budget_bytes_per_source":LIMIT,
       "studies":[source(pid) for pid in PDC]}
    path=pathlib.Path(opts.out);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(d,indent=2)+"\n")
    x=verify(d);q=pathlib.Path(opts.review_out);q.parent.mkdir(parents=True,exist_ok=True)
    q.write_text(json.dumps(x,indent=2)+"\n")
    print(json.dumps(x,indent=2))
    print("PDC_FULL_NUMERIC_QC",json.dumps(d["studies"]))
if __name__=="__main__":main()
