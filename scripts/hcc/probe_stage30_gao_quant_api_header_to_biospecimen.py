#!/usr/bin/env python3
"""Stage30: source-authoritative PDC quantDataMatrix column ↔ biospecimen IDs.

Queries public API and reads at most 2 MiB of each log2_ratio JSON response.
Parses only the first complete nested header row in memory. Then requests
official PDC biospecimen case/sample/aliquot mapping. Exports only counts:
no gene names, patient/sample/aliquot IDs, individual per-sample hashes,
file URLs or molecular values. PDC same-study metadata alias intersection
is NOT an independent RNA OEP000321 join or full quantitative QC.
"""
import argparse,datetime as dt,hashlib,json,pathlib,re,urllib.request
from probe_stage30_gao_catalog import pdc_query,PDC
from probe_stage30_gao_direct_quant_endpoint import ENDPOINT,MAX_BYTES,inspect_first_complete_rows

def get_quant_header(pid):
    query='{quantDataMatrix(pdc_study_id:"'+pid+'",data_type:"log2_ratio")}'
    req=urllib.request.Request(ENDPOINT,data=json.dumps({"query":query}).encode(),
        headers={"Content-Type":"application/json","User-Agent":"HCC-Stage30-Quant-Header-Biopsy-QC/1"},method="POST")
    with urllib.request.urlopen(req,timeout=45) as rsp:
        raw=rsp.read(MAX_BYTES+1)
        status=rsp.status
    if status!=200:raise ValueError("PDC quant header HTTP status invalid")
    limited=raw[:MAX_BYTES]
    text=limited.decode("utf-8","replace")
    m=re.search(r'"quantDataMatrix"\s*:\s*\[',text[:120000])
    if not m:raise ValueError("missing quantDataMatrix top-level array")
    start=m.end()
    while start<len(text) and text[start] in " \t\r\n,":
        start+=1
    if start>=len(text) or text[start]!="[":
        raise ValueError("first quant element is not header list")
    decoder=json.JSONDecoder()
    header,end=decoder.raw_decode(text,start)
    if not isinstance(header,list) or not 2<len(header)<1500:
        raise ValueError("invalid or unavailable source quant header")
    if any(not isinstance(x,str) for x in header):
        raise ValueError("expected string quant column labels")
    # First cell is label for gene/protein feature, not a sample identifier.
    ids=[v.strip() for v in header[1:]]
    if len(ids)!=len(set(ids)) or any(not x for x in ids):
        raise ValueError("nonunique or empty quant names")
    anonymous=inspect_first_complete_rows(text,start,max_rows=5)
    return {"status":"VERIFIED_PDC_QUANT_JSON_HEADER",
       "http_status":status,"read_capped_at_bytes":MAX_BYTES,
       "response_was_larger_than_byte_cap":len(raw)>MAX_BYTES,
       "source_prefix_sha256":hashlib.sha256(limited).hexdigest(),
       "first_header_field_count_n":len(header),
       "candidate_column_names_n":len(ids),
       "first_header_gene_label_present":header[0].strip().lower() in ("gene/aliquot","gene","protein"),
       "bounded_postheader_rows_seen_n":anonymous["first_data_rows_inspected_n"],
       "bounded_postheader_numeric_cells_n":anonymous["first_data_rows_finite_numeric_cells_n"],
       "bounded_postheader_widths_consistent":anonymous["all_inspected_rows_same_width"],
       "per_column_or_per_gene_labels_exported":False,
       "raw_numeric_cells_exported":False}, ids

def official_biospecimen(pid):
    q='{biospecimenPerStudy(pdc_study_id:"'+pid+'"){case_id case_submitter_id sample_id sample_submitter_id aliquot_id aliquot_submitter_id sample_type}}'
    data,_=pdc_query(q)
    rows=data.get("biospecimenPerStudy")
    if not isinstance(rows,list) or not 100<=len(rows)<=10000:
        raise ValueError("PDC biospecimen count unexpected")
    bynamespace={k:{} for k in ("case_id","case_submitter_id","sample_id",
                                "sample_submitter_id","aliquot_id","aliquot_submitter_id")}
    for r in rows:
        case=str(r.get("case_id") or "").strip()
        sample=str(r.get("sample_id") or "").strip()
        aliquot=str(r.get("aliquot_id") or "").strip()
        if not case or not sample or not aliquot:
            raise ValueError("official PDC biospecimen missing UUID")
        typ=str(r.get("sample_type") or "").lower()
        label=("Normal" if "normal" in typ or "non-tumor" in typ else
          "Tumor" if "tumor" in typ or "primary" in typ or "malignant" in typ else "OtherOrUnknown")
        for k in bynamespace:
            val=r.get(k)
            if not val:continue
            key=str(val).strip()
            # Case IDs legitimately occur in multiple tumour/normal samples.
            # Retain ALL records per key; never infer that a case ID uniquely
            # determines an aliquot.
            record={"case":case,"sample":sample,"aliquot":aliquot,"type":label}
            bynamespace[k].setdefault(key,[])
            if record not in bynamespace[k][key]:
                bynamespace[k][key].append(record)
    return bynamespace

def score(header,namespace):
    rawset=set(header)
    common=rawset&set(namespace)
    matches=[m for k in common for m in namespace[k]]
    ci={x.casefold() for x in header}&{x.casefold() for x in namespace}
    return {"exact_source_header_alias_n":len(common),
        "casefold_source_header_alias_n":len(ci),
        "distinct_exact_matched_cases_n":len({m["case"] for m in matches}),
        "distinct_exact_matched_samples_n":len({m["sample"] for m in matches}),
        "distinct_exact_matched_aliquots_n":len({m["aliquot"] for m in matches}),
        "exact_matched_Tumor_alias_n":sum(any(m["type"]=="Tumor" for m in namespace[k]) for k in common),
        "exact_matched_Normal_alias_n":sum(any(m["type"]=="Normal" for m in namespace[k]) for k in common),
        "exact_matched_OtherOrUnknown_alias_n":sum(any(m["type"]=="OtherOrUnknown" for m in namespace[k]) for k in common),
        "alias_keys_with_multiple_biospecimen_records_n":sum(len(namespace[k])>1 for k in common)}

def study(pid):
    out={"pdc_study_id":pid,"status":"ACCESS_HOLD",
        "RNA_sample_identity_verified":False,
        "full_quantitative_matrix_verified":False,
        "individual_case_sample_aliquot_or_gene_IDs_exported":False}
    try:
        src,ids=get_quant_header(pid)
        nss=official_biospecimen(pid)
        measured={k:score(ids,v) for k,v in nss.items()}
        best=max(measured,key=lambda k:measured[k]["exact_source_header_alias_n"])
        out.update(src)
        out.update({"status":"PDC_DIRECT_QUANT_HEADER_VS_BIOSPECIMEN_COMPARED",
            "source_namespace_match_counts":measured,
            "best_matching_official_namespace":best,
            "best_matching_header_columns_n":measured[best]["exact_source_header_alias_n"],
            "PDC_source_aliquot_namespace_verified":best in ("aliquot_id","aliquot_submitter_id") and measured[best]["exact_source_header_alias_n"]>0,
            "PDC_biospecimen_metadata_connected_to_header_values_only":True,
            "RNA_sample_identity_verified":False,
            "full_quantitative_matrix_verified":False,
            "individual_case_sample_aliquot_or_gene_IDs_exported":False})
    except Exception as e:
        out["failure_category"]=type(e).__name__
    return out

def verify(v):
    if v.get("schema")!="hcc-stage30-direct-pdc-header-biospecimen-match-v1" or {s.get("pdc_study_id") for s in v.get("studies",[])}!=set(PDC):
        raise ValueError("source studies invalid")
    for x in v["studies"]:
        if x.get("status") not in ("ACCESS_HOLD","PDC_DIRECT_QUANT_HEADER_VS_BIOSPECIMEN_COMPARED"):
            raise ValueError("status invalid")
        for forbidden in ("RNA_sample_identity_verified","full_quantitative_matrix_verified",
                          "individual_case_sample_aliquot_or_gene_IDs_exported"):
            if x.get(forbidden) is not False:raise ValueError("patient or assay overclaim")
        if x["status"]=="PDC_DIRECT_QUANT_HEADER_VS_BIOSPECIMEN_COMPARED":
            n=x.get("candidate_column_names_n")
            if type(n)!=int or not 1<=n<=1500:raise ValueError("numeric shape wrong")
            if not isinstance(x.get("source_namespace_match_counts"),dict):
                raise ValueError("source namespace missing")
            for y in x["source_namespace_match_counts"].values():
                if not 0<=y["exact_source_header_alias_n"]<=n:
                    raise ValueError("source alias count exceeds header width")
            if x.get("per_column_or_per_gene_labels_exported") is not False or x.get("raw_numeric_cells_exported") is not False:
                raise ValueError("quant source leak")
    return {"schema":"hcc-stage30-direct-quant-biospecimen-postflight-v1",
        "direct_PDC_quant_header_and_official_biospecimen_count_verified_study_n":
          sum(x["status"]=="PDC_DIRECT_QUANT_HEADER_VS_BIOSPECIMEN_COMPARED" for x in v["studies"]),
        "PDC_direct_quant_columns_verified_as_official_aliquot_namespace_study_n":
          sum(x.get("PDC_source_aliquot_namespace_verified") is True for x in v["studies"]),
        "G30_B":"PARTIAL_PDC_DIRECT_QUANT_HEADER_ONLY_FULL_MATRIX_HOLD",
        "G30_C":"HOLD_NODE_RNA_EXPRESSION_MATRIX",
        "G30_D":"STOP_NO_INDEPENDENT_SAME_PATIENT_RNA_PROTEIN_JOIN",
        "external_biological_validation_complete":False}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True);ap.add_argument("--review-out",required=True)
    a=ap.parse_args()
    doc={"schema":"hcc-stage30-direct-pdc-header-biospecimen-match-v1",
         "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
         "studies":[study(p) for p in PDC]}
    dest=pathlib.Path(a.out);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(doc,indent=2)+"\n")
    summary=verify(doc)
    target=pathlib.Path(a.review_out);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))
    print("PDC_DIRECT_QUANT_ALIAS_AGGREGATES",json.dumps([
        {k:x.get(k) for k in ("pdc_study_id","status","failure_category",
             "first_header_field_count_n","first_header_gene_label_present",
             "candidate_column_names_n","bounded_postheader_rows_seen_n",
             "bounded_postheader_numeric_cells_n","bounded_postheader_widths_consistent",
             "best_matching_official_namespace","best_matching_header_columns_n",
             "source_namespace_match_counts","PDC_source_aliquot_namespace_verified")}
        for x in doc["studies"]]))
if __name__=="__main__":main()
