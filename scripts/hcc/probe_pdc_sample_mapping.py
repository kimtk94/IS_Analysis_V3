#!/usr/bin/env python3
"""PDC internal quantitative-matrix header vs sample-table correspondence audit.

Only public file prefixes / tiny sample metadata are read; counts are emitted, no
patient/aliquot/sample identifiers or signed URLs are stored or displayed.
Not an RNA-vs-protein join and not a clinical outcome validation.
"""
import argparse
import csv
import datetime as dt
import io
import json
import pathlib
import re

from probe_pdc_quant_headers import PDC_FILES, get_file_link, bounded_prefix

def normalize(x):
    return re.sub(r"[^a-z0-9]", "", str(x).lower())

def candidate_value(x):
    x=normalize(x)
    return x if len(x)>=5 and any(ch.isalpha() for ch in x) else ""

def read_table(raw, row_cap=1000):
    text=raw.decode("utf-8-sig",errors="replace")
    reader=csv.reader(io.StringIO(text),delimiter="\t")
    rows=[]
    for i,row in enumerate(reader):
        if i>=row_cap:break
        rows.append(row)
    return rows

def schema_summary(matrix, metadata):
    """Compare public matrix header columns with values from one sample metadata column."""
    if not matrix or not metadata:
        return {"status":"EMPTY_TABLE"}
    headers=matrix[0]
    possible={candidate_value(x) for x in headers}
    possible.discard("")
    meta_headers=metadata[0]
    scores=[]
    for j,label in enumerate(meta_headers):
        names={candidate_value(row[j]) for row in metadata[1:] if j<len(row)}
        names.discard("")
        matched=names & possible
        category=next((word for word in ("sample","case","patient","aliquot","channel",
                  "plex","run","batch","id") if word in label.lower()),"other")
        scores.append({"metadata_column_index":j,"column_role":category,
                       "distinct_nonempty_identifier_like_values":len(names),
                       "exact_normalized_header_matches":len(matched)})
    scores.sort(key=lambda x:(-x["exact_normalized_header_matches"],x["metadata_column_index"]))
    return {
        "status":"CORRESPONDENCE_SCREENED",
        "matrix_header_columns":len(headers),
        "matrix_header_unique_normalized_identifier_like":len(possible),
        "metadata_columns":len(meta_headers),
        "metadata_data_rows":max(0,len(metadata)-1),
        "top_5_metadata_fields_by_identifier_overlap":scores[:5],
        "any_exact_normalized_overlap":any(x["exact_normalized_header_matches"]>0 for x in scores),
        "full_expression_matrix_join_verified":False,
        "rna_protein_same_patient_join_verified":False
    }

def get_prefix(study,filename):
    _,url=get_file_link(study,filename)
    if not url:raise ValueError("missing signed URL")
    payload,code,ranged=bounded_prefix(url)
    return payload

def probe():
    records={}
    for study,files in PDC_FILES.items():
        result={"study":study,"matrix_file":files[0],"sample_file":files[1],
                "status":"UNVERIFIED"}
        try:
            matrix=read_table(get_prefix(study,files[0]),row_cap=2)
            metadata=read_table(get_prefix(study,files[1]),row_cap=1000)
            result["schema_comparison"]=schema_summary(matrix,metadata)
            result["status"]="PUBLIC_MATRIX_METADATA_CORRESPONDENCE_SCREENED"
        except Exception as ex:
            result["error_type"]=type(ex).__name__
        records[study]=result
    return records

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",default="results/hcc/feasibility/pdc_sample_mapping.json")
    args=parser.parse_args()
    result={"checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
        "records":probe(),
        "scope":"No patient IDs, sample IDs, protein expression values or signed URLs written",
        "warnings":[
            "Matching one metadata label to matrix header is only a naming sanity check.",
            "Patient/biopsy RNA-to-protein join is NOT validated by this audit."
        ]}
    path=pathlib.Path(args.output)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
