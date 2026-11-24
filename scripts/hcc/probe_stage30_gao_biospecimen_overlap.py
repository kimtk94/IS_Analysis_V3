#!/usr/bin/env python3
"""Read-only official PDC 198/199 biospecimen case/sample/aliquot overlap.

Patient/aliquot labels are exclusively in runner process memory, never saved,
printed, hashed individually, or published. PDC proteome-vs-phospho metadata
overlap does NOT prove numerical source abundance or RNA availability.
"""
import argparse
import collections
import datetime as dt
import json
import pathlib
from probe_stage30_gao_catalog import PDC,pdc_query

QUERY='{ biospecimenPerStudy(pdc_study_id:"%s") { case_id case_submitter_id sample_id sample_submitter_id aliquot_id aliquot_submitter_id sample_type } }'

def normalize_kind(s):
    x=str(s or "").strip().lower()
    if "normal" in x or "non-tumor" in x or "non tumor" in x:
        return "Normal"
    if "tumor" in x or "malignant" in x or "primary" in x:
        return "Tumor"
    if "blood" in x:return "Blood"
    return "OtherOrUnknown"

def accession(v):
    if v is None:return ""
    return str(v).strip()

def records(pid):
    try:
        d,sha=pdc_query(QUERY%pid)
        raw=d.get("biospecimenPerStudy")
        if not isinstance(raw,list):raise ValueError("unexpected source biospecimen records")
        if not 0<len(raw)<=10000:raise ValueError("unexpected biospecimen item count")
        cases=set();samples=set();aliquots=set()
        categories=collections.Counter();category_cases=collections.defaultdict(set)
        case_sample_types=collections.defaultdict(set)
        missing={"case":0,"sample":0,"aliquot":0}
        for r in raw:
            if not isinstance(r,dict):continue
            # Use stable submitter name if present; fallback to PDC UUID.
            ca=accession(r.get("case_submitter_id") or r.get("case_id"))
            sa=accession(r.get("sample_submitter_id") or r.get("sample_id"))
            al=accession(r.get("aliquot_submitter_id") or r.get("aliquot_id"))
            for label,value in (("case",ca),("sample",sa),("aliquot",al)):
                if not value:missing[label]+=1
            if ca:cases.add(ca)
            if sa:samples.add(sa)
            if al:aliquots.add(al)
            typ=normalize_kind(r.get("sample_type"))
            categories[typ]+=1
            if ca:
                category_cases[typ].add(ca)
                case_sample_types[ca].add(typ)
        within_pairs=sum({"Normal","Tumor"}<=flags for flags in case_sample_types.values())
        summary={"status":"VERIFIED_SOURCE_BIOSPECIMEN_METADATA",
                 "pdc_study_id":pid,"response_sha256":sha,
                 "raw_biospecimen_row_n":len(raw),
                 "unique_case_n":len(cases),
                 "unique_sample_n":len(samples),
                 "unique_aliquot_n":len(aliquots),
                 "sample_type_rows":dict(sorted(categories.items())),
                 "case_counts_by_type":{k:len(v) for k,v in sorted(category_cases.items())},
                 "cases_with_both_tumor_and_normal_metadata_n":within_pairs,
                 "missing_identifier_rows":missing,
                 "actual_molecular_matrix_column_samples_verified":False,
                 "private_case_sample_aliquot_ids_exported":False}
        return summary,{"cases":cases,"samples":samples,"aliquots":aliquots,
                        "tumor":category_cases["Tumor"],"normal":category_cases["Normal"]}
    except Exception as exc:
        return {"status":"ACCESS_HOLD","pdc_study_id":pid,
            "error_type":type(exc).__name__,
            "actual_molecular_matrix_column_samples_verified":False,
            "private_case_sample_aliquot_ids_exported":False},None

def overlap(a,b):
    if a is None or b is None:
        return {"status":"SOURCE_MATCHING_NOT_AVAILABLE",
            "PDC_same_source_protein_phospho_metadata_only":True,
            "RNA_case_link_verified":False,
            "case_or_sample_ids_exported":False}
    return {"status":"PDC_PROTEIN_PHOSPHO_BIOSPECIMEN_AGGREGATES",
        "common_submitter_case_n":len(a["cases"]&b["cases"]),
        "common_submitter_sample_n":len(a["samples"]&b["samples"]),
        "common_submitter_aliquot_n":len(a["aliquots"]&b["aliquots"]),
        "case_has_tumor_in_both_studies_n":len(a["tumor"]&b["tumor"]),
        "case_has_normal_in_both_studies_n":len(a["normal"]&b["normal"]),
        "PDC_same_source_protein_phospho_metadata_only":True,
        "RNA_case_link_verified":False,
        "quant_matrix_header_link_verified":False,
        "case_or_sample_ids_exported":False}

def verify(x):
    if x.get("schema")!="hcc-stage30-gao-proteome-phospho-biospecimen-v1":
        raise ValueError("source schema")
    rows=x.get("studies")
    if not isinstance(rows,list) or {q.get("pdc_study_id") for q in rows}!=set(PDC):
        raise ValueError("PDC target drift")
    for q in rows:
        if q.get("private_case_sample_aliquot_ids_exported") is not False or q.get("actual_molecular_matrix_column_samples_verified") is not False:
            raise ValueError("sample identity evidence overclaim/leak")
        if q.get("status")=="VERIFIED_SOURCE_BIOSPECIMEN_METADATA":
            for key in ("unique_case_n","unique_sample_n","unique_aliquot_n"):
                n=q.get(key)
                if type(n) is not int or not 0<n<=10000:raise ValueError("unexpected unique ID count")
        elif q.get("status")!="ACCESS_HOLD":
            raise ValueError("invalid metadata status")
    o=x.get("proteome_phospho_metadata_intersection")
    if (not isinstance(o,dict) or o.get("PDC_same_source_protein_phospho_metadata_only") is not True
            or o.get("RNA_case_link_verified") is not False
            or o.get("case_or_sample_ids_exported") is not False):
        raise ValueError("unsupported RNA/genomic join")
    return {"schema":"hcc-stage30-pdc-biospecimen-overlap-postflight-v1",
        "PDC_biospecimen_metadata_verified_study_n":sum(q["status"]=="VERIFIED_SOURCE_BIOSPECIMEN_METADATA" for q in rows),
        "case_count_shared_between_protein_and_phospho_source_metadata_n":o.get("common_submitter_case_n"),
        "G30_A":"METADATA_CASE_CANDIDATES_ONLY",
        "G30_B":"HOLD_FULL_PROTEIN_GENE_BY_ALIQUOT_MATRIX",
        "G30_C":"HOLD_RNA_SOURCE_MATRIX",
        "G30_D":"STOP_NO_MATCHED_INDEPENDENT_RNA_PROTEOME_COHORT",
        "same_study_case_id_overlap_does_not_imply_independent_biological_validation":True,
        "biological_replication_completed":False}

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True);p.add_argument("--review-out",required=True)
    a=p.parse_args()
    source=[records(pid) for pid in PDC]
    x={"schema":"hcc-stage30-gao-proteome-phospho-biospecimen-v1",
       "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
       "studies":[v[0] for v in source],
       "proteome_phospho_metadata_intersection":overlap(source[0][1],source[1][1])}
    out=pathlib.Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(x,indent=2)+"\n",encoding="utf-8")
    reviewed=verify(x);dest=pathlib.Path(a.review_out);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(reviewed,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(reviewed,indent=2))
    print("PDC_BIOSPECIMEN_AGGREGATES",json.dumps(x["studies"]))
    print("PDC_METADATA_CASE_INTERSECTION",json.dumps(x["proteome_phospho_metadata_intersection"]))
if __name__=="__main__":main()
