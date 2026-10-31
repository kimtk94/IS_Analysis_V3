#!/usr/bin/env python3
"""Ng/Suter public processed RDS sample labels vs official Ng Supplement 1.

The 50 RDS labels are read from an ephemeral private file in the Actions runner.
No labels, patient IDs, rows, raw phenotype values or hashes are logged/uploaded.
Match strictly by normalized identity; do NOT fabricate fuzzy matches.
"""
import argparse
import collections
import datetime as dt
import hashlib
import io
import json
import pathlib
import re
import xml.etree.ElementTree as ET
import zipfile

from probe_ng_gao_public import (
    NG_URLS, fetch_bounded, _as_table_rows, _normalized_label, NS, MAX_ENTRY
)

OFFICIAL_SHA256="edf1d2afd3baaedd56edee8536c16f10e4a62cb6cdac929863a062671ee73d25"
SAMPLE_N=50
YES={"yes","y","true","1","x","available"}

def official_rows():
    raw=None
    for url in NG_URLS:
        try:
            data,_=fetch_bounded(url,limit=12*1024*1024,timeout=20)
            if hashlib.sha256(data).hexdigest()==OFFICIAL_SHA256:
                raw=data
                break
        except (OSError,ValueError):
            pass
    if raw is None:raise ValueError("source workbook not checksum-verifiable")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        shared=[]
        if "xl/sharedStrings.xml" in z.namelist():
            meta=z.getinfo("xl/sharedStrings.xml")
            if meta.file_size>MAX_ENTRY:raise ValueError("shared strings too large")
            strings=ET.fromstring(z.read("xl/sharedStrings.xml"))
            shared=["".join(x.text or "" for x in si.iter(NS+"t"))
                    for si in strings.findall(NS+"si")]
        for name in z.namelist():
            if not re.fullmatch(r"xl/worksheets/sheet\d+\.xml",name):
                continue
            if z.getinfo(name).file_size>MAX_ENTRY:continue
            doc=ET.fromstring(z.read(name))
            cells=doc.findall(".//"+NS+"sheetData/"+NS+"row")
            rows=_as_table_rows(cells,shared)
            for idx,row in enumerate(rows[:6]):
                columns={_normalized_label(v):col for col,v in row.items() if v}
                if {"tumorbiopsyid","patientid","rnaseq","proteomics"}.issubset(columns):
                    return rows[idx+1:],columns
    raise ValueError("no official sample metadata sheet")

def compare(private_labels,rows,columns):
    if len(private_labels)!=SAMPLE_N or len(set(private_labels))!=SAMPLE_N:
        raise ValueError("RDS sample axis is not 50 unique labels")
    if any(not re.fullmatch(r"[a-z0-9]+",x) for x in private_labels):
        raise ValueError("invalid normalized sample label")
    src={}
    duplicate=0
    etiology_by_biopsy={}
    eti_col=columns.get("underlyingliverdisease")
    for row in rows:
        tumor=_normalized_label(row.get(columns["tumorbiopsyid"],""))
        patient=_normalized_label(row.get(columns["patientid"],""))
        if not tumor:continue
        if eti_col is not None:
            etiology_by_biopsy[tumor]=str(row.get(eti_col,""))
        flags={k:str(row.get(columns[k],"")).strip().lower() in YES
               for k in ("rnaseq","proteomics")}
        if tumor in src:
            duplicate+=1
            if src[tumor]!=(patient,flags):raise ValueError("source alias collision")
        else:src[tumor]=(patient,flags)
    matched=[(key,src[key]) for key in private_labels if key in src]
    # Per-patient source-reported liver disease labels are preclinical candidate
    # grouping only. Overlapping etiologies are not merged into single-exposure.
    patient_source_group={}
    def categories(raw):
        value=raw.lower()
        found=set()
        if re.search(r"\bhcv\b|hepatitis c",value):found.add("HCV")
        if re.search(r"\bhbv\b|hepatitis b",value):found.add("HBV")
        if re.search(r"\bald\b|alcohol",value):found.add("ALD")
        if re.search(r"\bnafld\b|\bnash\b|steato|metabolic",value):found.add("NAFLD")
        return tuple(sorted(found)) if found else ("SOURCE_OTHER_OR_UNCLASSIFIED",)
    for key,(patient,flags) in matched:
        if patient:
            patient_source_group.setdefault(patient,set()).add(categories(etiology_by_biopsy.get(key,"")))
    group_counts=collections.Counter()
    multi_biopsy_label_conflict=0
    for patient,groups in patient_source_group.items():
        if len(groups)!=1:
            multi_biopsy_label_conflict+=1
            group_counts["CONFLICTING_SOURCE_LABELS"]+=1
        else:
            k="+".join(next(iter(groups)))
            group_counts[k]+=1
    presence={}
    for disease in ("ALD","HCV","HBV","NAFLD"):
        n=sum(disease in next(iter(groups)) for groups in patient_source_group.values() if len(groups)==1)
        presence[disease]=n if n>=5 else None
    safe_groups={key:(count if count>=5 else None)
                 for key,count in sorted(group_counts.items())}

    all_both=sum(all(flags.values()) for _,(_,flags) in matched)
    unique_patients={pat for _,(pat,_) in matched if pat}
    return {
        "public_source_sample_alias_n":len(src),
        "source_duplicate_same_alias_row_n":duplicate,
        "processed_rds_sample_label_n":len(private_labels),
        "exact_normalized_biopsy_label_matched_n":len(matched),
        "RDS_labels_unmatched_n":SAMPLE_N-len(matched),
        "matched_biopsies_with_both_assay_flags_n":all_both,
        "distinct_nonempty_matched_patient_labels_n":len(unique_patients),
        "matched_biopsies_without_patient_label_n":
            sum(not patient for _,(patient,_) in matched),
        "published_suter_biopsy_n":50,
        "published_suter_unique_patient_n":48,
        "source_patient_multiplicity_consistent_with_suter":
            len(matched)==50 and len(unique_patients)==48,
        "strict_same_biopsy_source_label_match":
            len(matched)==50 and all_both==50 and len(unique_patients)==48,
        "source_etiology_candidate_patients_per_group_suppress_lt5":safe_groups,
        "source_etiology_candidate_any_exposure_suppress_lt5":presence,
        "source_etiology_small_cell_suppressed_count":sum(v<5 for v in group_counts.values()),
        "conflicting_etiology_labels_across_biopsies_n":multi_biopsy_label_conflict,
        "etiology_clinical_codebook_verified":False,
        "source_TP_quant_matrix_same_sample_names":True,
        "scientific_etiology_coding_confirmed":False,
        "prognosis_event_coding_confirmed":False,
        "source_raw_patient_ids_or_clinical_values_logged":False
    }

def audit(private_path):
    result={"status":"UNVERIFIED","V02":"BLOCKED","sample_ids_or_patient_rows_in_output":False}
    try:
        p=pathlib.Path(private_path)
        if not p.is_file() or p.stat().st_size>20000:
            raise ValueError("private RDS label file missing or oversized")
        labels=[_normalized_label(line.strip()) for line in
                p.read_text(encoding="utf-8").splitlines() if line.strip()]
        rows,columns=official_rows()
        r=compare(labels,rows,columns)
        result.update(r)
        result["source_sha256_verified"]=True
        result["status"]=("EXACT_SOURCE_SAMPLE_ALIAS_AND_PATIENT_JOIN_PRECHECK"
                          if r["strict_same_biopsy_source_label_match"]
                          else "PARTIAL_OR_UNMATCHED_SAMPLE_ALIASES_HOLD")
        # A source label match is not a reviewed batch, normalization,
        # treatment, etiology or clinical-event codebook.
        result["V02"]="PARTIAL_SOURCE_ALIAS_ONLY" if r["strict_same_biopsy_source_label_match"] else "BLOCKED"
    except Exception as e:
        result["error_type"]=type(e).__name__
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--private-labels",required=True)
    p.add_argument("--out",required=True)
    args=p.parse_args()
    report={"checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
            "audit":audit(args.private_labels)}
    out=pathlib.Path(args.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2,ensure_ascii=False))
    if report["audit"]["status"]=="UNVERIFIED":
        raise SystemExit(1)
if __name__=="__main__":main()
