#!/usr/bin/env python3
"""Xing 2023 etiology-label feasibility, with source-checked population N=152.

This is a *codebook preflight*, not HBV/HCV/ALD/MASH disease classification.
Only numerical summaries / allowlisted group counts are released. No patient
IDs, original XLSX, uncensored individual field values or raw category terms.
"""
import argparse
import collections
import datetime as dt
import hashlib
import io
import json
import pathlib
import re
import statistics
import zipfile

import probe_pmc2026_cloud_proteomics as pmc
from probe_xing2023_clinical_id_join import get_file, sheet_rows
from probe_xing2023_clinical_layout import EXPECTED_SHA256
from probe_xing2023_clinical_completeness import (
    inspect_availability, strict_rows, classify_missing, finite_numeric, safe_label
)

COLS = {"candidate_HBV_related":8, "candidate_HCV_related":13,
        "candidate_underlying_liver_disease":16}
MIN_CATEGORY_REPORT_N = 5

def numeric_shape(rows, col):
    strings = [str(row.get(col, "")).strip() for row in rows]
    actual = [float(v) for v in strings
              if not classify_missing(v) and finite_numeric(v)]
    unique = set(actual)
    bins = {0.0,1.0}
    return {
        "present_nonmissing":sum(not classify_missing(v) for v in strings),
        "finite_numeric":len(actual),
        "missing_or_unclassified":len(strings)-len(actual),
        "unique_numeric_value_count":len(unique),
        "all_finite_values_binary_0_1":bool(actual) and unique.issubset(bins),
        "numeric_zero_count":sum(v==0 for v in actual),
        "numeric_one_count":sum(v==1 for v in actual),
        "numeric_nonbinary_count":sum(v not in bins for v in actual),
        "nonbinary_range_summary":(
            {"minimum":min(actual),"median":statistics.median(actual),"maximum":max(actual)}
            if actual and not unique.issubset(bins) else None),
        "clinical_positive_negative_classification_verified":False,
    }

def classify_literal_text(value):
    """Allow only explicit clinical *literal* group labels, not e.g. cirrhosis."""
    s=re.sub(r"\s+"," ",str(value or "").strip().upper())
    if classify_missing(s):
        return "MISSING"
    # Never serialize a raw clinical phrase or infer disease from a lab titre.
    markers=set()
    if re.search(r"\b(?:HBV|HEPATITIS B)\b",s):markers.add("HBV_LITERAL")
    if re.search(r"\b(?:HCV|HEPATITIS C)\b",s):markers.add("HCV_LITERAL")
    if re.search(r"\b(?:ALD|ALCOHOL(?:IC)?(?:[- ]RELATED)?(?: LIVER DISEASE)?)\b",s):markers.add("ALD_LITERAL")
    if re.search(r"\b(?:MASH|NASH|NAFLD|MASLD|STEATOHEPATITIS)\b",s):markers.add("METABOLIC_LITERAL")
    if len(markers)>1:return "MULTIPLE_LITERAL"
    if len(markers)==1:return next(iter(markers))
    return "UNCLASSIFIED_OR_OTHER"

def summarize_etiology(matched, headers):
    rows=list(matched.values())
    if len(rows)!=152:raise ValueError("expected 152 source-matched sample records")
    # Reject source layout change, but no lab assay can be interpreted as clinical
    # positivity merely because header mentions HBV or HCV.
    if "etiology_HBV" not in safe_label(headers[8]):
        raise ValueError("HBV-associated header no longer found")
    if "etiology_HCV" not in safe_label(headers[13]):
        raise ValueError("HCV-associated header no longer found")
    markers={name:numeric_shape(rows,col) for name,col in COLS.items() if col in (8,13)}
    joint=sum(all(not classify_missing(row.get(col,"")) and finite_numeric(row.get(col,""))
                  for col in (8,13)) for row in rows)
    categories=collections.Counter(classify_literal_text(row.get(16,"")) for row in rows)
    if sum(categories.values())!=152:raise ValueError("inconsistent category coverage")
    return {
        "matched_protein_to_clinical_rows":152,
        "candidate_source_columns":COLS,
        "marker_numeric_shapes":markers,
        "joint_two_marker_numeric_coverage_n":joint,
        "underlying_liver_disease_text_field":{
            "source_col":16,
            "n":152,
            "category_count_with_small_cells_suppressed":{
                key:(val if val>=MIN_CATEGORY_REPORT_N else None)
                for key,val in sorted(categories.items())
            },
            "number_of_suppressed_categories":sum(v<MIN_CATEGORY_REPORT_N for v in categories.values()),
            "raw_categories_logged":False
        },
        "etiology_patient_group_assignment_confirmed":False,
        "HBV_HCV_titer_thresholds_confirmed":False,
        "HBV_HCV_ALD_MASH_group_counts_confirmed":False,
        "etiology_interaction_power_calculable":False,
        "reason":"A detectable HBV/HCV worksheet header and numeric completeness do not establish infection status, etiology, or group power.",
        "source_patient_identifiers_or_values_logged":False,
        "V04":"PARTIAL",
        "V02":"BLOCKED"
    }

def audit():
    out={"status":"UNVERIFIED","V04":"PARTIAL","V02":"BLOCKED"}
    try:
        metadata,_,_=pmc.bounded(pmc.BASE+"metadata/"+pmc.PMCID+".1.json",
                                 pmc.MAX_META,timeout=20)
        workbooks={}
        for name in ("mmc2.xlsx","mmc3.xlsx"):
            raw=get_file(metadata,name)
            if hashlib.sha256(raw).hexdigest()!=EXPECTED_SHA256[name]:
                raise ValueError("source workbook checksum changed")
            workbooks[name]=raw
        with zipfile.ZipFile(io.BytesIO(workbooks["mmc2.xlsx"])) as c, \
             zipfile.ZipFile(io.BytesIO(workbooks["mmc3.xlsx"])) as p:
            rawheaders=sheet_rows(p,2,limit=2)
            protein_headers=[v for idx,v in sorted(rawheaders[1][1].items()) if idx>=4]
            clinical=sheet_rows(c,3,limit=200)
            info=inspect_availability(protein_headers,clinical)
            if not info["matching_audit_passed"] or info["matched_sample_label_pairs"]!=152:
                raise ValueError("source clinical/protein matching failed")
            headers={idx:" ".join(row.get(idx,"") for num,row in clinical if num<=2)
                     for idx in COLS.values()}
            result=summarize_etiology(strict_rows(protein_headers,clinical),headers)
        out.update(result)
        out["source_checksums_verified"]=True
        out["status"]="ETIOLOGY_CODEBOOK_PREFLIGHT_AUDITED"
    except Exception as exc:
        out["error_type"]=type(exc).__name__
    return out

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",default="results/hcc/feasibility/xing2023_etiology_codebook_preflight.json")
    args=parser.parse_args()
    doc={"checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),"audit":audit()}
    out=pathlib.Path(args.out)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(doc,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(doc,ensure_ascii=False,indent=2))
    if doc["audit"]["status"]!="ETIOLOGY_CODEBOOK_PREFLIGHT_AUDITED":
        raise SystemExit(1)

if __name__=="__main__":
    main()
