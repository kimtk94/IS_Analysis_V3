#!/usr/bin/env python3
"""Xing 2023: privacy-safe clinical covariate availability on matched proteomic pairs.

Checks official PMC-hosted mmc2/mmc3 SHA256 and 152/152 T/N row keys before
counting variable completeness. Never exports patient identifiers, raw clinical
values, source workbooks or row-wise sensitive results. Column semantic labels
are allowlisted patterns; endpoint meaning is NOT assumed from a column name.
"""
import argparse
import collections
import datetime as dt
import hashlib
import io
import json
import pathlib
import re
import zipfile

import probe_pmc2026_cloud_proteomics as pmc
from probe_xing2023_clinical_id_join import get_file, sheet_rows, canon_id
from probe_xing2023_clinical_layout import (
    EXPECTED_SHA256, paired_numeric_header_keys, audit_clinical_row_identity_layout
)

# Labels detected from publication worksheet headers, NOT endpoint codebook verification.
LABEL_PATTERNS = (
    ("RFS", r"\brfs\b|recurren(?:ce|t)|relapse|disease[- ]free"),
    ("OS", r"\bos\b|overall\s+survival|survival\s+(?:time|status)"),
    ("DSS", r"\bdss\b|disease[- ]specific\s+survival"),
    ("PFS", r"\bpfs\b|progression[- ]free"),
    ("follow_up", r"follow[- ]?up"),
    ("survival_duration", r"\btime\b|\bmonth\b|\byear\b|\bduration\b"),
    ("event_status", r"\bevent\b|\bstatus\b|\bdeath\b|\bcens"),
    ("etiology_HBV", r"\bhbv\b|hbsag|hepatitis\s*b"),
    ("etiology_HCV", r"\bhcv\b|hepatitis\s*c"),
    ("etiology_MASH", r"\bmash\b|\bnafld\b|\bnash\b|steato"),
    ("etiology_ALD", r"\bald\b|alcohol|ethanol"),
    ("etiology_other", r"\betiology\b|aetiology|cirrhosis|liver\s+disease"),
    ("BCLC", r"\bbclc\b"),
    ("TNM", r"\btnm\b"),
    ("stage", r"\bstage\b"),
    ("tumor_grade", r"\bgrade\b|differentiat"),
    ("MVI", r"\bmvi\b|vascular\s+invasion"),
    ("AFP", r"\bafp\b|alpha[- ]fetoprotein"),
    ("age", r"\bage\b"),
    ("sex", r"\bsex\b|\bgender\b"),
    ("therapy", r"treat|therap|resect|surgery|ablat|sorafenib|chemo"),
    ("tumor_size", r"tumou?r\s+size|tumou?r\s+diameter"),
    ("tumor_number", r"tumou?r\s+number|tumou?r\s+count|multiple\s+tumou?r"),
    ("subtype", r"\bsubtype\b|subgroup|cluster"),
)
MISSING_MARKERS = {"", "na", "n/a", "nan", "none", "null", ".", "-", "--",
                   "unknown", "not available", "not reported", "missing"}
HEADER_ROW_CEILING = 2
AUDIT_SOURCE = "PMC10783603 mmc2.xlsx sheet3; verified SHA256"

def safe_label(value):
    """Return only canonical family names, never copy a patient/source string."""
    s = re.sub(r"[\r\n\t]+", " ", str(value or "")).lower().strip()
    if not s or len(s) > 140 or pmc.classify_patient_label(s):
        return []
    return sorted({name for name, regex in LABEL_PATTERNS if re.search(regex, s, re.I)})

def classify_missing(raw):
    val = str(raw or "").strip()
    if val.lower() in MISSING_MARKERS:
        return True
    if val in ("#N/A", "#VALUE!", "#REF!", "#DIV/0!"):
        return True
    return False

def finite_numeric(raw):
    import math
    try:
        return math.isfinite(float(str(raw).strip()))
    except (TypeError, ValueError, OverflowError):
        return False

def strict_rows(protein_headers, rows):
    """Make a unique one-to-one row map; do not export mapping or row indices."""
    pair_keys = paired_numeric_header_keys(protein_headers)
    found = {}
    for _, cells in rows:
        t, n = cells.get(1, ""), cells.get(2, "")
        p_t = pmc.classify_patient_label(t)
        p_n = pmc.classify_patient_label(n)
        if not p_t or not p_n or p_t[1] != "T" or p_n[1] != "N":
            continue
        kt, kn = canon_id(p_t[0]), canon_id(p_n[0])
        if kt is None or kt != kn or kt not in pair_keys:
            continue
        if kt in found:
            raise ValueError("duplicate source sample key")
        found[kt] = cells
    if len(pair_keys) != 152 or len(found) != 152:
        raise ValueError("152 paired sample identities not established")
    if any(not any(str(v).strip() for c, v in cells.items() if c > 2)
           for cells in found.values()):
        raise ValueError("matched row has no clinical variables")
    return found

def inspect_availability(protein_headers, clinical_rows):
    layout = audit_clinical_row_identity_layout(protein_headers, clinical_rows)
    if not layout["structure_consistent_for_all_protein_pairs"]:
        raise ValueError("T/N sample-to-clinical structural mapping fails")
    matched = strict_rows(protein_headers, clinical_rows)
    header_values = {c: " ".join(str(fields.get(c, "")) for _, fields in clinical_rows
                                   if _ <= HEADER_ROW_CEILING).strip()
                     for c in range(3, 34)}
    metrics = []
    for col in range(3, 34):
        values = [fields.get(col, "") for fields in matched.values()]
        n_present = sum(not classify_missing(x) for x in values)
        n_numeric = sum(not classify_missing(x) and finite_numeric(x) for x in values)
        names = safe_label(header_values[col])
        metrics.append({
            "excel_column": col,
            "detected_header_families": names,
            "header_recognized": bool(names),
            "n_matched": 152,
            "n_nonmissing": n_present,
            "n_missing": 152-n_present,
            "n_numeric": n_numeric,
            "missingness_fraction": round((152-n_present)/152, 4),
            # No raw categorical value counts: codebook and cell semantics unknown.
            "variable_definition_reviewed": False,
            "endpoint_or_covariate_ready": False,
        })
    # Completeness only; not whether follow-up, event and censoring are valid.
    return {
        "source": AUDIT_SOURCE,
        "matched_sample_label_pairs": len(matched),
        "matching_audit_passed": True,
        "clinical_variable_column_count_audited": len(metrics),
        "matched_rows_with_any_covar_value": sum(
            any(not classify_missing(fields.get(c, "")) for c in range(3, 34))
            for fields in matched.values()
        ),
        "header_detection_scope": "literal workbook rows 1-2; allowlisted concepts only",
        "variables": metrics,
        "recognized_families": sorted({k for m in metrics for k in m["detected_header_families"]}),
        "endpoint_event_counts": None,
        "recurrence_or_OS_model_ready": False,
        "etiology_group_counts": None,
        "etiology_interaction_ready": False,
        "clinical_codebook_verified": False,
        "raw_patient_identifiers_or_phenotypes_logged": False,
        "rna_protein_quantitative_join": False,
        "V02": "BLOCKED",
        "interpretation": ("152 source protein T/N labels map to clinical rows. "
                           "Completeness only; category meaning, recurrence time/status "
                           "and clinical endpoints need source-codebook confirmation.")
    }

def audit():
    out = {"status": "UNVERIFIED", "V02": "BLOCKED",
           "original_xlsx_in_outputs": False}
    try:
        meta, _, _ = pmc.bounded(
            pmc.BASE + "metadata/" + pmc.PMCID + ".1.json",
            pmc.MAX_META, timeout=20)
        raw = {}
        for name in ("mmc2.xlsx", "mmc3.xlsx"):
            binary = get_file(meta, name)
            if hashlib.sha256(binary).hexdigest() != EXPECTED_SHA256[name]:
                raise ValueError("official workbook SHA mismatch")
            raw[name] = binary
        with zipfile.ZipFile(io.BytesIO(raw["mmc2.xlsx"])) as zc, \
             zipfile.ZipFile(io.BytesIO(raw["mmc3.xlsx"])) as zp:
            prot_rows = sheet_rows(zp, 2, limit=2)
            if len(prot_rows) < 2 or prot_rows[1][0] != 2:
                raise ValueError("protein header missing")
            protein_headers = [v for c, v in sorted(prot_rows[1][1].items()) if c >= 4]
            clin = sheet_rows(zc, 3, limit=200)
            result = inspect_availability(protein_headers, clin)
        out.update(result)
        out["status"] = "MATCHED_CLINICAL_COLUMN_COMPLETENESS_AUDITED"
    except Exception as exc:
        out["error_type"] = type(exc).__name__
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/hcc/feasibility/xing2023_clinical_completeness.json")
    args = ap.parse_args()
    doc = {"checked_utc": dt.datetime.now(dt.timezone.utc).isoformat(), "audit": audit()}
    dest = pathlib.Path(args.out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(doc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(doc, ensure_ascii=False, indent=2))
    # Do not silently label an inaccessible or changed source as a complete audit.
    if doc["audit"]["status"] != "MATCHED_CLINICAL_COLUMN_COMPLETENESS_AUDITED":
        raise SystemExit(1)

if __name__ == "__main__":
    main()
