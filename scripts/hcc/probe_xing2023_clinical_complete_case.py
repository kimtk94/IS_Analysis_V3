#!/usr/bin/env python3
"""Privacy-safe HCC Xing2023 complete-case eligibility audit.

Read PMC Cloud mmc2/mmc3 with pinned SHA256, match 152 T/N protein/clinical
sample IDs in RAM, aggregate candidate covariate complete-case counts only.
No row identifiers, raw phenotype values or original XLSX files are emitted.
Candidate RFS / OS numeric flags are NOT confirmed survival outcomes.
"""
import argparse
import datetime as dt
import hashlib
import io
import json
import pathlib
import zipfile

import probe_pmc2026_cloud_proteomics as pmc
from probe_xing2023_clinical_id_join import get_file, sheet_rows
from probe_xing2023_clinical_layout import EXPECTED_SHA256
from probe_xing2023_clinical_completeness import (
    strict_rows, classify_missing, finite_numeric, inspect_availability,
    safe_label
)

SAMPLE_COUNT = 152
# Exact source worksheet columns, *not* proof of semantics/valid category codes.
COLUMN = {
    "age": 3, "sex": 4, "AFP": 5, "HBV_marker": 8,
    "HCV_marker": 13, "MVI": 15, "tumor_size": 21,
    "TNM_descriptor": 22, "TNM_stage": 23, "BCLC": 24,
    "RFS_candidate_code": 29, "RFS_candidate_duration": 31,
    "survival_candidate_code": 32, "survival_candidate_duration": 33
}
# For physical numeric fields strict finite numeric required; categorical
# fields are presence-only until the clinical source codebook is checked.
NUMERIC_REQUIRED = {"age", "AFP", "HBV_marker", "HCV_marker",
                    "tumor_size", "RFS_candidate_code",
                    "RFS_candidate_duration", "survival_candidate_code",
                    "survival_candidate_duration"}
BINARY_REQUIRED = {"RFS_candidate_code", "survival_candidate_code"}
POSITIVE_REQUIRED = {"RFS_candidate_duration", "survival_candidate_duration"}

PROFILE_DEFS = {
    "endpoint_only_RFS": ["RFS_candidate_code", "RFS_candidate_duration"],
    "endpoint_only_survival": ["survival_candidate_code", "survival_candidate_duration"],
    "baseline_RFS": ["age", "sex", "tumor_size", "BCLC",
                     "RFS_candidate_code", "RFS_candidate_duration"],
    "baseline_survival": ["age", "sex", "tumor_size", "BCLC",
                          "survival_candidate_code", "survival_candidate_duration"],
    "clinical_extended_RFS": ["age", "sex", "AFP", "MVI", "tumor_size",
                              "TNM_stage", "BCLC", "RFS_candidate_code",
                              "RFS_candidate_duration"],
    "clinical_extended_survival": ["age", "sex", "AFP", "MVI", "tumor_size",
                                   "TNM_stage", "BCLC", "survival_candidate_code",
                                   "survival_candidate_duration"],
    "viral_marker_coverage": ["HBV_marker", "HCV_marker"],
    "clinical_extended_viral_RFS": ["age", "sex", "AFP", "HBV_marker",
                                   "HCV_marker", "MVI", "tumor_size",
                                   "TNM_stage", "BCLC",
                                   "RFS_candidate_code", "RFS_candidate_duration"],
    "clinical_extended_viral_survival": [
        "age", "sex", "AFP", "HBV_marker", "HCV_marker", "MVI",
        "tumor_size", "TNM_stage", "BCLC",
        "survival_candidate_code", "survival_candidate_duration"
    ],
}

def valid_value(row, field):
    value = row.get(COLUMN[field], "")
    if classify_missing(value):
        return False
    if field in NUMERIC_REQUIRED:
        if not finite_numeric(value):
            return False
        f = float(value)
        if field in BINARY_REQUIRED and f not in (0.0, 1.0):
            return False
        if field in POSITIVE_REQUIRED and f <= 0:
            return False
    return True

def profile_audit(matched):
    """Aggregate only; never include subject identifiers or individual rows."""
    rows = list(matched.values())
    n = len(rows)
    if n != SAMPLE_COUNT:
        raise ValueError("unexpected sample denominator")
    fields = {
        name: {"valid_n": sum(valid_value(row, name) for row in rows),
               "invalid_or_missing_n": sum(not valid_value(row, name) for row in rows),
               "validation_kind": ("binary_numeric_only_not_interpreted"
                                   if name in BINARY_REQUIRED
                                   else "positive_finite_numeric"
                                   if name in POSITIVE_REQUIRED
                                   else "finite_numeric"
                                   if name in NUMERIC_REQUIRED
                                   else "not_missing_only")}
        for name in COLUMN
    }
    profiles = {}
    for name, selected in PROFILE_DEFS.items():
        rows_ok = [row for row in rows if all(valid_value(row, f) for f in selected)]
        profiles[name] = {
            "candidate_complete_case_n": len(rows_ok),
            "excluded_n": n-len(rows_ok),
            "field_count": len(selected),
            "field_names": selected,
            "scientifically_eligible": False
        }
    # Restricted to *numeric code values*, NOT event outcomes. This checks
    # potential collider/completeness bias before using richer covariates.
    for endpoint, col in (("RFS", COLUMN["RFS_candidate_code"]),
                          ("survival", COLUMN["survival_candidate_code"])):
        for profile_name in (f"clinical_extended_{endpoint}",
                             f"clinical_extended_viral_{endpoint}"):
            selected = PROFILE_DEFS[profile_name]
            cc = [row for row in rows if all(valid_value(row, field) for field in selected)]
            code0 = sum(float(row[col]) == 0 for row in cc)
            code1 = sum(float(row[col]) == 1 for row in cc)
            profiles[profile_name]["uninterpreted_event_code_0_n"] = code0
            profiles[profile_name]["uninterpreted_event_code_1_n"] = code1
    return {
        "n_matched_structurally": n,
        "fields": fields,
        "profiles": profiles,
        "prognostic_coding_reviewed": False,
        "biological_etiology_assignments_reviewed": False,
        "survival_model_allowed": False,
        "notes": [
            "Complete-case counts reflect missing-value and primitive numeric constraints, not codebook-valid categories.",
            "0 and 1 are raw candidate endpoint codes, not established recurrence/death events.",
            "These profiles cannot be used for powered RFS/OS claims before event semantics and confounders are verified."
        ],
        "patient_ids_or_values_logged": False
    }

def audit():
    output = {"status": "UNVERIFIED", "V02": "BLOCKED", "V12": "PARTIAL"}
    try:
        meta, _, _ = pmc.bounded(
            pmc.BASE + "metadata/" + pmc.PMCID + ".1.json",
            pmc.MAX_META, timeout=20
        )
        raw = {}
        for name in ("mmc2.xlsx", "mmc3.xlsx"):
            binary = get_file(meta, name)
            if hashlib.sha256(binary).hexdigest() != EXPECTED_SHA256[name]:
                raise ValueError("unexpected official workbook checksum")
            raw[name] = binary
        with zipfile.ZipFile(io.BytesIO(raw["mmc2.xlsx"])) as zc, \
             zipfile.ZipFile(io.BytesIO(raw["mmc3.xlsx"])) as zp:
            protein_rows = sheet_rows(zp, 2, limit=2)
            if len(protein_rows) != 2 or protein_rows[1][0] != 2:
                raise ValueError("proteome header missing")
            protein_headers = [v for c, v in sorted(protein_rows[1][1].items()) if c >= 4]
            clinical_rows = sheet_rows(zc, 3, limit=200)
            # Fail closed if source structural pairing changed.
            availability = inspect_availability(protein_headers, clinical_rows)
            if not availability["matching_audit_passed"] or \
               availability["matched_sample_label_pairs"] != SAMPLE_COUNT:
                raise ValueError("clinical source alignment changed")
            matched = strict_rows(protein_headers, clinical_rows)
            # Assert recognizable source header context for core numeric and
            # viral fields, while avoiding raw-header recording.
            headers = {
                col: " ".join(str(row.get(col, "")) for ord_, row in clinical_rows
                              if ord_ <= 2)
                for col in COLUMN.values()
            }
            expected = {
                "age": "age", "AFP": "AFP", "HBV_marker": "etiology_HBV",
                "HCV_marker": "etiology_HCV", "tumor_size": "tumor_size",
                "RFS_candidate_code": "RFS"
            }
            for field, family in expected.items():
                if family not in safe_label(headers[COLUMN[field]]):
                    raise ValueError("expected source column header family missing")
            summary = profile_audit(matched)
        output.update(summary)
        output["status"] = "MATCHED_COMPLETE_CASE_PRESCREEN_AUDITED"
        output["source_data_sha256_verified"] = True
        output["source_publication"] = "PMC10783603 Xing 2023"
        output["event_semantics_source_verified"] = False
        output["rna_protein_patient_quant_join_verified"] = False
    except Exception as exc:
        output["error_type"] = type(exc).__name__
    return output

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/hcc/feasibility/xing2023_clinical_complete_case.json")
    args = ap.parse_args()
    doc = {"checked_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "audit": audit()}
    p = pathlib.Path(args.out)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(doc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(doc, indent=2, ensure_ascii=False))
    if doc["audit"]["status"] != "MATCHED_COMPLETE_CASE_PRESCREEN_AUDITED":
        raise SystemExit(1)

if __name__ == "__main__":
    main()
