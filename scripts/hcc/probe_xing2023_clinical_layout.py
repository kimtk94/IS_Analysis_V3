#!/usr/bin/env python3
"""Privacy-preserving Xing2023 clinical-key orientation audit.

Diagnose whether the *internal* clinicopathological worksheet in mmc3 or the
separate mmc2 worksheets contain an identity-key layout compatible with the
protein sample header IDs. No phenotype, raw sample/case ID, abundance or
workbook bytes are emitted, and candidate overlaps never certify an identity.
"""
import argparse
import collections
import datetime as dt
import hashlib
import io
import json
import pathlib
import zipfile
import xml.etree.ElementTree as ET

import probe_pmc2026_cloud_proteomics as pmc
from probe_xing2023_clinical_id_join import (
    canon_id, get_file, match_clinical_ids, sheet_rows, trailing_code
)

EXPECTED_SHA256 = {
    "mmc3.xlsx": "7da732a933eba7a84936b8033d37564fd27125e29d059f9d7e306cf456965e17",
    "mmc2.xlsx": "38e409ee57eae5621176bfd6a5538c1ff94751725f550f5472997b4cca08826b",
}
MAX_ROWS = 250
MAX_COLS = 307


def paired_numeric_header_keys(protein_headers):
    pairs = collections.defaultdict(set)
    for label in protein_headers:
        classified = pmc.classify_patient_label(label)
        if classified:
            key = canon_id(classified[0])
            if key:
                pairs[key].add(classified[1])
    return {key for key, types in pairs.items() if types == {"T", "N"}}


def transpose_rows(rows):
    """Inspect horizontal patient tables by treating each source row as a field."""
    cols = sorted({col for _, values in rows for col in values})
    return [
        (i + 4, {rownum: values[col] for rownum, values in rows if col in values})
        for i, col in enumerate(cols)
    ]


def score_sheet(protein_headers, rows):
    """All overlap metrics are exploratory, including full numeric matches."""
    pair_keys = paired_numeric_header_keys(protein_headers)
    vertical = match_clinical_ids(protein_headers, rows)
    horizontal = match_clinical_ids(protein_headers, transpose_rows(rows))
    header_types = collections.defaultdict(set)
    prefixed_numeric_candidates = set()
    for _, fields in rows:
        for val in fields.values():
            parsed = pmc.classify_patient_label(val)
            if parsed:
                k = canon_id(parsed[0])
                if k in pair_keys:
                    header_types[k].add(parsed[1])
            k, _ = trailing_code(val)
            if k in pair_keys:
                prefixed_numeric_candidates.add(k)

    return {
        "rows_inspected": len(rows),
        "max_col_index_observed": max((c for _, values in rows for c in values), default=0),
        "pair_key_count_in_protein_header": len(pair_keys),
        "vertical_best_numeric_overlap": vertical["best_clinical_column_match_count"],
        "horizontal_best_numeric_overlap": horizontal["best_clinical_column_match_count"],
        "vertical_best_prefixed_overlap": vertical["prefixed_id_strategy"]["best_prefixed_id_code_matches"],
        "horizontal_best_prefixed_overlap": horizontal["prefixed_id_strategy"]["best_prefixed_id_code_matches"],
        "matching_numeric_t_header_keys_in_cells": sum("T" in x for x in header_types.values()),
        "matching_numeric_n_header_keys_in_cells": sum("N" in x for x in header_types.values()),
        "matching_numeric_tn_header_keys_in_cells": sum(x == {"T", "N"} for x in header_types.values()),
        "any_prefixed_numeric_candidates": len(prefixed_numeric_candidates) > 0,
        "clinical_field_identity_confirmed": False,
        "patient_level_join_confirmed": False,
        "interpretation": "Candidate column/row pattern overlap only; not proof of patient identity."
    }



def audit_clinical_row_identity_layout(protein_headers, rows):
    """Require sample-label T/N correspondence on the SAME metadata record."""
    pair_keys = paired_numeric_header_keys(protein_headers)
    keyed_rows = collections.defaultdict(set)
    keyed_types = collections.defaultdict(set)
    occurrence = collections.Counter()
    label_columns = collections.Counter()
    matched_data_rows = 0
    rows_one_key = 0
    rows_same_tn_key = 0
    rows_with_other_populated_fields = 0
    total_data_rows = 0
    for rownum, fields in rows:
        if rownum <= 3:
            continue
        total_data_rows += 1
        labels = []
        matching_cols = set()
        for col, value in fields.items():
            parsed = pmc.classify_patient_label(value)
            if not parsed:
                continue
            key = canon_id(parsed[0])
            if key not in pair_keys:
                continue
            tag = parsed[1]
            labels.append((key, tag))
            matching_cols.add(col)
            keyed_rows[key].add(rownum)
            keyed_types[key].add(tag)
            occurrence[(key, tag)] += 1
            label_columns[(col, tag)] += 1
        if not labels:
            continue
        matched_data_rows += 1
        keys = {key for key, _ in labels}
        if len(keys) == 1:
            rows_one_key += 1
            if {tag for _, tag in labels} == {"T", "N"}:
                rows_same_tn_key += 1
        if any(str(v).strip() for col, v in fields.items() if col not in matching_cols):
            rows_with_other_populated_fields += 1

    both_in_one_row = sum(
        keyed_types[k] == {"T", "N"} and len(keyed_rows[k]) == 1
        and occurrence[(k, "T")] == 1 and occurrence[(k, "N")] == 1
        for k in pair_keys
    )
    return {
        "clinical_data_rows_scanned": total_data_rows,
        "rows_with_protein_matching_labels": matched_data_rows,
        "rows_having_exactly_one_protein_patient_key": rows_one_key,
        "rows_with_same_key_tumor_and_normal_labels": rows_same_tn_key,
        "unique_protein_keys_with_exactly_one_T_and_one_N_in_same_row": both_in_one_row,
        "protein_keys_with_any_matching_clinical_label": len(keyed_rows),
        "matched_clinical_rows_with_non_label_fields": rows_with_other_populated_fields,
        "repeated_t_or_n_sample_label_count": sum(max(0, n - 1) for n in occurrence.values()),
        "matched_label_column_counts": [
            {"column_ordinal": col, "assay": tag, "n": count}
            for (col, tag), count in sorted(label_columns.items())
        ],
        "structure_consistent_for_all_protein_pairs": (
            len(pair_keys) == 152 and both_in_one_row == 152
            and matched_data_rows == 152 and rows_same_tn_key == 152
            and rows_with_other_populated_fields == 152
        ),
        "clinical_variable_definitions_verified": False,
        "clinically_eligible_patient_count_confirmed": False,
        "rna_protein_join_confirmed": False,
        "interpretation": "Same-row sample-label structural link, NOT independently reviewed clinical variable definitions or an RNA join."
    }

def workbook_sheet_count(z):
    root = ET.fromstring(z.read("xl/workbook.xml"))
    return len(root.find(pmc.NS + "sheets") or [])


def audit_from_workbooks(proteome_raw, clinical_raw):
    """Parse authorized public workbooks in memory; emit safe aggregate counts."""
    with zipfile.ZipFile(io.BytesIO(proteome_raw)) as zp, \
         zipfile.ZipFile(io.BytesIO(clinical_raw)) as zc:
        header_rows = sheet_rows(zp, 2, limit=2)
        if len(header_rows) < 2 or header_rows[1][0] != 2:
            raise ValueError("missing protein worksheet header")
        protein_headers = [
            value for col, value in sorted(header_rows[1][1].items()) if col >= 4
        ]
        if len(protein_headers) != 304:
            raise ValueError("unexpected number of proteome sample header cells")
        pair_count = len(paired_numeric_header_keys(protein_headers))
        if pair_count != 152:
            raise ValueError("proteome header T/N pair count changed")
        checks = []
        for name, archive, indices in (
            ("mmc3_internal", zp, (3,)),
            ("mmc2_external_clinical", zc, range(1, workbook_sheet_count(zc) + 1)),
        ):
            for index in indices:
                path = "xl/worksheets/sheet%d.xml" % index
                if path not in archive.namelist():
                    continue
                rows = sheet_rows(archive, index, limit=MAX_ROWS)
                record = {"source": name, "sheet_ordinal": index,
                          **score_sheet(protein_headers, rows)}
                if name == "mmc2_external_clinical" and index == 3:
                    record["row_identity_structure"] = audit_clinical_row_identity_layout(protein_headers, rows)
                checks.append(record)
    return {
        "protein_tn_header_pair_count": pair_count,
        "sheet_diagnostics": checks,
        "scan_limits": {"max_rows_per_sheet": MAX_ROWS, "max_columns": MAX_COLS},
        "patient_identifiers_or_phenotypes_logged": False,
        "protein_clinical_join_confirmed": False,
        "rna_protein_join_confirmed": False,
        "V02": "BLOCKED",
        "next_gate": "Verify official patient-key codebook and one-to-one mapping before clinical join promotion."
    }


def audit():
    result = {"status": "UNVERIFIED", "V02": "BLOCKED"}
    try:
        meta, _, _ = pmc.bounded(
            pmc.BASE + "metadata/" + pmc.PMCID + ".1.json", pmc.MAX_META, timeout=20
        )
        raw = {}
        for name in ("mmc3.xlsx", "mmc2.xlsx"):
            binary = get_file(meta, name)
            if hashlib.sha256(binary).hexdigest() != EXPECTED_SHA256[name]:
                raise ValueError("unexpected workbook SHA256")
            raw[name] = binary
        result.update(audit_from_workbooks(raw["mmc3.xlsx"], raw["mmc2.xlsx"]))
        result["status"] = "INTERNAL_AND_EXTERNAL_CLINICAL_LAYOUT_AUDITED"
    except Exception as exc:
        result["error_type"] = type(exc).__name__
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out",
        default="results/hcc/feasibility/xing2023_clinical_layout_audit.json"
    )
    args = parser.parse_args()
    output = {"checked_utc": dt.datetime.now(dt.timezone.utc).isoformat(), "audit": audit()}
    path = pathlib.Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(output, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
