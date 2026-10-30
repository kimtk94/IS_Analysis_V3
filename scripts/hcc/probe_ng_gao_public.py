#!/usr/bin/env python3
"""HCC Ng/Gao public-source probe: bounded supplement inspection + public PDC metadata.
No EGA, raw reads, restricted genomes, patient-level output, or fabricated joins.
"""
import argparse
import datetime as dt
import hashlib
import io
import json
import pathlib
import re
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

MAX_DOWNLOAD = 12 * 1024 * 1024
MAX_ENTRY = 18 * 1024 * 1024
SUPPLEMENT = "41467_2022_29960_MOESM4_ESM.xlsx"
NG_URLS = (
    "https://pmc.ncbi.nlm.nih.gov/articles/instance/9068765/bin/" + SUPPLEMENT,
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC9068765/bin/" + SUPPLEMENT,
    "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-022-29960-8/MediaObjects/" + SUPPLEMENT,
)
PDC_ENDPOINT = "https://pdc.cancer.gov/graphql"
PDC_QUERY = (
    '{ proteome: study(pdc_study_id: "PDC000198") { pdc_study_id study_name cases_count aliquots_count analytical_fraction } '
    'phosphoproteome: study(pdc_study_id: "PDC000199") { pdc_study_id study_name cases_count aliquots_count analytical_fraction } }'
)
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
REL = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def fetch_bounded(url, post_data=None, timeout=18, limit=MAX_DOWNLOAD):
    headers = {"User-Agent": "HCC-Multiomics-Feasibility/2.0 (public supplementary metadata)",
               "Accept": "application/octet-stream, application/json, */*"}
    if post_data is not None:
        headers["Content-Type"] = "application/json"
    request = urllib.request.Request(url, data=post_data, headers=headers,
                                     method="POST" if post_data is not None else "GET")
    with urllib.request.urlopen(request, timeout=timeout) as response:
        data = response.read(limit + 1)
        status = response.status
    if len(data) > limit:
        raise ValueError("max download exceeded; refusing source")
    return data, status


def _cell_value(cell, shared):
    kind = cell.attrib.get("t")
    if kind == "inlineStr":
        return "".join(node.text or "" for node in cell.iter(NS + "t"))
    raw = cell.find(NS + "v")
    value = raw.text if raw is not None else ""
    if kind == "s" and value:
        i = int(value)
        return shared[i] if i < len(shared) else ""
    return value or ""


def _normalized_label(s):
    return re.sub(r"[^a-z0-9]+", "", str(s).lower())


def _column_index(col):
    ix = 0
    for ch in col:
        ix = ix * 26 + (ord(ch) - ord("A") + 1)
    return ix - 1


def _as_table_rows(rows, shared):
    """Build a sparse XLSX sheet with A1/cell reference alignment."""
    result = []
    for row in rows:
        values = {}
        for fallback, cell in enumerate(row.findall(NS + "c")):
            m = re.match(r"([A-Z]+)", cell.attrib.get("r", ""))
            idx = _column_index(m.group(1)) if m else fallback
            if idx < 500:
                values[idx] = _cell_value(cell, shared)
        result.append(values)
    return result


def _metadata_pair_summary(rows):
    """Summarize SAME-BIOPSY flags, not expression-matrix validation."""
    header_row = None
    mapping = {}
    for i, row in enumerate(rows[:6]):
        names = {_normalized_label(v): col for col, v in row.items() if v}
        if {"tumorbiopsyid", "patientid", "rnaseq", "proteomics"}.issubset(names):
            header_row = i
            mapping = names
            break
    if header_row is None:
        return {"status": "HEADER_NOT_FOUND", "matrix_join_verified": False}
    keys = ("rnaseq", "proteomics", "phosphoproteomics")
    hist = {key: {} for key in keys}
    true_values = {"yes", "y", "true", "1", "x", "available"}
    false_values = {"no", "n", "false", "0", "na", "n/a", "-", ""}
    tumor_rows = 0
    patients = set()
    tumors = set()
    flagged = []
    unknown_flags = {}
    etiology_all = {}
    etiology_paired = {}
    stage_paired = {}
    for row in rows[header_row + 1:]:
        tumor = row.get(mapping["tumorbiopsyid"], "").strip()
        patient = row.get(mapping["patientid"], "").strip()
        if not tumor:
            continue
        tumor_rows += 1
        disease = row.get(mapping["underlyingliverdisease"], "").strip() if "underlyingliverdisease" in mapping else ""
        stage = row.get(mapping["bclc"], "").strip() if "bclc" in mapping else ""
        etiology_all[disease] = etiology_all.get(disease, 0) + 1
        tumors.add(tumor)
        if patient:
            patients.add(patient)
        statuses = {}
        for key in keys:
            v = row.get(mapping[key], "").strip() if key in mapping else ""
            hist[key][v] = hist[key].get(v, 0) + 1
            if v.lower() in true_values:
                statuses[key] = True
            elif v.lower() in false_values:
                statuses[key] = False
            else:
                statuses[key] = None
                unknown_flags[key] = unknown_flags.get(key, 0) + 1
        if statuses.get("rnaseq") is True and statuses.get("proteomics") is True:
            etiology_paired[disease] = etiology_paired.get(disease, 0) + 1
            stage_paired[stage] = stage_paired.get(stage, 0) + 1
            flagged.append((tumor, patient, statuses.get("phosphoproteomics")))
    pair_patients = {p for _, p, _ in flagged if p}
    return {
        "status": "PUBLIC_SAMPLE_METADATA_INSPECTED",
        "tumor_rows_with_biopsy_id": tumor_rows,
        "unique_tumor_biopsy_ids": len(tumors),
        "unique_patient_ids": len(patients),
        "rna_and_proteome_yes_biopsy_rows": len(flagged),
        "rna_and_proteome_yes_unique_patient_ids": len(pair_patients),
        "rna_and_proteome_and_phosphoproteome_yes_biopsy_rows":
            sum(1 for _, _, value in flagged if value is True),
        "underlying_liver_disease_raw_counts_all_biopsies": etiology_all,
        "underlying_liver_disease_raw_counts_paired_biopsies": etiology_paired,
        "bclc_raw_counts_paired_biopsies": stage_paired,
        "etiology_warning": "Original categories may be overlapping, nonexclusive and missing. These are biopsy—not independent patient—counts. DO NOT infer two group Ns before patient-level disease codebook validation.",
        "modality_flag_value_counts": hist,
        "uninterpreted_modality_flag_rows": unknown_flags,
        "matrix_join_verified": False,
        "note": "Patient & biopsy identifiers came from a PUBLIC sample sheet and were only counted in memory. No expression/protein matrices joined.",
    }


def inspect_xlsx(data):
    """Public cohort metadata only, no patient/sample identifiers in outputs."""
    if not zipfile.is_zipfile(io.BytesIO(data)):
        raise ValueError("download not an XLSX ZIP workbook")
    with zipfile.ZipFile(io.BytesIO(data)) as book:
        sheets = [n for n in book.namelist() if re.fullmatch(r"xl/worksheets/sheet\d+\.xml", n)]
        if not sheets:
            raise ValueError("no xlsx worksheets")
        def safe_read(name):
            info = book.getinfo(name)
            if info.file_size > MAX_ENTRY:
                raise ValueError("xlsx XML entry too large: " + name)
            return book.read(name)
        shared = []
        if "xl/sharedStrings.xml" in book.namelist():
            root = ET.fromstring(safe_read("xl/sharedStrings.xml"))
            shared = ["".join(n.text or "" for n in si.iter(NS + "t"))
                      for si in root.findall(NS + "si")]
        summary = []
        for name in sheets[:20]:
            root = ET.fromstring(safe_read(name))
            xmlrows = root.findall(".//" + NS + "sheetData/" + NS + "row")
            rows = _as_table_rows(xmlrows, shared)
            colnames = [v for v in rows[0].values() if v] if rows else []
            summary.append({
                "sheet": name, "rows_in_xml": len(rows),
                "header_labels": colnames[:80],
                "modality_sample_metadata": _metadata_pair_summary(rows),
                "expression_matrix_verified": False,
            })
        return {
            "sheets": summary, "worksheet_count": len(sheets),
            "interpretation": "Public sample metadata only. Counts of modality flags are NOT matched RNA/protein expression matrices.",
        }


def probe_ng():
    attempts = []
    for url in NG_URLS:
        entry = {"url": url}
        try:
            payload, code = fetch_bounded(url)
            entry.update({"status": "FETCHED_PUBLIC_XLSX", "http_status": code,
                          "download_bytes": len(payload),
                          "sha256": hashlib.sha256(payload).hexdigest(),
                          "xlsx_schema": inspect_xlsx(payload)})
            attempts.append(entry)
            break
        except (urllib.error.HTTPError, urllib.error.URLError,
                TimeoutError, ValueError, zipfile.BadZipFile, ET.ParseError,
                OSError) as ex:
            entry.update({"status": "UNVERIFIED", "error": str(ex)[:200]})
            attempts.append(entry)
    return {"accession": "Ng2022-SupplementaryData1", "attempts": attempts,
            "pairing_status": "NOT_VERIFIED", "patient_join_n": None}


def probe_gao():
    entry = {"accessions": ["PDC000198", "PDC000199"],
             "endpoint": PDC_ENDPOINT, "metadata_status": "UNVERIFIED",
             "patient_join_n": None, "pairing_status": "NOT_VERIFIED"}
    try:
        payload, code = fetch_bounded(PDC_ENDPOINT, json.dumps({"query": PDC_QUERY}).encode(), limit=500000)
        response = json.loads(payload)
        if response.get("errors"):
            raise ValueError("GraphQL: " + str(response["errors"])[:200])
        data = response.get("data") or {}
        studies = {}
        for key, values in data.items():
            if isinstance(values, dict):
                values = [values]
            if isinstance(values, list):
                studies[key] = {
                    "record_count": len(values),
                    "records": [{k: v.get(k) for k in (
                        "pdc_study_id", "study_name", "cases_count",
                        "aliquots_count", "analytical_fraction")}
                        for v in values[:3] if isinstance(v, dict)]
                }
        if not studies or all(x["record_count"] == 0 for x in studies.values()):
            raise ValueError("PDC returned no matching study records")
        entry.update({"metadata_status": "FETCHED_PUBLIC_PDC_GRAPHQL",
                      "http_status": code, "studies": studies})
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError,
            OSError, ValueError) as ex:
        entry["error"] = str(ex)[:240]
    return entry


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results/hcc/feasibility/pair_sources_probe.json")
    args = parser.parse_args()
    report = {"checked_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
              "source_evidence_only": True, "ng": probe_ng(), "gao": probe_gao(),
              "note": "A downloaded supplement or PDC project record is NOT a validated sample-level RNA/protein intersection."}
    p = pathlib.Path(args.output)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
