#!/usr/bin/env python3
"""Bounded PUBLIC manifest audit for HCC quantitative-matrix feasibility.

Only file names, sizes, and workbook *shapes* are emitted. Never fetch restricted
EGA RNA, PRIDE raw mass spectra, individual case IDs, or hospital information.
A returned file name is NOT proof its table can be joined to a transcriptome.
"""
import argparse
import collections
import datetime as dt
import io
import json
import pathlib
import re
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

PDC = "https://pdc.cancer.gov/graphql"
PRIDE = "https://www.ebi.ac.uk/pride/ws/archive/v3"
SPRINGER = ("https://media.springernature.com/original/springer-static/esm/"
            "art%3A10.1038%2Fs41467-022-29960-8/MediaObjects/")
NG_FILES = {
    "Supplementary_Data_2": "41467_2022_29960_MOESM5_ESM.xlsx",
    "Supplementary_Data_13": "41467_2022_29960_MOESM16_ESM.xlsx",
    "Supplementary_Data_14": "41467_2022_29960_MOESM17_ESM.xlsx",
    "Source_Data": "41467_2022_29960_MOESM19_ESM.xlsx",
}
PDC_IDS = ("PDC000198", "PDC000199")
PRIDE_IDS = ("PXD025705", "PXD025836")
ZIP_LIMIT = 7 * 1024 * 1024
XML_LIMIT = 24 * 1024 * 1024
JSON_LIMIT = 2 * 1024 * 1024
N = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def read_bounded(url, *, limit, data=None, timeout=17):
    headers = {"User-Agent": "HCC-Feasibility-Public-Audit/2.1", "Accept": "application/json, */*"}
    if data is not None:
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, headers=headers,
                                 method="POST" if data is not None else "GET")
    with urllib.request.urlopen(req, timeout=timeout) as response:
        payload = response.read(limit + 1)
    if len(payload) > limit:
        raise ValueError("response exceeded audited byte cap")
    return payload


def request_graphql(query):
    payload = read_bounded(PDC, limit=JSON_LIMIT,
                           data=json.dumps({"query": query}).encode("utf-8"))
    answer = json.loads(payload)
    if answer.get("errors"):
        raise ValueError("GraphQL validation: " + str(answer["errors"])[:300])
    if not isinstance(answer.get("data"), dict):
        raise ValueError("GraphQL has no result data")
    return answer["data"]


def probe_pdc():
    report = {}
    for acc in PDC_IDS:
        item = {"accession": acc, "source": PDC, "status": "UNVERIFIED"}
        try:
            q = ('{ study(pdc_study_id:"' + acc +
                 '"){study_id pdc_study_id study_name cases_count aliquots_count '
                 'filesCount{data_category file_type files_count}} }')
            result = request_graphql(q).get("study")
            study = result[0] if isinstance(result, list) and result else (
                result if isinstance(result, dict) else None)
            if not study:
                raise ValueError("study not found")
            item["study_id"] = study.get("study_id")
            item["repository_cases_count"] = study.get("cases_count")
            item["repository_aliquots_count"] = study.get("aliquots_count")
            item["file_type_counts"] = study.get("filesCount", [])[:100]
            item["status"] = "STUDY_FILETYPE_MANIFEST_FETCHED"
            # Limit file records; do not retrieve spectra or signed download URLs.
            q2 = ('{getPaginatedFiles(pdc_study_id:"' + acc +
                  '" offset:0 limit:75){total files{file_name file_type '
                  'file_format data_category} }}')
            records = request_graphql(q2).get("getPaginatedFiles") or {}
            if isinstance(records, list):
                records = records[0] if records else {}
            files = records.get("files") or []
            item["files_total_reported"] = records.get("total")
            item["files_preview_count"] = len(files)
            item["files_preview"] = [
                {k: x.get(k) for k in ("file_name", "file_type", "file_format", "data_category")}
                for x in files[:75] if isinstance(x, dict)
            ]
            item["status"] = "PAGED_FILE_MANIFEST_FETCHED"
        except Exception as exc:
            item["error"] = str(exc)[:250]
        report[acc] = item
    return report


def probe_pride():
    output = {}
    for acc in PRIDE_IDS:
        item = {"accession": acc, "status": "UNVERIFIED", "source": PRIDE}
        try:
            url = f"{PRIDE}/projects/{acc}/files?pageSize=100&page=0"
            payload = read_bounded(url, limit=JSON_LIMIT)
            data = json.loads(payload)
            if isinstance(data, dict):
                files = (data.get("_embedded") or {}).get("files") or data.get("files") or data.get("content") or []
            else:
                files = data
            if not isinstance(files, list):
                raise ValueError("unexpected PRIDE file list response")
            item["status"] = "FIRST_PAGE_FILE_MANIFEST_FETCHED"
            item["file_rows_first_page"] = len(files)
            counts = collections.Counter()
            candidates = []
            for f in files:
                if not isinstance(f, dict):
                    continue
                kind = f.get("fileCategory")
                label = kind.get("value", str(kind)) if isinstance(kind, dict) else str(kind)
                counts[label] += 1
                name = f.get("fileName") or f.get("fileNamePath") or f.get("name") or ""
                size = f.get("fileSizeBytes")
                if re.search(r"protein|quant|matrix|expression|normaliz|swath|sdrf|result|report|evidence|peptide", name, re.I):
                    candidates.append({"name": name, "bytes": size, "category": label})
            item["category_count_first_page"] = dict(counts)
            item["candidate_files_first_page"] = candidates[:60]
            item["note"] = "First 100 records only. Missing file not proof of absence."
        except Exception as exc:
            item["error"] = str(exc)[:250]
        output[acc] = item
    return output


def workbook_shape(payload):
    if not zipfile.is_zipfile(io.BytesIO(payload)):
        raise ValueError("not an XLSX archive")
    results = []
    with zipfile.ZipFile(io.BytesIO(payload)) as zf:
        paths = sorted(n for n in zf.namelist() if re.fullmatch(r"xl/worksheets/sheet\d+\.xml", n))
        for path in paths[:25]:
            info = zf.getinfo(path)
            if info.file_size > XML_LIMIT:
                results.append({"sheet_xml": path, "uncompressed_xml_bytes": info.file_size,
                                "status": "SKIPPED_XML_TOO_LARGE"})
                continue
            root = ET.fromstring(zf.read(path))
            dimension = root.find(N + "dimension")
            row_nodes = root.findall(".//" + N + "sheetData/" + N + "row")
            width_first = []
            for row in row_nodes[:2]:
                width_first.append(len(row.findall(N + "c")))
            results.append({"sheet_xml": path, "shape_ref": dimension.attrib.get("ref") if dimension is not None else None,
                            "xml_row_count": len(row_nodes), "first_row_cell_counts": width_first,
                            "uncompressed_xml_bytes": info.file_size})
    return results


def probe_ng_workbooks():
    reports = {}
    for name, filename in NG_FILES.items():
        url = SPRINGER + filename
        item = {"name": name, "file": filename, "source": url, "status": "UNVERIFIED"}
        try:
            data = read_bounded(url, limit=ZIP_LIMIT)
            item["bytes"] = len(data)
            item["sheets"] = workbook_shape(data)
            item["status"] = "PUBLIC_XLSX_SHAPE_VERIFIED"
            item["quantitative_matrix_join_status"] = "NOT_VERIFIED"
        except Exception as exc:
            item["error"] = str(exc)[:250]
        reports[name] = item
    return reports


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results/hcc/feasibility/quant_matrix_manifest.json")
    args = parser.parse_args()
    report = {
        "utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "audit_scope": "Public metadata/file list and workbook geometry; never expression values or patient IDs.",
        "ng_supplements": probe_ng_workbooks(),
        "gao_pdc": probe_pdc(),
        "ng_pride": probe_pride(),
        "RNA_PROTEIN_MATRIX_JOIN": "NOT_VERIFIED",
    }
    output = pathlib.Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    # Concise run logs; full public file-name listing goes into Actions artifact.
    concise = {"ng": {k: {"status": v["status"], "sheets": v.get("sheets", []), "error": v.get("error")}
                       for k, v in report["ng_supplements"].items()},
               "pdc": {k: {"status": v["status"], "files_total": v.get("files_total_reported"),
                           "file_types": v.get("file_type_counts", []), "error": v.get("error")}
                       for k, v in report["gao_pdc"].items()},
               "pride": {k: {"status": v["status"], "page_files": v.get("file_rows_first_page"),
                             "candidate_names": [f["name"] for f in v.get("candidate_files_first_page", [])[:15]],
                             "error": v.get("error")}
                         for k, v in report["ng_pride"].items()},
               "matrix_join": "NOT_VERIFIED"}
    print(json.dumps(concise, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
