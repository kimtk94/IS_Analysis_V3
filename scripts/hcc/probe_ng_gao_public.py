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
    '{ study(pdc_study_id: "PDC000198") { pdc_study_id study_name } '
    'study2: study(pdc_study_id: "PDC000199") { pdc_study_id study_name } }'
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


def inspect_xlsx(data):
    """Inspect only sheet dimensions + first 10 nonempty rows. Never print values from sample rows."""
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
            rows = root.findall(".//" + NS + "sheetData/" + NS + "row")
            preview = []
            for row in rows[:10]:
                cells = row.findall(NS + "c")
                preview.append([_cell_value(c, shared) for c in cells[:50]])
            probable_headers = []
            for row in preview[:4]:
                # Only field names; suppress values that resemble patient or sample IDs.
                keep = [v for v in row if any(x in v.lower() for x in (
                    "patient", "sample", "etiolog", "viral", "alcohol", "hbv", "hcv",
                    "proteom", "transcript", "rna", "seq", "clinical", "stage",
                    "biopsy", "mutation", "tissue", "tumor", "id", "status",
                    "library", "bclc")) and len(v) < 85 and "@" not in v]
                if keep:
                    probable_headers = keep[:30]
                    break
            summary.append({"sheet": name, "rows_in_xml": len(rows),
                            "first_rows_width": [len(v) for v in preview],
                            "candidate_column_labels": probable_headers,
                            "file_level_pair_count": None})
        return {"sheets": summary, "worksheet_count": len(sheets),
                "interpretation": "XLSX sheet discovery only. No RNA/protein sample identity join or patient N verified."}


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
        entry.update({"metadata_status": "FETCHED_PUBLIC_PDC_GRAPHQL",
                      "http_status": code,
                      "studies": {
                          k: {"pdc_study_id": v.get("pdc_study_id"), "study_name": v.get("study_name")}
                          for k, v in data.items() if isinstance(v, dict)
                      }})
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
