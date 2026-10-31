#!/usr/bin/env python3
"""PDC Gao Stage33: inspect CDAP .sample.txt and measure SOURCE-derived mapping.

No individual sample IDs, headers, matrix values, filenames or signed URLs are
persisted, logged or transferred. Does not certify joins from row position.
Document: https://pdc.cancer.gov/data-dictionary/CDAP_ProteinReports_description_20180914.pdf
"""
import argparse
import collections
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import pathlib
import re
import urllib.request

from probe_stage30_gao_catalog import PDC, pdc_query
from probe_stage30_gao_quant_header import safe_public_https
from probe_stage30_gao_quant_api_header_to_biospecimen import get_quant_header, official_biospecimen

MAX_SOURCE = 2 * 1024 * 1024
MAX_ROWS = 1500
RATIO_PATTERN = re.compile(r"^\s*(?:tmt[_\-]?)?(\d{3}[nc]?)\s*/\s*(?:tmt[_\-]?)?(\d{3}[nc]?|1\.0)\s*$", re.I)
REPORTER = re.compile(r"^(?:TMT[_\-]?)?1\d{2}[NC]?$", re.I)

def discover(pid, graphql=pdc_query):
    q = ('{filesPerStudy(pdc_study_id:' + json.dumps(pid) +
         ',data_category:"Protein Assembly",offset:0,limit:30)'
         '{file_name data_category file_format file_size md5sum signedUrl{url}}}')
    data, sha = graphql(q)
    rows = data.get("filesPerStudy")
    if not isinstance(rows, list) or len(rows) > 30:
        raise ValueError("PDC source file schema drift")
    if any(not isinstance(r, dict) or r.get("data_category") != "Protein Assembly" for r in rows):
        raise ValueError("PDC category not honored")
    found = []
    for r in rows:
        if not re.search(r"\.sample\.(?:txt|tsv)(?:\.gz)?$", str(r.get("file_name") or ""), re.I):
            continue
        urls = r.get("signedUrl") or []
        if isinstance(urls, dict):
            urls = [urls]
        url = next((x.get("url") for x in urls if isinstance(x, dict)
                    and safe_public_https(x.get("url"))), None)
        found.append({"url": url, "size": r.get("file_size"),
                      "md5": r.get("md5sum"), "gz": str(r.get("file_name")).lower().endswith(".gz")})
    return found, sha

def read_bytes(item, opener=urllib.request.urlopen):
    u = item["url"]
    if not safe_public_https(u):
        raise ValueError("sample companion signed URL not allowlisted")
    request = urllib.request.Request(u, headers={
        "Range": "bytes=0-2097151",
        "User-Agent": "HCC-Stage33-PDC-CDAP-sample-companion/1"})
    with opener(request, timeout=30) as response:
        body = response.read(MAX_SOURCE + 1)
    if len(body) > MAX_SOURCE:
        raise ValueError("source is larger than bounded read")
    # Incomplete Range responses cannot be asserted to be a complete sample mapping.
    declared = item.get("size")
    complete = isinstance(declared, (str, int)) and str(declared).isdigit() and int(declared) == len(body)
    if not complete:
        raise ValueError("source is partial or file_size unavailable")
    recorded_md5 = item.get("md5")
    md5_ok = bool(re.fullmatch(r"[a-fA-F0-9]{32}", str(recorded_md5 or ""))) and (
        hashlib.md5(body).hexdigest().lower() == str(recorded_md5).lower())
    if not md5_ok:
        raise ValueError("download bytes fail source MD5")
    if item["gz"]:
        with gzip.GzipFile(fileobj=io.BytesIO(body)) as compressed:
            body = compressed.read(MAX_SOURCE + 1)
        if len(body) > MAX_SOURCE:
            raise ValueError("decompressed sample file exceeds limit")
    return body

def parse_sample_table(raw):
    txt = raw.decode("utf-8-sig")
    lines = [x for x in txt.splitlines() if x.strip() and not x.lstrip().startswith("#")]
    if not lines:
        raise ValueError("source .sample.txt empty")
    table = list(csv.reader(lines, delimiter="\t"))
    if len(table) < 2 or len(table) > MAX_ROWS + 1:
        raise ValueError("unexpected source experiment design row count")
    hdr = [x.strip() for x in table[0]]
    if len(hdr) != len(set(hdr)):
        raise ValueError("duplicate source headings")
    if not {"AnalyticalSample", "Ratios"}.issubset(hdr):
        raise ValueError("source is not CDAP sample.txt experimental design")
    if any(len(x) != len(hdr) for x in table[1:]):
        raise ValueError("nonrectangular source companion")
    reporters = [h for h in hdr if REPORTER.fullmatch(h)]
    if not reporters:
        raise ValueError("no source TMT reporter columns")
    rows = [dict(zip(hdr, x)) for x in table[1:]]
    return rows, reporters

def report_candidates(label):
    # Only source-derived CDAP report conventions (not run/plex guesses).
    v = label.strip()
    if not v or v.upper() == "POOL":
        return set()
    return {v, v + " Log Ratio", v + " Unshared Log Ratio"}

def compare(rows, reporters, quant_headers, biospecimens):
    qset = set(quant_headers)
    all_known = {}
    for ns in ("aliquot_id", "aliquot_submitter_id"):
        for alias, records in biospecimens.get(ns, {}).items():
            all_known.setdefault(alias, set()).update(r["aliquot"] for r in records)
    source_labels = []
    numerator_labels = []
    parsed = 0
    unknown_ratios = 0
    reference_labels = 0
    ratios_per_row = []
    for row in rows:
        for reporter in reporters:
            val = (row.get(reporter) or "").strip()
            if val:
                source_labels.append(val)
                reference_labels += (val.upper() == "POOL")
        specs = [x.strip() for x in (row.get("Ratios") or "").split(",") if x.strip()]
        ratios_per_row.append(len(specs))
        for spec in specs:
            m = RATIO_PATTERN.fullmatch(spec)
            if not m:
                unknown_ratios += 1
                continue
            num = m.group(1).upper()
            candidates = [h for h in reporters if h.upper().replace("TMT_", "").replace("TMT-", "") == num]
            if len(candidates) != 1:
                unknown_ratios += 1
                continue
            label = (row.get(candidates[0]) or "").strip()
            if label:
                numerator_labels.append(label)
                parsed += 1
    # Compare source sample-ID columns to true official IDs; avoid aliases
    # appearing in unrelated clinical fields or arbitrary first-row order.
    official = [all_known.get(label, set()) for label in numerator_labels]
    unique_known = [next(iter(m)) for m in official if len(m) == 1]
    label_to_aliquots = collections.defaultdict(set)
    for label, ids in zip(numerator_labels, official):
        label_to_aliquots[label].update(ids)
    ambiguous_ids = sum(len(x) > 1 for x in official)
    alias_counts = collections.Counter()
    matched_by_header = collections.defaultdict(set)
    for label, ids in zip(numerator_labels, official):
        for cand in report_candidates(label):
            if cand in qset:
                alias_counts[cand] += 1
                matched_by_header[cand].update(ids)
    unique_headers_with_one_aliquot = {
        k: next(iter(v)) for k,v in matched_by_header.items()
        if len(v) == 1 and alias_counts[k] == 1
    }
    return {
        "source_analytical_sample_rows_n": len(rows),
        "source_reporter_column_n": len(reporters),
        "source_nonblank_reporter_assignments_n": len(source_labels),
        "source_reporter_reference_literal_pool_n": reference_labels,
        "source_ratio_specifications_n": sum(ratios_per_row),
        "source_ratios_with_parseable_nonempty_numerator_n": parsed,
        "source_ratios_unparsed_or_unassigned_n": unknown_ratios,
        "source_ratio_labels_distinct_n": len(set(numerator_labels)),
        "source_ratio_labels_repeated_n": len(numerator_labels) - len(set(numerator_labels)),
        "source_ratio_labels_matching_unique_official_aliquot_n": len(unique_known),
        "source_ratio_labels_with_ambiguous_official_aliquot_n": ambiguous_ids,
        "source_ratio_distinct_unique_official_aliquots_n": len(set(unique_known)),
        "source_ratio_headers_with_direct_CDAP_label_match_n": len(matched_by_header),
        "source_ratio_headers_unique_official_aliquot_n": len(unique_headers_with_one_aliquot),
        "source_ratio_headers_distinct_unique_aliquots_n": len(set(unique_headers_with_one_aliquot.values())),
        "quant_columns_n": len(quant_headers),
        "candidate_full_mapping": len(quant_headers) == len(unique_headers_with_one_aliquot) == 330
            and len(set(unique_headers_with_one_aliquot.values())) == 330,
        "mapping_proven": False,  # source alias shape alone cannot certify a link
        "per_sample_or_header_values_exported": False,
    }

def one(pid):
    out = {"pdc_study_id": pid, "status": "ACCESS_HOLD",
           "sample_ids_or_matrix_header_values_exported": False,
           "signed_urls_or_file_names_exported": False,
           "source_to_biospecimen_join_certified": False,
           "RNA_to_protein_patient_join_certified": False}
    try:
        found, digest = discover(pid)
        out["source_metadata_sha256"] = digest
        out["cdap_sample_file_candidates_n"] = len(found)
        if len(found) != 1:
            out["status"] = "REVIEW_AMBIGUOUS_OR_MISSING_SAMPLE_COMPANION"
            return out
        raw = read_bytes(found[0])
        out["source_sample_file_sha256"] = hashlib.sha256(raw).hexdigest()
        rows, reporters = parse_sample_table(raw)
        quant_meta, qh = get_quant_header(pid)
        if quant_meta.get("status") != "VERIFIED_PDC_QUANT_JSON_HEADER":
            raise ValueError("unverified quant column source")
        biospecimens = official_biospecimen(pid)
        out.update(compare(rows, reporters, qh, biospecimens))
        out["status"] = "SOURCE_CDAP_SAMPLE_TABLE_CROSSCHECKED"
    except Exception as exc:
        out["error_class"] = type(exc).__name__
    return out

def review(payload):
    if payload.get("schema") != "hcc-stage33-CDAP-original-sample-bridge-v1":
        raise ValueError("Stage33 schema drift")
    evidence = payload.get("studies")
    if not isinstance(evidence, list) or {x.get("pdc_study_id") for x in evidence} != set(PDC):
        raise ValueError("study identity drift")
    for x in evidence:
        for flag in ("sample_ids_or_matrix_header_values_exported",
                     "signed_urls_or_file_names_exported",
                     "source_to_biospecimen_join_certified",
                     "RNA_to_protein_patient_join_certified"):
            if x.get(flag) is not False:
                raise ValueError("unsupported patient claim or unsafe export")
        if x.get("status") == "SOURCE_CDAP_SAMPLE_TABLE_CROSSCHECKED":
            if x.get("quant_columns_n") != 330:
                raise ValueError("quant matrix width drift")
            if x.get("mapping_proven") is not False or x.get("per_sample_or_header_values_exported") is not False:
                raise ValueError("premature mapping declaration")
        elif x.get("status") not in ("ACCESS_HOLD", "REVIEW_AMBIGUOUS_OR_MISSING_SAMPLE_COMPANION"):
            raise ValueError("invalid stage33 status")
    return {
        "schema": "hcc-stage33-CDAP-original-sample-bridge-review-v1",
        "CDAP_source_structural_comparisons_completed_n":
            sum(x["status"] == "SOURCE_CDAP_SAMPLE_TABLE_CROSSCHECKED" for x in evidence),
        "G30_B": "HOLD_UNTIL_330_SOURCE_COLUMN_ALIQUOT_LINKS_CERTIFIED",
        "G30_C": "HOLD_NODE_RNA_MATRIX",
        "G30_D": "STOP_NO_MATCHED_EXTERNAL_RNA_PROTEIN_CASES",
        "biological_inference_allowed": False,
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--review-out", required=True)
    a = ap.parse_args()
    payload = {"schema": "hcc-stage33-CDAP-original-sample-bridge-v1",
               "checked_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
               "studies": [one(pid) for pid in PDC]}
    verdict = review(payload)
    for path, x in ((a.out, payload), (a.review_out, verdict)):
        p = pathlib.Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(x, indent=2) + "\n")
    print("HCC_STAGE33_AGGREGATES", json.dumps(payload))
    print("HCC_STAGE33_REVIEW", json.dumps(verdict))

if __name__ == "__main__":
    main()
