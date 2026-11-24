#!/usr/bin/env python3
"""Direct PDC study-ID file category discovery, independent of Stage31 UUID lookup.

Read-only metadata. Save anonymous counts/digests only. A GCT filename or
technical API success is never proof of an aliquot or patient-level join.
"""
import argparse
import collections
import datetime as dt
import json
import pathlib
import re
from probe_stage30_gao_catalog import PDC, pdc_query

CATEGORIES = ("Other Metadata", "Protein Assembly")
MAX_FILES = 30

def classify(name):
    s = str(name or "").lower()
    return {
        "gct_suffix": bool(re.search(r"\.gct(?:\.gz)?$", s)),
        "cdap_sample_suffix": bool(re.search(r"\.sample\.(?:txt|tsv)(?:\.gz)?$", s)),
        "sample_design_keyword": bool(re.search(r"sample.?map|sample.?sheet|experimental.?design|plex.?map", s)),
        "quant_matrix_keyword": bool(re.search(r"quant|abundan|matrix|protein|phospho", s)),
    }

def query(pid, category, fetch=pdc_query):
    if pid not in PDC or category not in CATEGORIES:
        raise ValueError("unapproved source query")
    result = {
        "pdc_study_id": pid, "category": category,
        "status": "ACCESS_HOLD", "source_identifiers_or_names_exported": False,
        "signed_urls_exported": False, "sample_to_aliquot_join_verified": False,
        "patient_level_multiomics_verified": False,
    }
    # Do not first call study(pdc_study_id): one failed UUID lookup must not
    # obscure whether filesPerStudy can list metadata for the accession.
    q = ('{filesPerStudy(pdc_study_id:' + json.dumps(pid) +
         ',data_category:' + json.dumps(category) +
         ',offset:0,limit:' + str(MAX_FILES) +
         '){file_name data_category file_type file_format md5sum}}')
    try:
        data, digest = fetch(q)
        rows = data.get("filesPerStudy")
        if not isinstance(rows, list) or len(rows) > MAX_FILES:
            raise ValueError("unexpected filesPerStudy schema")
        if any(not isinstance(x, dict) or x.get("data_category") != category for x in rows):
            raise ValueError("PDC category filter ignored")
        flags, formats = collections.Counter(), collections.Counter()
        recorded_md5s = 0
        for row in rows:
            flags.update({k: int(v) for k, v in classify(row.get("file_name")).items()})
            formats[str(row.get("file_format") or "unknown")[:60]] += 1
            recorded_md5s += bool(re.fullmatch(r"[0-9a-fA-F]{32}", str(row.get("md5sum") or "")))
        result.update({
            "status": "OFFICIAL_CATEGORY_PAGE_VERIFIED",
            "first_page_files_n": len(rows),
            "page_limit_n": MAX_FILES,
            "gct_filename_candidates_n": flags["gct_suffix"],
            "cdap_sample_filename_candidates_n": flags["cdap_sample_suffix"],
            "sample_design_filename_candidates_n": flags["sample_design_keyword"],
            "quant_matrix_filename_candidates_n": flags["quant_matrix_keyword"],
            "format_counts": dict(sorted(formats.items())),
            "source_response_sha256": digest,
            "valid_md5_metadata_n": recorded_md5s,
            "category_exhaustive": False,  # no count or last-page proof
        })
    except Exception as exc:
        # Never print raw provider errors: they may expose URLs or identifiers.
        result["error_class"] = type(exc).__name__
    return result

def verify(obj):
    if obj.get("schema") != "hcc-stage32-pdc-direct-file-category-v1":
        raise ValueError("schema drift")
    x = obj.get("evidence")
    if not isinstance(x, dict) or x.get("pdc_study_id") not in PDC or x.get("category") not in CATEGORIES:
        raise ValueError("source drift")
    for name in ("source_identifiers_or_names_exported", "signed_urls_exported",
                 "sample_to_aliquot_join_verified", "patient_level_multiomics_verified"):
        if x.get(name) is not False:
            raise ValueError("prohibited identifier exposure or inference")
    if x.get("status") == "OFFICIAL_CATEGORY_PAGE_VERIFIED":
        n = x.get("first_page_files_n")
        if type(n) is not int or not 0 <= n <= MAX_FILES:
            raise ValueError("impossible page")
        for k in ("gct_filename_candidates_n", "cdap_sample_filename_candidates_n",
                  "sample_design_filename_candidates_n", "quant_matrix_filename_candidates_n",
                  "valid_md5_metadata_n"):
            if type(x.get(k)) is not int or not 0 <= x[k] <= n:
                raise ValueError("invalid candidate count")
        if x.get("category_exhaustive") is not False:
            raise ValueError("unproved exhaustive inventory")
    elif x.get("status") != "ACCESS_HOLD":
        raise ValueError("unsupported status")
    return {
        "schema": "hcc-stage32-pdc-file-category-postflight-v1",
        "source_category_page_technical_success": x["status"] == "OFFICIAL_CATEGORY_PAGE_VERIFIED",
        "PDC_330_quant_header_to_aliquot_mapping": "HOLD_UNVERIFIED",
        "NODE_RNA_matrix_access": "HOLD_UNVERIFIED",
        "external_paired_RNA_protein_analysis": "STOP",
        "new_biological_claims_allowed": False,
    }

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--study", choices=PDC, required=True)
    p.add_argument("--category", choices=CATEGORIES, required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--review-out", required=True)
    a = p.parse_args()
    obj = {"schema": "hcc-stage32-pdc-direct-file-category-v1",
           "checked_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "evidence": query(a.study, a.category)}
    post = verify(obj)
    for name, content in ((a.out, obj), (a.review_out, post)):
        path = pathlib.Path(name)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(content, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("HCC_STAGE32_AGGREGATES", json.dumps(obj, sort_keys=True))
    print("HCC_STAGE32_POSTFLIGHT", json.dumps(post, sort_keys=True))

if __name__ == "__main__":
    main()
