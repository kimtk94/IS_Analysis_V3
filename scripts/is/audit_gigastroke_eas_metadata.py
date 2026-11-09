#!/usr/bin/env python3
"""Read-only GIGASTROKE EAS GWAS Catalog v2 metadata audit.

Does NOT inspect full summary-statistics files, alleles or QC, and does NOT
measure biological replication. Fails closed when registry unavailable.
"""
import argparse
import json
import os
import pathlib
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone

BASE = "https://www.ebi.ac.uk/gwas/rest/api/v2/studies"
CONTROL_N = 237242
EXPECTED = [
    {"id": "GCST90104544", "outcome": "any stroke", "cases": 27413},
    {"id": "GCST90104545", "outcome": "ischemic stroke", "cases": 19032},
    {"id": "GCST90104546", "outcome": "cardioembolic stroke", "cases": 926},
    {"id": "GCST90104547", "outcome": "large artery stroke", "cases": 1735},
    {"id": "GCST90104548", "outcome": "small vessel stroke", "cases": 5532},
]


def validate_constants():
    assert len({x["id"] for x in EXPECTED}) == 5
    totals = [264655, 256274, 238168, 238977, 242774]
    assert [x["cases"] + CONTROL_N for x in EXPECTED] == totals
    return True


def find_trait(node):
    if isinstance(node, dict):
        for key in ("disease_trait", "diseaseTrait", "reported_trait", "reportedTrait"):
            value = node.get(key)
            if isinstance(value, str) and value.strip():
                return value.strip()
        for val in node.values():
            found = find_trait(val)
            if found:
                return found
    if isinstance(node, list):
        for val in node:
            found = find_trait(val)
            if found:
                return found
    return None


def fetch_json(url):
    req = urllib.request.Request(
        url,
        headers={"Accept": "application/json", "User-Agent": "IS-Analysis-V3-ResearchMetadataAudit/1.0"},
    )
    last_error = None
    for attempt in range(2):
        try:
            with urllib.request.urlopen(req, timeout=20) as resp:
                if resp.status != 200:
                    raise RuntimeError(f"HTTP {resp.status}")
                return json.load(resp), None
        except (urllib.error.URLError, TimeoutError, ValueError, RuntimeError) as ex:
            last_error = f"{type(ex).__name__}: {str(ex)[:200]}"
            if attempt == 0:
                time.sleep(1)
    return None, last_error


def audit_one(item):
    accession = item["id"]
    url = f"{BASE}/{accession}"
    result = dict(item, controls=CONTROL_N, metadata_url=url, status="SOURCE_UNAVAILABLE")
    payload, error = fetch_json(url)
    if payload is None:
        result["error"] = error
        return result
    raw = json.dumps(payload, ensure_ascii=False).lower()
    trait = find_trait(payload)
    result["catalog_trait"] = trait
    if accession.lower() not in raw or "36180795" not in raw:
        result["status"] = "METADATA_MISMATCH"
        result["error"] = "Expected accession and PMID 36180795 not present in response"
    elif not trait:
        result["status"] = "METADATA_INCOMPLETE"
        result["error"] = "No reported trait field found: cannot verify phenotype"
    elif item["outcome"] == "any stroke":
        result["status"] = "METADATA_MATCH" if "stroke" in trait.lower() and "ischemic" not in trait.lower() else "REVIEW_TRAIT"
    elif item["outcome"].split()[0] in trait.lower():
        result["status"] = "METADATA_MATCH"
    else:
        result["status"] = "REVIEW_TRAIT"
        result["error"] = f"Expected {item['outcome']}; received {trait}"
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="artifacts/is_eas_gwas_catalog")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    validate_constants()
    if args.self_test:
        print("SELF_TEST_PASS: reference IDs and literature case-count arithmetic")
        return 0
    out = pathlib.Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows = [audit_one(item) for item in EXPECTED]
    audit = {
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": "GWAS Catalog metadata REST v2",
        "source_note": "Metadata only. No proof of availability or validity of full summary statistics.",
        "sample_count_provenance": "Published literature (not API-verified): doi:10.3390/jcdd11080247",
        "rows": rows,
    }
    (out / "metadata_audit.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    headers = ["Study", "Outcome", "Published cases", "Source metadata", "Status"]
    lines = ["| " + " | ".join(headers) + " |", "|---|---|---:|---|---|"]
    for x in rows:
        lines.append(
            f"| {x['id']} | {x['outcome']} | {x['cases']:,} | "
            f"{x.get('catalog_trait') or '-'} | {x['status']} |"
        )
    lines.extend(["", "Metadata checks only; study counts are literature-sourced, not verified from downloaded GWAS files."])
    summary = "\n".join(lines) + "\n"
    (out / "summary.md").write_text(summary, encoding="utf-8")
    print(summary)
    step_summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if step_summary:
        with open(step_summary, "a", encoding="utf-8") as handle:
            handle.write("## GIGASTROKE EAS metadata audit (not an omics analysis)\n" + summary)
    statuses = [x["status"] for x in rows]
    if not all(x == "METADATA_MATCH" for x in statuses):
        print("HOLD: not all catalog metadata matches verified", file=sys.stderr)
        return 2
    print("METADATA_AUDIT_PASS: only metadata. G0 full-file gate still NOT VERIFIED.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
