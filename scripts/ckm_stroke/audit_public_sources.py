#!/usr/bin/env python3
"""Bounded, fail-closed official-source inventory probe for CKM-stroke GWAS.

Only inspects official landing pages and binary response headers. It NEVER
downloads full GWAS, asserts genome build from URLs, or treats reachability
as genotype/statistical QC.
"""
import argparse
import csv
import datetime
import json
import os
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import Request, urlopen

REQUIRED = {"trait", "ancestry", "source", "kind", "url", "needle", "expected_build"}
ALLOWED_HOSTS = {"www.ebi.ac.uk", "ckdgen.imbi.uni-freiburg.de", "csg.sph.umich.edu"}
KINDS = {"page", "file"}

def read_inventory(path):
    with path.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or not REQUIRED.issubset(reader.fieldnames):
            raise ValueError("source inventory columns missing")
        items = list(reader)
    if not items:
        raise ValueError("empty inventory")
    seen = set()
    for item in items:
        key = (item["trait"], item["ancestry"])
        if key in seen:
            raise ValueError("duplicate source record " + str(key))
        seen.add(key)
        if item["kind"] not in KINDS and item["kind"] != "pending":
            raise ValueError("unknown source kind " + item["kind"])
        if item["kind"] != "pending":
            parts = urlparse(item["url"])
            if parts.scheme != "https" or parts.hostname not in ALLOWED_HOSTS:
                raise ValueError("source not allowlisted: " + item["url"])
    return items

def check_one(item, opener=urlopen, timeout=18, max_bytes=120000):
    if item["kind"] == "pending":
        return {"status": "HOLD", "reason": "official source not yet selected"}
    url = item["url"]
    try:
        if item["kind"] == "page":
            req = Request(url, headers={"User-Agent": "CKM-Stroke-SourceAudit/1.0"})
            with opener(req, timeout=timeout) as resp:
                data = resp.read(max_bytes + 1)
                final = urlparse(resp.geturl())
                if final.hostname not in ALLOWED_HOSTS:
                    return {"status": "HOLD", "reason": "redirect outside trusted source"}
                if len(data) > max_bytes:
                    return {"status": "HOLD", "reason": "landing page exceeds bounded read"}
                text = data.decode("utf-8", errors="replace")
                marker = item["needle"]
                if marker and marker.lower() not in text.lower():
                    return {"status": "HOLD", "reason": "expected accession/trait marker not found"}
                return {"status": "PAGE_VERIFIED", "http_status": resp.status,
                        "bytes_read": len(data), "final_host": final.hostname}
        req = Request(url, method="HEAD", headers={"User-Agent": "CKM-Stroke-SourceAudit/1.0"})
        with opener(req, timeout=timeout) as resp:
            final = urlparse(resp.geturl())
            if final.hostname not in ALLOWED_HOSTS:
                return {"status": "HOLD", "reason": "binary redirected outside trusted source"}
            if not (200 <= resp.status < 300):
                return {"status": "HOLD", "reason": "binary HEAD response not OK", "http_status": resp.status}
            return {"status": "BINARY_REACHABLE", "http_status": resp.status,
                    "content_length": resp.headers.get("Content-Length", ""),
                    "accept_ranges": resp.headers.get("Accept-Ranges", ""),
                    "final_host": final.hostname,
                    "notice": "no full download, hash, decompression or allele/build QC"}
    except (OSError, HTTPError, URLError, TimeoutError, ValueError) as exc:
        return {"status": "SOURCE_UNAVAILABLE", "reason": str(exc)[:280]}

def audit(rows, opener=urlopen):
    details = []
    for item in rows:
        result = check_one(item, opener=opener)
        details.append({k: item[k] for k in ("trait", "ancestry", "source", "kind",
                                            "url", "expected_build")} | result)
    counts = {status: sum(x["status"] == status for x in details) for status
              in ("PAGE_VERIFIED", "BINARY_REACHABLE", "HOLD", "SOURCE_UNAVAILABLE")}
    return {"schema": "ckm-stroke-source-audit-v1",
            "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "source_commit": os.getenv("GITHUB_SHA", ""),
            "items": details, "counts": counts,
            "scientific_status": "HOLD",
            "source_readiness": "HOLD" if counts["HOLD"] or counts["SOURCE_UNAVAILABLE"]
                                else "METADATA_ONLY",
            "notice": "HTTP reachability and page markers never establish summary-statistic correctness"}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory", default="config/ckm_stroke/public_sources.tsv")
    parser.add_argument("--out", default="artifacts/ckm_stroke/public_source_audit.json")
    opts = parser.parse_args()
    rows = read_inventory(Path(opts.inventory))
    report = audit(rows)
    dest = Path(opts.out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"counts": report["counts"], "source_readiness": report["source_readiness"]},
                     ensure_ascii=False))
    for x in report["items"]:
        print(x["trait"], x["ancestry"], x["status"], x.get("reason", ""))
    summary = os.getenv("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a", encoding="utf-8") as fh:
            fh.write("\n### Official-source bounded audit\n\n| Trait | Ancestry | Status |\n|---|---|---|\n")
            for x in report["items"]:
                fh.write(f'| {x["trait"]} | {x["ancestry"]} | {x["status"]} |\n')
            fh.write("\nScientific status: **HOLD** (no primary GWAS inference)\n")
    return 0  # report HOLD as evidence, not a false-negative CI failure

if __name__ == "__main__":
    raise SystemExit(main())
