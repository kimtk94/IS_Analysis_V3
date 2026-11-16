#!/usr/bin/env python3
"""Bounded P1 GWAS source landing-page audit (metadata only, never promotes sources).

The output is an evidence triage, not complete summary-stat, LDSC or MR input QC.
"""
import argparse
import csv
import datetime
import json
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request, urlopen

ALLOWED = {
    "giant-consortium.web.broadinstitute.org",
    "gwas.mrcieu.ac.uk",
    "opengwas.io",
    "www.diagram-consortium.org",
}
COLUMNS = {
    "trait", "priority", "study", "ancestry", "landing_url",
    "metadata_source_url", "reported_build", "reported_sample_size",
    "published_md5", "page_marker", "access_and_qc_notes", "status",
}
MAX_READ = 256000


def manifest(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not COLUMNS.issubset(set(reader.fieldnames or [])):
            raise ValueError("candidate inventory missing columns")
        rows = list(reader)
    if not rows or set(r["trait"] for r in rows) != {"BMI", "SBP", "T2D"}:
        raise ValueError("inventory must contain BMI, SBP, and T2D")
    seen = set()
    for row in rows:
        key = (row["trait"], row["priority"])
        if key in seen or row["status"] != "CANDIDATE_ONLY" or row["ancestry"] != "EUR":
            raise ValueError("duplicate, approved, or non-EUR source")
        seen.add(key)
        if row["priority"] not in ("primary", "secondary"):
            raise ValueError("unexpected priority")
        for url_col in ("landing_url", "metadata_source_url"):
            url = urlparse(row[url_col])
            if url.scheme != "https" or not url.netloc:
                raise ValueError("candidate URL must be HTTPS")
        if urlparse(row["landing_url"]).hostname not in ALLOWED:
            raise ValueError("nonapproved landing host")
        if not row["page_marker"]:
            raise ValueError("empty page marker")
    return rows


def probe(row, opener=urlopen, byte_cap=MAX_READ, timeout=12):
    url = row["landing_url"]
    req = Request(url, headers={"User-Agent": "CKM-IS-P1-CandidateAudit/1.0"})
    try:
        with opener(req, timeout=timeout) as response:
            final = urlparse(response.geturl())
            if final.scheme != "https" or final.hostname not in ALLOWED:
                return {"status": "HOLD_REDIRECT", "final_host": final.hostname}
            data = response.read(byte_cap + 1)
            if len(data) > byte_cap:
                return {"status": "HOLD_PAGE_TOO_LARGE", "bytes_read": len(data)}
            txt = " ".join(data.decode("utf-8", errors="replace").split()).lower()
            expected = " ".join(row["page_marker"].split()).lower()
            matched = expected in txt
            return {
                "status": "LANDING_MARKER_FOUND" if matched else "HOLD_MARKER_NOT_FOUND",
                "http_status": getattr(response, "status", None),
                "final_host": final.hostname,
                "bytes_read": len(data),
                "marker_found": matched,
                "ready_for_ldsc": False,
                "original_binary_downloaded": False,
            }
    except Exception as exc:
        return {"status": "SOURCE_UNAVAILABLE", "error_type": type(exc).__name__,
                "detail": str(exc)[:180], "ready_for_ldsc": False}


def audit(rows, opener=urlopen):
    observations = []
    for row in rows:
        result = probe(row, opener=opener)
        observations.append({
            "trait": row["trait"], "priority": row["priority"],
            "study": row["study"], "landing_url": row["landing_url"],
            "inventory_status": row["status"], **result,
        })
    return {
        "schema": "ckm-stroke-p1-candidate-landing-audit-v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "research_status": "HOLD",
        "source_approved": False,
        "ready_for_ldsc": False,
        "observations": observations,
        "limitations": "Landing markers are not original file headers, checksums, ancestry/LD match, or verified sample overlap.",
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--inventory", default="config/ckm_stroke/p1_gwas_candidates.tsv")
    p.add_argument("--out", required=True)
    args = p.parse_args()
    result = audit(manifest(args.inventory))
    dest = Path(args.out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    for row in result["observations"]:
        print(row["trait"], row["priority"], row["status"], row.get("detail", ""))
    print("SCIENTIFIC_HOLD; no source promoted")


if __name__ == "__main__":
    main()
