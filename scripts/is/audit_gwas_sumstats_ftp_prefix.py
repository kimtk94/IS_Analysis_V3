#!/usr/bin/env python3
"""Read-only, bounded discovery of original GIGASTROKE GWAS Catalog FTP files.

This script identifies a remote file and samples a small prefix for schema QC.
It does NOT download/verify the entire file, reproduce source checksums, check
genome-wide statistics, establish cohort diagnostic adjudication, or perform
MAGMA/LDSC. Every output explicitly marks those steps as HOLD.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
from html.parser import HTMLParser
import json
import os
from pathlib import Path
import re
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
import zlib

ROOT = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
EXPECTED = {
    "EUR": {"accession": "GCST90104540", "phenotype": "ischemic_stroke"},
    "EAS": {"accession": "GCST90104545", "phenotype": "ischemic_stroke"},
}
MAX_HTML = 512 * 1024
MAX_PREFIX = 512 * 1024
MAX_EXPANDED = 512 * 1024
TIMEOUT = 25


def accession_directory(accession):
    match = re.fullmatch(r"GCST(\d{6,9})", accession)
    if match is None:
        raise ValueError("Invalid GWAS Catalog accession")
    n = int(match.group(1))
    lo = ((n - 1) // 1000) * 1000 + 1
    hi = lo + 999
    width = len(match.group(1))
    return (
        f"{ROOT}/GCST{lo:0{width}d}-GCST{hi:0{width}d}/"
        f"{accession}/"
    )


class Listing(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        if tag.lower() == "a":
            address = dict(attrs).get("href", "")
            if address:
                self.links.append(address)


def extract_directory_names(html, base):
    parser = Listing()
    parser.feed(html)
    names = set()
    for href in parser.links:
        absolute = urllib.parse.urljoin(base, href)
        parsed = urllib.parse.urlparse(absolute)
        if parsed.scheme != "https" or parsed.netloc != "ftp.ebi.ac.uk":
            continue
        if urllib.parse.urlparse(absolute).path.rsplit("/", 1)[0] != urllib.parse.urlparse(base).path.rstrip("/"):
            continue
        name = urllib.parse.unquote(parsed.path.rsplit("/", 1)[-1])
        if name and name not in ("..", ".") and "/" not in name:
            names.add(name)
    return sorted(names)


def classify_filename(name):
    low = name.lower()
    if low.endswith((".tsv.gz", ".txt.gz", ".sumstats.gz", ".tsv", ".txt")):
        if low.endswith("-meta.yaml") or low.endswith(".tbi"):
            return "OTHER"
        return "DATA"
    if low.endswith((".yaml", ".yml")):
        return "META"
    return "OTHER"


def choose_data_file(names, accession):
    data = [x for x in names if classify_filename(x) == "DATA"]
    # Original folder listing only. Never infer harmonised file or full QC.
    data.sort(key=lambda x: (
        0 if x.startswith(accession + ".") else 1,
        0 if x.endswith(".gz") else 1,
        x,
    ))
    return data[0] if data else None


def request_url(url, *, headers=None, cap=MAX_PREFIX, method="GET"):
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "IS-Analysis-V3-GwasCatalog-Ftp-BoundedAudit/1.0",
                 "Accept": "*/*", **(headers or {})},
        method=method,
    )
    with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
        final = response.geturl()
        parsed = urllib.parse.urlparse(final)
        if parsed.scheme != "https" or parsed.netloc not in {"ftp.ebi.ac.uk"}:
            raise ValueError("Unexpected FTP redirect destination")
        content = response.read(cap) if method == "GET" else b""
        return {
            "status": response.status,
            "headers": {key.lower(): value for key, value in response.headers.items()},
            "data": content,
            "url": final,
        }


def sniff_text(prefix, filename):
    if filename.lower().endswith(".gz") or prefix.startswith(b"\x1f\x8b"):
        try:
            decompressor = zlib.decompressobj(16 + zlib.MAX_WBITS)
            raw = decompressor.decompress(prefix, MAX_EXPANDED)
        except zlib.error as exc:
            return None, "INVALID_GZIP_PREFIX: " + str(exc)[:100]
        if not raw:
            return None, "NO_GZIP_TEXT"
        return raw.decode("utf-8", errors="replace"), None
    return prefix[:MAX_EXPANDED].decode("utf-8", errors="replace"), None


def inspect_table_header(text):
    if not text:
        return {"header_status": "UNAVAILABLE"}
    # Skip comments only when explicit # line, no truncation assumptions.
    lines = [l for l in text.splitlines() if l.strip() and not l.startswith("##") and not l.startswith("#")]
    if not lines:
        return {"header_status": "UNAVAILABLE"}
    first = lines[0]
    delimiter = "\t" if "\t" in first else ("," if "," in first else None)
    fields = next(csv.reader([first], delimiter=delimiter)) if delimiter else first.split()
    fields = [v.strip().strip('"').lstrip("#") for v in fields]
    norm = {re.sub(r"[^a-z0-9]", "", field.lower()) for field in fields}
    groups = {
        "p_value": {"pvalue", "pval", "p", "pv", "pvaluefixed", "pvaluemeta", "lp"},
        "variant": {"rsid", "snp", "variantid", "chromosome", "chr", "chrom", "hmchrom"},
        "position": {"basepairlocation", "bp", "pos", "position", "hmpos"},
        "effect": {"beta", "effect", "effectsize", "or", "oddsratio", "logor"},
        "se": {"standarderror", "se", "stderr", "sebeta", "sefixed"},
        "effect_allele": {"effectallele", "ea", "a1", "allele1", "alt"},
        "other_allele": {"otherallele", "nea", "a2", "allele2", "ref"},
        "frequency": {"effectallelefrequency", "eaf", "af", "freq", "maf"},
        "sample_size": {"samplesize", "n", "ntotal", "ncase", "ncases", "ncontrol", "ncontrols", "neff"},
    }
    detected = {key: bool(norm.intersection(aliases)) for key, aliases in groups.items()}
    essential = ["p_value", "variant", "effect", "se", "effect_allele", "other_allele"]
    # A raw file can be pre-GWAS-SSF: missing fields means REVIEW, not a data fabrication.
    status = "HEADER_CORE_FIELDS_PRESENT" if all(detected.get(k) for k in essential) else "HEADER_REVIEW"
    return {
        "header_status": status,
        "column_count": len(fields),
        "columns": fields[:90],
        "field_groups": detected,
        "data_row_sampled": bool(len(lines) >= 2),
    }


def verify_accession(ancestry, spec):
    accession = spec["accession"]
    base = accession_directory(accession)
    report = {
        "ancestry": ancestry, "accession": accession,
        "phenotype": spec["phenotype"], "directory": base,
        "file_discovery": "SOURCE_UNAVAILABLE",
        "remote_file": None,
        "prefix_status": "NOT_ATTEMPTED",
        "whole_file_sha256_verified": False,
        "full_sumstats_qc": "HOLD",
        "magma_readiness": "HOLD",
        "ldsc_seg_readiness": "HOLD",
        "disclaimer": "Bounded source inspection only, not full file download, sample count or scientific QC.",
    }
    try:
        listing = request_url(base, cap=MAX_HTML)
        if listing["status"] != 200:
            raise ValueError("Listing HTTP status " + str(listing["status"]))
        names = extract_directory_names(listing["data"].decode("utf-8", errors="replace"), base)
        report["directory_listing"] = names[:100]
        chosen = choose_data_file(names, accession)
        if not chosen:
            report["file_discovery"] = "LISTED_NO_ORIGINAL_DATA"
            return report
        url = urllib.parse.urljoin(base, chosen)
        report["remote_file"] = url
        report["file_discovery"] = "REMOTE_FILE_LISTED"
        meta = [x for x in names if classify_filename(x) == "META"]
        report["metadata_filenames"] = meta[:20]
        try:
            head = request_url(url, method="HEAD")
            report["remote_file_http_status"] = head["status"]
            report["remote_content_length_bytes"] = head["headers"].get("content-length")
            report["remote_etag"] = head["headers"].get("etag")
            report["remote_last_modified"] = head["headers"].get("last-modified")
        except (urllib.error.URLError, ValueError, TimeoutError) as exc:
            report["head_error"] = f"{type(exc).__name__}: {str(exc)[:200]}"
        try:
            data = request_url(url, headers={"Range": f"bytes=0-{MAX_PREFIX - 1}"}, cap=MAX_PREFIX)
            report["range_http_status"] = data["status"]
            report["sampled_compressed_bytes"] = len(data["data"])
            report["prefix_sha256_not_whole_file"] = hashlib.sha256(data["data"]).hexdigest()
            decoded, error = sniff_text(data["data"], chosen)
            if error:
                report["prefix_status"] = "PREFIX_UNREADABLE"
                report["prefix_error"] = error
            else:
                report["header"] = inspect_table_header(decoded)
                report["prefix_status"] = report["header"]["header_status"]
        except (urllib.error.URLError, ValueError, TimeoutError) as exc:
            report["prefix_status"] = "SOURCE_UNAVAILABLE"
            report["prefix_error"] = f"{type(exc).__name__}: {str(exc)[:200]}"
    except (urllib.error.URLError, ValueError, TimeoutError) as exc:
        report["error"] = f"{type(exc).__name__}: {str(exc)[:200]}"
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", default="artifacts/is_eas_gwas_source_discovery")
    args = parser.parse_args()
    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    reports = [verify_accession(a, s) for a, s in EXPECTED.items()]
    document = {
        "audit_time_utc": datetime.now(timezone.utc).isoformat(),
        "type": "REMOTE_BOUNDED_PREFIX_AUDIT_ONLY",
        "full_download_completed": False,
        "full_source_integrity_checked": False,
        "phenotype_cohort_audit_completed": False,
        "reports": reports,
    }
    (output / "source_discovery.json").write_text(
        json.dumps(document, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    lines = [
        "# GIGASTROKE EUR/EAS source discovery (bounded prefix only)",
        "", "| Ancestry | Study | Remote discovery | Header | Full source QC |",
        "|---|---|---|---|---|",
    ]
    for row in reports:
        lines.append("| {ancestry} | {accession} | {discovery} | {header} | HOLD |".format(
            ancestry=row["ancestry"], accession=row["accession"],
            discovery=row["file_discovery"], header=row["prefix_status"]
        ))
    lines.extend(["", "## Located files / sampled header evidence", ""])
    for row in reports:
        lines.append(f"**{row[\'ancestry\']}**: {row.get(\'remote_file\') or \'not found\'}")
        lines.append(f"Remote size (HEAD): {row.get(\'remote_content_length_bytes\', \'unknown\')} bytes; sampled bytes: {row.get(\'sampled_compressed_bytes\', 0)}")
        h = row.get("header", {})
        lines.append(f"Header columns: {\', \'.join(h.get(\'columns\', []))}")
        lines.append(f"Detected fields: {json.dumps(h.get(\'field_groups\', {}), sort_keys=True)}")
        lines.append("")
    lines.extend([
        "", "This step does not verify complete bytes, sample counts, genome-wide",
        "allele distributions, LD, or MAGMA/LDSC-SEG readiness. Treat those as HOLD.",
    ])
    summary = "\n".join(lines) + "\n"
    (output / "source_discovery.md").write_text(summary, encoding="utf-8")
    print(summary)
    step_summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if step_summary:
        with open(step_summary, "a", encoding="utf-8") as stream:
            stream.write(summary)
    # Source unavailability is not a software failure; evidence must still show HOLD.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
