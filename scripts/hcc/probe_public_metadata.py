#!/usr/bin/env python3
"""Bounded probe of *public* GEO supplemental metadata; never downloads raw reads."""
import argparse
import collections
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import pathlib
import urllib.request

# Only small public supplements; no SRA/EGA/raw FASTQ or protected metadata.
SOURCES = [
    ("GSE14520", "GSE14520_Extra_Supplement.txt.gz",
     "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/suppl/GSE14520_Extra_Supplement.txt.gz"),
    ("GSE149614", "GSE149614_HCC.metadata.updated.txt.gz",
     "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.metadata.updated.txt.gz"),
]
COMPRESSED_LIMIT = 2 * 1024 * 1024
DECOMPRESSED_LIMIT = 40 * 1024 * 1024
KEYWORDS = ("patient", "subject", "sample", "etiolog", "hbv", "hcv", "mash",
            "recurr", "relapse", "survival", "follow", "event", "time", "tissue")


def bounded_download(url, timeout=20):
    request = urllib.request.Request(url, headers={"User-Agent": "HCC-Audit/1.0 (public metadata only)"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        compressed = response.read(COMPRESSED_LIMIT + 1)
    if len(compressed) > COMPRESSED_LIMIT:
        raise ValueError("compressed metadata exceeds 2MiB limit")
    with gzip.GzipFile(fileobj=io.BytesIO(compressed)) as gz:
        raw = gz.read(DECOMPRESSED_LIMIT + 1)
    if len(raw) > DECOMPRESSED_LIMIT:
        raise ValueError("decompressed metadata exceeds 40MiB limit")
    return raw, hashlib.sha256(compressed).hexdigest()


def inspect(raw):
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = raw.decode("latin-1")
    nonempty = [line for line in text.splitlines() if line.strip()]
    if not nonempty:
        return {"row_count": 0, "column_count": 0, "columns": [], "candidate_columns": [],
                "format_warning": "empty metadata"}
    first = nonempty[0]
    delimiter = "\t" if first.count("\t") >= first.count(",") else ","
    parsed = csv.reader(nonempty, delimiter=delimiter)
    headers = next(parsed)
    if len(headers) > 300:
        raise ValueError("more than 300 columns in metadata")
    headers = [x.strip() for x in headers]
    nonnull = collections.Counter()
    n = 0
    for cells in parsed:
        n += 1
        if n > 250000:
            raise ValueError("metadata exceeds 250000 rows")
        for j, value in enumerate(cells[:len(headers)]):
            if value.strip():
                nonnull[j] += 1
    selected = [i for i, name in enumerate(headers) if any(s in name.lower() for s in KEYWORDS)]
    return {
        "row_count": n, "column_count": len(headers),
        "columns": headers[:80],
        "candidate_columns": [{"name": headers[i], "nonempty_rows": nonnull[i]} for i in selected[:60]],
        "format_warning": "" if len(headers) > 1 else "single-column or unstructured supplement; manual parse needed"
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results/hcc/feasibility/geo_metadata_probe.json")
    args = parser.parse_args()
    records = []
    for accession, filename, url in SOURCES:
        entry = {"accession": accession, "source_file": filename, "url": url}
        try:
            payload, sha = bounded_download(url)
            entry.update({"download_status": "FETCHED_PUBLIC_SUPPLEMENT",
                          "sha256_compressed": sha, "decompressed_bytes": len(payload)})
            entry.update(inspect(payload))
        except Exception as ex:
            # Distinguish inaccessible network/server from a missing study; do not fail other cohorts.
            entry.update({"download_status": "UNVERIFIED", "error": str(ex)[:240]})
        records.append(entry)
    report = {
        "utc_checked": dt.datetime.now(dt.timezone.utc).isoformat(),
        "interpretation": "Schema inspection only, NOT verified patient pairing, clinical endpoint, or etiology counts.",
        "records": records
    }
    p = pathlib.Path(args.output)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"records": [
        {"accession": x["accession"], "status": x["download_status"],
         "rows": x.get("row_count"), "candidate_columns": x.get("candidate_columns", [])[:12],
         "error": x.get("error")} for x in records]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
