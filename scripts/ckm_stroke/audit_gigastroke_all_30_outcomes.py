#!/usr/bin/env python3
"""Verify GIGASTROKE 30 accession metadata and official GRCh37 files.

This is an *outcome registry*, not allele harmonisation, LDSC, MR, or
independent replication. Fail closed for missing ancestry or phenotype data.
Only publisher-hosted GWAS Catalog resources are queried.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
from pathlib import Path
import time
import urllib.error
import urllib.parse
import urllib.request

ACCESSIONS = tuple(f"GCST{i}" for i in range(90104534, 90104564))
API_BASE = "https://www.ebi.ac.uk/gwas/rest/api/studies/"
FTP_BASE = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
            "GCST90104001-GCST90105000/")
SOURCE_TRAITS = ("stroke", "cerebral", "infarction", "ischaemic", "ischemic")
MAX_RESPONSE_BYTES = 100_000
REQUEST_TIMEOUT = 25

def request(url, method="GET", max_bytes=MAX_RESPONSE_BYTES):
    req = urllib.request.Request(
        url, method=method,
        headers={"User-Agent": "TPMI-43321-GIGASTROKE-outcome-registry/1.0",
                 "Accept": "application/json" if method == "GET" else "*/*"})
    with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT) as res:
        if res.status not in (200, 206):
            raise ValueError(f"unexpected HTTP {res.status}")
        resolved = urllib.parse.urlsplit(res.geturl())
        if (resolved.scheme != "https" or resolved.hostname not in
                {"www.ebi.ac.uk", "ftp.ebi.ac.uk"}):
            raise ValueError("unapproved host/redirect")
        data = res.read(max_bytes + 1) if method == "GET" else b""
        if len(data) > max_bytes:
            raise ValueError("oversize metadata response")
        return res.status, dict(res.headers), data

def safe_text(x):
    if isinstance(x, str):
        return x[:750]
    if isinstance(x, dict):
        return {str(k): safe_text(v) for k, v in list(x.items())[:18]}
    if isinstance(x, list):
        return [safe_text(v) for v in x[:12]]
    return x

def extract_study_metadata(raw):
    if not isinstance(raw, dict):
        raise ValueError("GWAS REST response not an object")
    # Never infer phenotype or ancestry from study accession numeric order.
    fields = ("diseaseTrait", "initialSampleSize", "replicationSampleSize",
              "discoverySampleAncestry", "replicationSampleAncestry",
              "ancestry", "reportedTrait", "trait", "publicationInfo",
              "genotypingTechnology", "studyAccession",
              "accessionId", "diseaseTraitId")
    return {
        "response_keys": list(raw)[:40],
        "api_evidence": {k: safe_text(raw[k]) for k in fields if k in raw},
        "raw_text_excerpt": json.dumps(
            {k: raw[k] for k in fields if k in raw},
            ensure_ascii=False, default=str)[:2200],
    }

def single(accession):
    result = {"accession": accession,
              "metadata_status": "HOLD_METADATA",
              "official_grch37_file_status": "HOLD_FILE",
              "format_and_allele_schema_verified": False,
              "sample_independence_verified": False,
              "snp_coverage_verified": False,
              "science": "HOLD"}
    url = API_BASE + accession
    try:
        code, headers, data = request(url)
        raw = json.loads(data.decode("utf-8-sig"))
        result.update(extract_study_metadata(raw))
        result["metadata_status"] = "PASS_REST_JSON"
        result["metadata_url"] = url
    except Exception as exc:
        result["metadata_error"] = (
            f"{type(exc).__name__}: {exc}")[:260]
    file_url = f"{FTP_BASE}{accession}/{accession}_buildGRCh37.tsv.gz"
    result["publisher_grch37_url"] = file_url
    try:
        code, headers, _ = request(file_url, method="HEAD")
        result["official_grch37_file_status"] = "PASS_HEAD"
        result["publisher_file_bytes"] = headers.get("Content-Length")
        result["publisher_file_content_type"] = headers.get("Content-Type")
    except Exception as exc:
        result["publisher_file_error"] = (
            f"{type(exc).__name__}: {exc}")[:260]
    return result

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
        rows = list(pool.map(single, ACCESSIONS))
    passed_metadata = sum(x["metadata_status"] == "PASS_REST_JSON" for x in rows)
    passed_file = sum(x["official_grch37_file_status"] == "PASS_HEAD" for x in rows)
    result = {
        "schema": "GIGASTROKE-30-outcome-registry-v1",
        "source": "Mishra et al., Nature 2022, GWAS Catalog",
        "target_TPMI_phenocode": "433.21",
        "target_rare_lead": "rs187759247",
        "nominal_ancestry_or_case_counts_not_inferred": True,
        "data": rows,
        "candidate_count": len(rows),
        "accessible_catalog_metadata": passed_metadata,
        "accessible_grch37_files": passed_file,
        "technical_status": ("PASS_REGISTRY_PARTIAL" if passed_metadata or passed_file
                            else "HOLD_NO_SOURCE"),
        "scientific_status": "HOLD",
        "reason": ("Study-specific phenotype, ancestry, sample sizes, publisher "
                   "checksums, full SNP coverage and cohort overlap must be checked "
                   "before assigning an outcome to an analysis."),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({
        "technical_status": result["technical_status"],
        "candidate_count": len(rows),
        "accessible_catalog_metadata": passed_metadata,
        "accessible_grch37_files": passed_file,
        "rows": [{k: v for k, v in r.items() if k in
                  ("accession", "metadata_status", "official_grch37_file_status",
                   "api_evidence", "metadata_error", "publisher_file_error")}
                 for r in rows],
        "scientific_status": "HOLD"
    }, ensure_ascii=False))

if __name__ == "__main__":
    main()
