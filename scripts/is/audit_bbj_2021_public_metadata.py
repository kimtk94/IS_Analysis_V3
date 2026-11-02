#!/usr/bin/env python3
"""Read-only public metadata cross-check for BBJ 2021 ischemic stroke GWAS.
This never downloads person-level data, never adopts N into EAS5 and never runs coloc.
"""
import argparse
import hashlib
import html
from html.parser import HTMLParser
import json
import re
import urllib.parse
import urllib.request
from pathlib import Path

SOURCES = {
    "pheweb": ("https://pheweb.jp/pheno/IS", "pheweb.jp"),
    "opengwas": ("https://gwas.mrcieu.ac.uk/datasets/ebi-a-GCST90018644/", "gwas.mrcieu.ac.uk"),
}
EXPECTED = {"ncase": 22664, "ncontrol": 152022, "n_total": 174686,
            "study_accession": "GCST90018644", "build": "HG19/GRCh37"}


class PlainText(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.parts = []
    def handle_data(self, data):
        self.parts.append(data)
    def value(self):
        return re.sub(r"\s+", " ", html.unescape(" ".join(self.parts))).strip()


def inspect_source_html(source_name, page):
    if source_name not in SOURCES:
        raise ValueError("Source not allowlisted")
    text = PlainText()
    text.feed(page)
    val = text.value()
    if source_name == "opengwas" and "Loading..." in val and "GCST90018644" not in val:
        return {"source": source_name, "source_url": SOURCES[source_name][0],
                "status": "HOLD", "reason": "CLIENT_RENDERED_PAGE_NO_ACCESSION_METADATA",
                "public_page_only": True}
    if source_name == "pheweb":
        checks = {
            "correct_phenotype": bool(re.search(r"\b(?:IS:\s*)?Ischemic stroke\b", val, re.I)),
            "cases": bool(re.search(r"\b22,?664\s+cases\b", val, re.I)),
            "controls": bool(re.search(r"\b152,?022\s+controls\b", val, re.I)),
            "primary_study_citation": "Sakaue" in val and "2021" in val,
        }
    else:
        checks = {
            "correct_phenotype": "Ischemic stroke" in val,
            "correct_accession": "GCST90018644" in val,
            "cases": bool(re.search(r"\b22,?664\b", val)),
            "controls": bool(re.search(r"\b152,?022\b", val)),
            "correct_population": "East Asian" in val,
            "correct_build": bool(re.search(r"(?:HG19|GRCh37)", val, re.I)),
        }
    return {"source": source_name, "source_url": SOURCES[source_name][0],
            "content_sha256": hashlib.sha256(page.encode("utf-8")).hexdigest(),
            "checks": checks,
            "status": "METADATA_MATCH" if all(checks.values()) else "HOLD",
            "public_page_only": True}


def fetch_public_html(name, max_bytes=750000):
    if name not in SOURCES:
        raise ValueError("Source not allowlisted")
    url, host = SOURCES[name]
    request = urllib.request.Request(url, headers={"User-Agent": "IS-Research-Metadata-Audit/1.0"})
    with urllib.request.urlopen(request, timeout=18) as response:
        final = urllib.parse.urlsplit(response.url)
        allowed_hosts = {host}
        if name == "opengwas":
            allowed_hosts.add("opengwas.io")  # official successor, still HTML only
        if final.scheme != "https" or final.hostname not in allowed_hosts:
            raise ValueError("Redirect changed trusted origin")
        content = response.read(max_bytes + 1)
        if len(content) > max_bytes:
            raise ValueError("Public page exceeds bounded fetch")
        ctype = response.headers.get_content_type()
        if ctype not in ("text/html", "text/plain"):
            raise ValueError("Unexpected page content type")
    return content.decode("utf-8", errors="replace")


def check_sources(fetcher=fetch_public_html):
    results = {}
    for name in SOURCES:
        try:
            results[name] = inspect_source_html(name, fetcher(name))
        except (OSError, ValueError, UnicodeError) as exc:
            results[name] = {"source": name, "source_url": SOURCES[name][0],
                             "status": "HOLD",
                             "failure_class": type(exc).__name__,
                             "public_page_only": True}
    all_match = all(x["status"] == "METADATA_MATCH" for x in results.values())
    return {
        "schema": "is-bbj2021-public-metadata-crosswalk-v1",
        "source_results": results, "expected_public_study_metadata": EXPECTED,
        "public_source_metadata_agreement": "PASS" if all_match else "HOLD",
        "EAS5_archive_sha256_verified": False,
        "EAS5_original_GWAS_source_reapproved": False,
        "effect_allele_primary_file_manifest_verified": False,
        "case_control_N_adopted_into_coloc": False,
        "independent_GIGASTROKE_replication_verified": False,
        "scientific_status": "HOLD",
        "note": "Official phenotype and accession metadata are not proof of archived EAS5 file identity or harmonized effect alleles."
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    result = check_sources()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("BBJ_PUBLIC_METADATA_AUDIT=" + json.dumps(
        {k: v for k, v in result.items() if k != "source_results"}, sort_keys=True))
    print("BBJ_PUBLIC_SOURCE_STATUS=" + json.dumps(
        {key: {"status":val["status"],"checks":val.get("checks",{}),
               "failure_class":val.get("failure_class"), "reason":val.get("reason")}
         for key,val in result["source_results"].items()}, sort_keys=True))


if __name__ == "__main__":
    main()
