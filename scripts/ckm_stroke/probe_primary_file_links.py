#!/usr/bin/env python3
"""Inspect links from public GWAS primary distribution pages without downloading binaries."""
import argparse
import json
from html.parser import HTMLParser
from urllib.parse import urljoin, urlsplit
from urllib.request import Request, urlopen
from pathlib import Path

PAGES = {
    "BMI_GIANT": "https://giant-consortium.web.broadinstitute.org/index.php/GIANT_consortium_data_files",
    "SBP_GRASP": "https://www.nhlbi.nih.gov/grasp/full-results",
    "T2D_DIAGRAM": "https://www.diagram-consortium.org/downloads.html",
}
TARGETS = {
    "BMI_GIANT": ("Updated Meta-analysis Locke et al + UKBiobank 2018",
                  "README Summary Statistics for Yengo et al 2018"),
    "SBP_GRASP": ("UKB-ICBPmeta750k_SBPsummaryResults.txt.gz",
                  "Evangelou_README_BPGWAS.txt"),
    "T2D_DIAGRAM": ("T2DGGI GWAS EUR ancestry meta-analysis summary statistics",),
}
MAX_HTML = 750000


class LinkParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
        self.href = None
        self.words = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            self.href = dict(attrs).get("href")
            self.words = []

    def handle_data(self, data):
        if self.href is not None:
            self.words.append(data)

    def handle_endtag(self, tag):
        if tag == "a" and self.href:
            self.links.append((" ".join(" ".join(self.words).split()), self.href))
            self.href = None
            self.words = []


def pull_html(url, opener=urlopen, limit=MAX_HTML):
    req = Request(url, headers={"User-Agent": "CKM-IS-academic-source-audit/1.0"})
    with opener(req, timeout=18) as resp:
        final = resp.geturl()
        if urlsplit(final).scheme != "https" or urlsplit(final).hostname != urlsplit(url).hostname:
            raise ValueError("unapproved source redirect")
        raw = resp.read(limit + 1)
        if len(raw) > limit:
            raise ValueError("HTML size cap exceeded")
        return raw.decode("utf-8", "replace"), final


def inspect_one(key, opener=urlopen):
    html, final = pull_html(PAGES[key], opener)
    p = LinkParser()
    p.feed(html)
    results = []
    for needle in TARGETS[key]:
        candidates = [
            {"text": text, "href": urljoin(final, href)}
            for text, href in p.links
            if needle.lower() in text.lower()
        ]
        results.append({
            "target": needle,
            "count": len(candidates),
            "links": candidates[:4],
            "status": "ANCHOR_FOUND" if len(candidates) == 1 else "HOLD_NO_UNIQUE_ANCHOR",
        })
    return {
        "name": key, "source_url": final, "html_bytes": len(html.encode()),
        "links_checked": results,
        "requires_user_terms": key == "T2D_DIAGRAM",
        "binary_downloaded": False,
        "ready_for_ldsc": False,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    rows = []
    for key in PAGES:
        try:
            rows.append(inspect_one(key))
        except Exception as exc:
            rows.append({"name": key, "status": "SOURCE_UNAVAILABLE",
                         "reason": str(exc)[:250], "ready_for_ldsc": False})
    result = {
        "schema": "ckm-official-source-distribution-link-audit-v1",
        "scientific_status": "HOLD", "sources_approved": False,
        "no_binary_download": True, "results": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
