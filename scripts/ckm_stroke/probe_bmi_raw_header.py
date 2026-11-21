#!/usr/bin/env python3
"""Read a bounded prefix of the official GIANT 2018 corrected BMI GWAS gzip.

No full-file transfer, no archive storage, no approval for LDSC or MR.
"""
from __future__ import annotations

import argparse
import json
import zlib
from pathlib import Path
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

SOURCE = ("https://giant-consortium.web.broadinstitute.org/images/c/c8/"
          "Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz")
README = ("https://giant-consortium.web.broadinstitute.org/images/0/01/"
          "README_summary_statistics_Yengo_et_al_2018.txt")
HOST = "giant-consortium.web.broadinstitute.org"
EXPECTED = ("CHR", "POS", "SNP", "Tested_Allele", "Other_Allele",
            "Freq_Tested_Allele_in_HRS", "BETA", "SE", "P", "N")
MAX_PREFIX = 512 * 1024


def fetch_readme(opener=urlopen):
    req = Request(README, headers={"User-Agent": "CKM-Stroke-GIANT-BMI-QC/1.0"})
    with opener(req, timeout=20) as response:
        if urlsplit(response.geturl()).hostname != HOST:
            raise ValueError("README redirects outside official host")
        data = response.read(12001)
        if len(data) > 12000:
            raise ValueError("README exceeds read cap")
    txt = data.decode("utf-8-sig")
    checks = {"hg19": "Genome build hg19" in txt, "HRS_allele_frequency": "Health and Retirement Study" in txt}
    for col in EXPECTED:
        checks[col] = any(line.strip().startswith(col + ":") for line in txt.splitlines())
    return checks


def probe_gzip(opener=urlopen):
    req = Request(SOURCE, headers={"User-Agent": "CKM-Stroke-GIANT-BMI-QC/1.0",
                                   "Range": "bytes=0-524287"})
    raw_bytes = 0
    visible = b""
    with opener(req, timeout=30) as response:
        final = urlsplit(response.geturl())
        if final.scheme != "https" or final.hostname != HOST:
            raise ValueError("source redirected outside pinned host")
        status = getattr(response, "status", 200)
        if status not in (200, 206):
            raise ValueError("unexpected HTTP status")
        decomp = zlib.decompressobj(16 + zlib.MAX_WBITS)
        while raw_bytes < MAX_PREFIX and b"\n" not in visible:
            data = response.read(min(16384, MAX_PREFIX - raw_bytes))
            if not data:
                break
            raw_bytes += len(data)
            visible += decomp.decompress(data, 65536)
            if len(visible) > 65536:
                raise ValueError("uncompressed header prefix too large")
        meta = {
            "http_status": status,
            "content_length": response.headers.get("Content-Length"),
            "content_range": response.headers.get("Content-Range"),
            "compressed_prefix_bytes_read": raw_bytes,
        }
    if b"\n" not in visible:
        raise ValueError("header not found within bounded prefix")
    header = visible.split(b"\n", 1)[0].decode("utf-8-sig", "replace").rstrip("\r").split("\t")
    expected = list(EXPECTED)
    missing = [x for x in expected if x not in header]
    meta.update({"header_columns": header[:40],
                 "missing_required_columns": missing,
                 "header_has_expected_columns": not missing})
    return meta


def audit(opener=urlopen):
    output = {"schema": "ckm-giant-bmi-2018-prefix-audit-v1",
              "dataset": "Yengo2018 corrected GIANT+UKBB EUR BMI",
              "official_source": SOURCE,
              "official_readme": README,
              "raw_md5_verified": False, "full_gwas_downloaded": False,
              "approved_for_ldsc": False, "scientific_status": "HOLD"}
    try:
        readme_checks = fetch_readme(opener)
        prefix = probe_gzip(opener)
        output.update({"readme_checks": readme_checks, **prefix})
        good = all(readme_checks.values()) and prefix["header_has_expected_columns"]
        output["technical_header_status"] = "PASS_HEADER_ONLY" if good else "HOLD_SCHEMA_DISAGREEMENT"
    except Exception as exc:
        output["technical_header_status"] = "SOURCE_UNAVAILABLE_OR_INVALID"
        output["error_type"] = type(exc).__name__
        output["error_detail"] = str(exc)[:180]
    return output


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    result = audit()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    if result["technical_header_status"] != "PASS_HEADER_ONLY":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
