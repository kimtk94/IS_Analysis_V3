#!/usr/bin/env python3
"""Verify EUR GRCh37 LD Score reference mirror, without retaining its 33 MB archive.

Pins Zenodo DOI, file size, published MD5, archive structure, per-chromosome
LD score gzip integrity, and count coverage. This is reference input QC only;
never promotes GWAS exposure/outcome to LDSC scientific GO.
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import io
import json
import math
import re
import tarfile
import urllib.request
from pathlib import Path, PurePosixPath
from urllib.parse import urlsplit

RECORD_API = "https://zenodo.org/api/records/8182036"
DOI = "10.5281/zenodo.8182036"
FILE_URL = "https://zenodo.org/api/records/8182036/files/eur_w_ld_chr.tar.gz/content"
NAME = "eur_w_ld_chr.tar.gz"
SIZE = 33_357_890
EXPECTED_MD5 = "e2f16343c4cfaa76caa7d0c03d26b489"
MAX_ARCHIVE_BYTES = 40_000_000
MAX_ARCHIVE_MEMBERS = 130
MAX_INNER_BYTES = 40_000_000
CHROMOSOMES = set(range(1, 23))
LD_SCORE_PATTERN = re.compile(r"^eur_w_ld_chr/([1-9]|1[0-9]|2[0-2])\.l2\.ldscore\.gz$")
M_PATTERN = re.compile(r"^eur_w_ld_chr/([1-9]|1[0-9]|2[0-2])\.l2\.M$")
M_5_50_PATTERN = re.compile(r"^eur_w_ld_chr/([1-9]|1[0-9]|2[0-2])\.l2\.M_5_50$")


def read_record(opener=urllib.request.urlopen):
    request = urllib.request.Request(
        RECORD_API, headers={"User-Agent": "CKM-Stroke-LDSC-Reference-QC/1.0", "Accept": "application/json"}
    )
    with opener(request, timeout=30) as response:
        raw = response.read(65537)
        if len(raw) > 65536:
            raise ValueError("Zenodo record metadata exceeds cap")
        metadata = json.loads(raw)
    files = [x for x in metadata.get("files", []) if x.get("key") == NAME]
    if len(files) != 1:
        raise ValueError("missing unique pinned Zenodo file")
    f = files[0]
    if (str(metadata.get("id")) != "8182036" or metadata.get("doi") != DOI
        or f.get("checksum") != "md5:" + EXPECTED_MD5
        or f.get("size") != SIZE or f.get("links", {}).get("self") != FILE_URL):
        raise ValueError("Zenodo published metadata differs from pinned expected")
    return {"record": metadata["id"], "doi": DOI, "file": NAME,
            "published_md5": EXPECTED_MD5, "published_size": SIZE,
            "license_status": "SEE_ZENODO_RECORD",
            "publisher": "Zenodo archive mirror of legacy Alkes 1000G EUR LD scores"}


class DigestReader:
    def __init__(self, response, limit):
        self.response = response
        self.limit = limit
        self.nbytes = 0
        self.md5 = hashlib.md5(usedforsecurity=False)
        self.sha256 = hashlib.sha256()

    def read(self, size=-1):
        count = 1_048_576 if size is None or size < 0 else min(size, 1_048_576)
        buf = self.response.read(count)
        self.nbytes += len(buf)
        if self.nbytes > self.limit:
            raise ValueError("compressed archive byte cap exceeded")
        self.md5.update(buf)
        self.sha256.update(buf)
        return buf


def safe_member(name: str):
    path = PurePosixPath(name)
    if name.startswith(("/", "\\")) or ".." in path.parts or "\\" in name:
        raise ValueError("unsafe tar path")
    if not path.parts or path.parts[0] != "eur_w_ld_chr":
        raise ValueError("tar entry outside pinned root")


def gzip_header_and_rows(stream):
    with gzip.GzipFile(fileobj=stream, mode="rb") as gz:
        head = gz.readline(8192)
        if not head or len(head) >= 8192:
            raise ValueError("invalid LD score gzip header")
        header = head.decode("utf-8-sig", errors="strict").strip().split()
        if not {"CHR", "SNP", "BP", "L2"}.issubset(set(header)):
            raise ValueError("LD score header missing expected CHR SNP BP L2 columns")
        lines = 0
        inflated = len(head)
        while True:
            part = gz.read(1_048_576)
            if not part:
                break
            inflated += len(part)
            if inflated > MAX_INNER_BYTES:
                raise ValueError("LD score member exceeds decompression cap")
            lines += part.count(b"\n")
        if lines < 1000:
            raise ValueError("LD score file too few rows")
        return {"row_count_lower_bound": lines, "header": header[:16]}


def scan_tar(response, *, max_bytes=MAX_ARCHIVE_BYTES):
    if getattr(response, "status", 200) != 200:
        raise ValueError("archive transfer must be HTTP 200")
    if urlsplit(response.geturl()).hostname != "zenodo.org":
        raise ValueError("Zenodo archive redirected to unapproved host")
    length = response.headers.get("Content-Length")
    if length and (not length.isdecimal() or int(length) != SIZE):
        raise ValueError("archive Content-Length differs from published size")
    hashed = DigestReader(response, max_bytes)
    coverage = {"ldscore": set(), "M": set(), "M_5_50": set()}
    n_members = 0
    entries = []
    ld_counts = {}
    with tarfile.open(fileobj=hashed, mode="r|gz") as tar:
        for member in tar:
            n_members += 1
            if n_members > MAX_ARCHIVE_MEMBERS:
                raise ValueError("too many archive members")
            safe_member(member.name)
            if not (member.isfile() or member.isdir()):
                raise ValueError("symlinks and special archive members not supported")
            if member.isdir():
                continue
            if member.size > MAX_INNER_BYTES:
                raise ValueError("oversized archive member")
            matches = [(k, p.fullmatch(member.name)) for k, p in
                       (("ldscore", LD_SCORE_PATTERN), ("M", M_PATTERN),
                        ("M_5_50", M_5_50_PATTERN))]
            slot = next(((name, int(m.group(1))) for name, m in matches if m is not None), None)
            if slot is not None:
                k, chrom = slot
                if chrom in coverage[k]:
                    raise ValueError("duplicate chromosome reference entry")
                coverage[k].add(chrom)
                if k == "ldscore":
                    f = tar.extractfile(member)
                    if f is None:
                        raise ValueError("cannot inspect compressed LD score member")
                    ld_counts[str(chrom)] = gzip_header_and_rows(f)
                elif k in ("M", "M_5_50"):
                    f = tar.extractfile(member)
                    if f is None:
                        raise ValueError("cannot inspect M member")
                    data = f.read(256)
                    value = data.strip().split()
                    if len(value) != 1 or not math.isfinite(float(value[0])) or float(value[0]) <= 0:
                        raise ValueError("invalid reference M positive scalar")
            else:
                entries.append({"name": member.name[:120], "size": member.size})
    while True:
        trailing = hashed.read(1024 * 1024)
        if not trailing:
            break
    if hashed.nbytes != SIZE:
        raise ValueError("incomplete archive transfer")
    md5 = hashed.md5.hexdigest()
    if md5 != EXPECTED_MD5:
        raise ValueError("archive does not match Zenodo published MD5")
    missing = {k: sorted(CHROMOSOMES - v) for k, v in coverage.items()}
    # Standard LDSC uses .l2.M_5_50 by default. The .l2.M alternative is
    # needed only with --not-M-5-50; this mirror does not provide that option.
    required_missing = {k: missing[k] for k in ("ldscore", "M_5_50")}
    status = ("PASS_ZENODO_ARCHIVE_AND_22_CHR_LD_REFERENCE"
              if all(not v for v in required_missing.values())
              else "HOLD_MISSING_REQUIRED_LD_COMPONENT")
    return {
        "technical_status": status,
        "archive_sha256": hashed.sha256.hexdigest(),
        "archive_md5": md5,
        "compressed_bytes": hashed.nbytes,
        "members_seen": n_members,
        "chromosome_sets": {k: sorted(v) for k, v in coverage.items()},
        "missing_chromosomes": missing,
        "required_missing_chromosomes": required_missing,
        "ldsc_default_M_5_50_available": not missing["M_5_50"],
        "ldsc_not_M_5_50_option_supported": not missing["M"],
        "archived_chr6_old_must_be_excluded": any("6_old" in x["name"] for x in entries),
        "ldscore_headers": ld_counts,
        "noncanonical_entries": entries[:40],
        "ldscore_archive_retained": False,
        "archive_members_extracted": False,
        "ldsc_h2_rg_performed": False,
        "scientific_status": "HOLD",
        "ready_for_ldsc": False,
    }


def audit(opener=urllib.request.urlopen):
    result = {
        "schema": "ckm-eur-grch37-ldsc-zenodo-archive-qc-v1",
        "checked_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "mirrored_from": "Original Alkes Group 1000G EUR eur_w_ld_chr",
        "ancestry": "EUR",
        "genome_build": "GRCh37 (original archive lineage, verify variant-level coordinate consistency)",
        "source_approved": False, "ready_for_ldsc": False, "scientific_status": "HOLD"
    }
    try:
        result["provenance"] = read_record(opener)
        request = urllib.request.Request(
            FILE_URL, headers={"User-Agent": "CKM-Stroke-LDSC-Reference-QC/1.0"})
        with opener(request, timeout=60) as response:
            result["integrity"] = scan_tar(response)
    except Exception as e:
        result["integrity"] = {
            "technical_status": "HOLD_REFERENCE_DOWNLOAD_OR_SCHEMA",
            "error": type(e).__name__, "reason": str(e)[:240],
            "scientific_status": "HOLD", "ready_for_ldsc": False,
        }
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = audit()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n",
                        encoding="utf-8")
    record = report.get("integrity", {})
    print(json.dumps({
        "status": record.get("technical_status"),
        "sha256": record.get("archive_sha256"),
        "bytes": record.get("compressed_bytes"),
        "missing": record.get("missing_chromosomes"),
        "scientific_status": report["scientific_status"],
    }))
    if record.get("technical_status") != "PASS_ZENODO_ARCHIVE_AND_22_CHR_LD_REFERENCE":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
