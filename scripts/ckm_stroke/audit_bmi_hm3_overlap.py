#!/usr/bin/env python3
"""Bounded GIANT BMI duplicate and 1000G EUR HapMap3 allele-compatibility audit.

One streaming read of each pinned public input; no raw GWAS or reference files
retained. The evidence is QC of identity/overlap, *not* valid LDSC h2/rg.
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import json
import re
import urllib.request
from collections import Counter
from pathlib import Path
from urllib.parse import urlsplit

from scripts.ckm_stroke.audit_bmi_full_stream import (
    SOURCE as BMI_URL, EXPECTED as BMI_COLUMNS, BoundedHashReader,
    MAX_BYTES as BMI_LIMIT,
)

HM3_REF_COMMIT = "e9f9739c9c74ad5b4c653a6409a703cc005c62dc"
HM3_URL = ("https://raw.githubusercontent.com/MRCIEU/ldsc/" + HM3_REF_COMMIT +
           "/GRCh37/data/hm3/w_hm3.snplist")
EXPECTED_BMI_SHA256 = "0d6ed0ea97870916b830ccae349df94ca6f3cc68c025c79e784729af7f7136a4"
HM3_MAX_BYTES = 25_000_000
HM3_MAX_LINES = 2_000_000
BMI_MAX_LINES = 12_000_000
NUC = {"A": 0, "C": 1, "G": 2, "T": 3}
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
SAMPLES_MAX = 6


def get_url(url: str, opener=urllib.request.urlopen):
    req = urllib.request.Request(url, headers={"User-Agent": "CKM-IS-HM3-QC/1.0"})
    res = opener(req, timeout=60)
    parsed = urlsplit(res.geturl())
    if (parsed.scheme != "https" or parsed.hostname not in
            {"raw.githubusercontent.com", "giant-consortium.web.broadinstitute.org"}):
        res.close()
        raise ValueError("original source redirected outside allowlisted host")
    if getattr(res, "status", 200) != 200:
        res.close()
        raise ValueError("full stream requires HTTP 200")
    return res


def allele_pair(a: str, b: str) -> int:
    if a not in NUC or b not in NUC or a == b:
        raise ValueError("invalid SNP allele pair")
    return (min(NUC[a], NUC[b]) << 2) | max(NUC[a], NUC[b])


def orientation(a: str, b: str, reference: int) -> str:
    try:
        current = allele_pair(a, b)
        complement = allele_pair(COMPLEMENT[a], COMPLEMENT[b])
    except KeyError:
        return "allele_invalid"
    if current == reference and complement == reference:
        return "palindromic_ambiguous"
    if current == reference:
        return "direct_alleles"
    if complement == reference:
        return "strand_complement_possible"
    return "allele_mismatch"


def hm3_catalog(opener=urllib.request.urlopen):
    with get_url(HM3_URL, opener) as res:
        length = res.headers.get("Content-Length")
        if length and (not length.isdecimal() or int(length) > HM3_MAX_BYTES):
            raise ValueError("HM3 reference over byte cap")
        md = hashlib.sha256()
        count_bytes = 0
        rows = 0
        duplicates = 0
        conflicts = 0
        conflicted_ids: set[str] = set()
        values: dict[str, int] = {}
        buf = b""
        while True:
            chunk = res.read(1024 * 1024)
            if not chunk:
                break
            count_bytes += len(chunk)
            if count_bytes > HM3_MAX_BYTES:
                raise ValueError("HM3 reference exceeds byte cap")
            md.update(chunk)
            buf += chunk
            lines = buf.split(b"\n")
            buf = lines.pop()
            for line in lines:
                rows += 1
                if rows > HM3_MAX_LINES:
                    raise ValueError("HM3 row cap exceeded")
                if rows == 1:
                    if line.rstrip(b"\r") != b"SNP\tA1\tA2":
                        raise ValueError("unexpected HM3 header")
                    continue
                parts = line.rstrip(b"\r").decode("utf-8", "strict").split("\t")
                if len(parts) != 3:
                    raise ValueError("HM3 reference row has wrong fields")
                snp, a, b = parts
                if len(snp) > 128 or not snp:
                    raise ValueError("HM3 invalid SNP ID")
                key = allele_pair(a, b)
                if snp in conflicted_ids:
                    duplicates += 1
                elif snp in values:
                    duplicates += 1
                    if values[snp] != key:
                        conflicts += 1
                        conflicted_ids.add(snp)
                        del values[snp]
                else:
                    values[snp] = key
        if buf:
            raise ValueError("HM3 last row missing newline or partial transfer")
        if length and int(length) != count_bytes:
            raise ValueError("HM3 incomplete Content-Length")
    if not values:
        raise ValueError("empty HM3 reference")
    return values, {
        "ref_url": HM3_URL,
        "ref_commit": HM3_REF_COMMIT,
        "ref_sha256": md.hexdigest(),
        "compressed": False,
        "ref_bytes": count_bytes,
        "ref_lines": rows,
        "unique_rsids": len(values),
        "ref_unique_rsids_before_exclusion": len(values) + len(conflicted_ids),
        "duplicate_rsids": duplicates,
        "conflicting_ref_alleles": conflicts,
        "conflicted_unique_rsids_excluded": len(conflicted_ids),
        "excluded_conflicted_rsid_examples": sorted(conflicted_ids)[:SAMPLES_MAX],
        "source_status": ("PASS_REFERENCE_WITH_CONFLICT_EXCLUSIONS" if conflicted_ids
                          else "REF_BYTES_AND_SCHEMA_PASS"),
    }


def inspect_bmi(reference: dict[str, int], opener=urllib.request.urlopen):
    rsids: set[str] = set()
    positions: set[int] = set()
    snp_keys: set[int] = set()
    ctr = Counter()
    chrom = Counter()
    samples: dict[str, list[str]] = {"duplicate_rsids": [], "duplicate_coordinates": [],
                                     "duplicate_allele_keys": [], "allele_mismatch": []}
    with get_url(BMI_URL, opener) as resp:
        reader = BoundedHashReader(resp, BMI_LIMIT)
        length = resp.headers.get("Content-Length")
        if length and (not length.isdecimal() or int(length) > BMI_LIMIT):
            raise ValueError("BMI source byte cap exceeded")
        with gzip.GzipFile(fileobj=reader, mode="rb") as gz:
            head = gz.readline(8192).decode("utf-8-sig").rstrip("\r\n").split("\t")
            if tuple(head) != BMI_COLUMNS:
                raise ValueError("BMI original header mismatch")
            for row_index, line in enumerate(gz, 1):
                if row_index > BMI_MAX_LINES:
                    raise ValueError("BMI row cap exceeded")
                cells = line.decode("utf-8", "strict").rstrip("\r\n").split("\t")
                if len(cells) != 10:
                    raise ValueError("BMI row column mismatch")
                ch, bp, snp, a, b, *_ = cells
                try:
                    c = int(ch)
                    p = int(bp)
                except ValueError as exc:
                    raise ValueError("BMI invalid chr/bp") from exc
                if not (1 <= c <= 22 and 1 <= p <= 300_000_000):
                    raise ValueError("BMI nonautosomal or invalid coordinate")
                pair = allele_pair(a, b)
                coordinate = (c << 29) | p
                allele_key = (coordinate << 4) | pair
                ctr["rows"] += 1
                chrom[str(c)] += 1
                if snp in rsids:
                    ctr["duplicate_rsids"] += 1
                    if len(samples["duplicate_rsids"]) < SAMPLES_MAX:
                        samples["duplicate_rsids"].append(snp)
                else:
                    rsids.add(snp)
                if coordinate in positions:
                    ctr["duplicate_coordinates"] += 1
                    if len(samples["duplicate_coordinates"]) < SAMPLES_MAX:
                        samples["duplicate_coordinates"].append(f"{c}:{p}")
                else:
                    positions.add(coordinate)
                if allele_key in snp_keys:
                    ctr["duplicate_allele_keys"] += 1
                    if len(samples["duplicate_allele_keys"]) < SAMPLES_MAX:
                        samples["duplicate_allele_keys"].append(f"{c}:{p}:{a}/{b}")
                else:
                    snp_keys.add(allele_key)
                ref_pair = reference.get(snp)
                if ref_pair is not None:
                    ctr["hm3_rsids_matched_rows"] += 1
                    label = orientation(a, b, ref_pair)
                    ctr[label] += 1
                    if label == "allele_mismatch" and len(samples["allele_mismatch"]) < SAMPLES_MAX:
                        samples["allele_mismatch"].append(snp)
        if length and reader.count != int(length):
            raise ValueError("BMI stream truncated")
    digest = reader.sha256.hexdigest()
    if digest != EXPECTED_BMI_SHA256:
        raise ValueError("BMI compressed original SHA256 changed from verified source")
    ctr["unique_rsids"] = len(rsids)
    ctr["unique_coordinates"] = len(positions)
    ctr["unique_allele_keys"] = len(snp_keys)
    ctr["hm3_unique_rsids_intersection"] = len(rsids.intersection(reference))
    ctr["hm3_reference_unique_rsids"] = len(reference)
    ctr["hm3_ref_intersection_percent"] = round(
        100 * ctr["hm3_unique_rsids_intersection"] / len(reference), 4
    )
    ctr["hm3_bmi_intersection_percent"] = round(
        100 * ctr["hm3_unique_rsids_intersection"] / len(rsids), 4
    )
    return {
        "source_bmi_sha256": digest,
        "bmi_original_bytes": reader.count,
        "counts": dict(ctr),
        "chromosome_counts": dict(sorted(chrom.items(), key=lambda x: int(x[0]))),
        "limited_examples": samples,
        "technical_status": "PASS_DUPLICATE_AND_HM3_OVERLAP_QUANTIFIED",
        "limitations": "rsID/allele compatibility only; no MAF in full cohort, not LDSC munge/h2/rg.",
    }


def audit(opener=urllib.request.urlopen):
    output = {
        "schema": "ckm-bmi-giant-2018-hm3-overlap-v1",
        "date_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source_approved_for_ldsc": False, "scientific_status": "HOLD",
        "no_original_archive_retained": True,
        "HRS_frequency_is_cohort_EAF": False,
        "ldsc_h2_rg_performed": False,
    }
    try:
        ref, metadata = hm3_catalog(opener)
        output["reference"] = metadata
        if metadata["source_status"] not in ("REF_BYTES_AND_SCHEMA_PASS",
                                              "PASS_REFERENCE_WITH_CONFLICT_EXCLUSIONS"):
            raise ValueError("HapMap3 catalog unresolved allele conflict")
        output["result"] = inspect_bmi(ref, opener)
    except Exception as exc:
        output["result"] = {"technical_status": "HOLD_SOURCE_OR_SCHEMA",
                            "error_type": type(exc).__name__, "error_detail": str(exc)[:200]}
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    a = parser.parse_args()
    report = audit()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"technical_status": report["result"]["technical_status"],
                      "reference": report.get("reference", {}),
                      "result": report["result"].get("counts", {}),
                      "scientific_status": report["scientific_status"]}, ensure_ascii=False))
    if report["result"]["technical_status"] != "PASS_DUPLICATE_AND_HM3_OVERLAP_QUANTIFIED":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
