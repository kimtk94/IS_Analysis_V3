#!/usr/bin/env python3
"""Normalize public CKDGen/GLGC GWAS with exact loss classification.

Technical transformations only: no build conversion, allele harmonization,
or inference. Every excluded source row is counted by its first failing check.
"""
import argparse
import csv
import datetime
import gzip
import hashlib
import json
import math
import os
from pathlib import Path

MAP = {
    "ckdgen2019": {"CHR": "Chr", "POS": "Pos_b37", "EA": "Allele1",
                   "NEA": "Allele2", "BETA": "Effect", "SE": "StdErr", "P": "P-value"},
    "glgc2021": {"CHR": "CHROM", "POS": "POS_b37", "EA": "ALT",
                 "NEA": "REF", "BETA": "EFFECT_SIZE", "SE": "SE", "P": "pvalue"},
}
OUT_HEADER = ("CHR", "POS", "EA", "NEA", "BETA", "SE", "P")
AUTOSOMES = {str(i) for i in range(1, 23)}
EXCLUDE_KEYS = ("malformed_row", "missing_value", "non_autosomal", "invalid_position",
                "non_snp", "invalid_numeric", "non_finite", "p_underflow",
                "invalid_p", "invalid_se")
MISSING = {"", "NA", "N/A", ".", "NULL", "NONE"}


def open_text(path, mode="rt"):
    return gzip.open(path, mode, encoding="utf-8-sig", newline="") if str(path).endswith(".gz") else open(
        path, mode, encoding="utf-8-sig", newline="")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def classify(row, colmap):
    if not row or not isinstance(row, dict):
        return "malformed_row", None
    try:
        source = {k: row[v] for k, v in colmap.items()}
    except KeyError:
        return "malformed_row", None
    if any(v is None or not isinstance(v, str) for v in source.values()):
        return "malformed_row", None
    if any(v.strip().upper() in MISSING for v in source.values()):
        return "missing_value", None
    chrom = source["CHR"].strip().removeprefix("chr")
    if chrom not in AUTOSOMES:
        return "non_autosomal", None
    try:
        pos = int(source["POS"])
    except ValueError:
        return "invalid_position", None
    if pos <= 0:
        return "invalid_position", None
    ea, nea = source["EA"].strip().upper(), source["NEA"].strip().upper()
    if any(len(allele) != 1 or allele not in "ACGT" for allele in (ea, nea)) or ea == nea:
        return "non_snp", None
    try:
        beta, se, p = map(float, (source["BETA"], source["SE"], source["P"]))
    except ValueError:
        return "invalid_numeric", None
    if not all(map(math.isfinite, (beta, se, p))):
        return "non_finite", None
    if p == 0:
        return "p_underflow", None
    if not (0 < p <= 1):
        return "invalid_p", None
    if se <= 0:
        return "invalid_se", None
    return None, {"CHR": chrom, "POS": pos, "EA": ea, "NEA": nea,
                  "BETA": beta, "SE": se, "P": p}


def ckdgen_rows(input_fh):
    fields = input_fh.readline().strip().split()
    def iterator():
        for line in input_fh:
            vals = line.split()
            yield dict(zip(fields, vals)) if len(vals) == len(fields) else {}
    return fields, iterator()


def normalize(src, dest, style, expected_sha256=""):
    src, dest = Path(src), Path(dest)
    if style not in MAP:
        raise ValueError("unsupported GWAS source")
    if src.resolve() == dest.resolve():
        raise ValueError("input and output must differ")
    if not src.is_file():
        raise FileNotFoundError(src)
    if dest.suffix != ".gz":
        raise ValueError("output must end in .gz")
    original_digest = sha256(src)
    if expected_sha256 and original_digest.lower() != expected_sha256.lower():
        raise ValueError("source SHA256 mismatch")
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".partial")
    excluded = {key: 0 for key in EXCLUDE_KEYS}
    examples = {key: [] for key in EXCLUDE_KEYS}
    count = good = 0
    try:
        with open_text(src) as input_fh, gzip.open(tmp, "wt", newline="", encoding="utf-8") as output_fh:
            if style == "ckdgen2019":
                fields, reader = ckdgen_rows(input_fh)
            else:
                reader = csv.DictReader(input_fh, delimiter="\t")
                fields = reader.fieldnames
            missing_headers = set(MAP[style].values()) - set(fields or [])
            if missing_headers:
                raise ValueError("source headers missing: " + str(sorted(missing_headers)))
            writer = csv.DictWriter(output_fh, fieldnames=OUT_HEADER, delimiter="\t")
            writer.writeheader()
            for row in reader:
                count += 1
                reason, normalized = classify(row, MAP[style])
                if reason:
                    excluded[reason] += 1
                    if len(examples[reason]) < 3:
                        examples[reason].append({"row": count, "rsid": row.get("RSID", row.get("rsID", "")),
                                                 "raw_values": {key: row.get(value) for key, value in MAP[style].items()}})
                    continue
                writer.writerow(normalized)
                good += 1
        if good == 0:
            raise ValueError("no usable autosomal SNPs")
        tmp.replace(dest)
    finally:
        if tmp.exists():
            tmp.unlink()
    excluded_total = sum(excluded.values())
    if count != good + excluded_total:
        raise AssertionError("inconsistent row accounting")
    return {
        "status": "HOLD" if excluded_total else "PASS",
        "scientific_status": "HOLD", "source": style,
        "input_file": str(src), "input_sha256": original_digest,
        "output_file": str(dest), "output_sha256": sha256(dest),
        "input_rows": count, "output_rows": good, "excluded_total": excluded_total,
        "excluded": excluded, "excluded_examples": {k: v for k, v in examples.items() if v},
        "declared_genome_build": "GRCh37",
        "note": "Original complete stream parsed; source variants removed with reason counts. Build, sample-overlap and independent checksum require separate validation.",
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--source", required=True, choices=sorted(MAP))
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--report", required=True)
    p.add_argument("--expected-sha256", default="")
    args = p.parse_args()
    try:
        result = normalize(args.input, args.output, args.source, args.expected_sha256)
    except (OSError, ValueError, AssertionError, EOFError) as exc:
        result = {"status": "FAIL", "scientific_status": "NOT_RUN", "error": str(exc)}
    result["time_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    result["source_commit"] = os.getenv("GITHUB_SHA", "")
    path = Path(args.report)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result))
    return int(result["status"] == "FAIL")


if __name__ == "__main__":
    raise SystemExit(main())
