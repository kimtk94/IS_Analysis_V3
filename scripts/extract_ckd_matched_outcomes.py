#!/usr/bin/env python3
"""Extract CKD GWAS rows that match a reviewed rsID instrument set.

This stage does not harmonize allele orientation. It preserves each source's
published effect allele, other allele, beta/log-OR, SE and P value so that
allele harmonisation is explicit and auditable in the downstream MR stage.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from pathlib import Path

SCHEMAS = {
    "stanzick_egfr": {
        "delimiter": "\t", "rsid": "RSID", "ea": "Allele1", "oa": "Allele2",
        "beta": "Effect", "se": "StdErr.GC", "p": "P.value.GC", "eaf": "Freq1",
        "n": "n", "chr": "chr", "pos": "pos", "effect_scale": "beta",
    },
    "wuttke_ckd": {
        "delimiter": None, "rsid": "RSID", "ea": "Allele1", "oa": "Allele2",
        "beta": "Effect", "se": "StdErr", "p": "P-value", "eaf": "Freq1",
        "n": "n_total_sum", "chr": "Chr", "pos": "Pos_b37", "effect_scale": "log_odds",
    },
    "wuttke_bun": {
        "delimiter": None, "rsid": "RSID", "ea": "Allele1", "oa": "Allele2",
        "beta": "Effect", "se": "StdErr", "p": "P-value", "eaf": "Freq1",
        "n": "n_total_sum", "chr": "Chr", "pos": "Pos_b37", "effect_scale": "beta",
    },
    "teumer_uacr": {
        "delimiter": None, "rsid": "RSID", "ea": "Allele1", "oa": "Allele2",
        "beta": "Effect", "se": "StdErr", "p": "P-value", "eaf": "Freq1",
        "n": "n_total_sum", "chr": "Chr", "pos": "Pos_b37", "effect_scale": "beta",
    },
    "gorski_egfrcys": {
        "delimiter": ",", "rsid": "rsID", "ea": "allele1", "oa": "allele2",
        "beta": "beta", "se": "se", "p": "pval", "eaf": "freqA1",
        "n": "N", "chr": None, "pos": None, "effect_scale": "beta",
    },
    "chen_egfr": {
        "delimiter": None, "rsid": "SNP", "ea": "Allele1", "oa": "Allele2",
        "beta": "Effect", "se": "StdErr", "p": "P", "eaf": None,
        "n": None, "chr": "CHR", "pos": "BP", "effect_scale": "beta",
    },
    "chen_bun": {
        "delimiter": None, "rsid": "SNP", "ea": "Allele1", "oa": "Allele2",
        "beta": "Effect", "se": "StdErr", "p": "P", "eaf": None,
        "n": None, "chr": "CHR", "pos": "BP", "effect_scale": "beta",
    },
}


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.suffix == ".gz" else path.open("r", encoding="utf-8", errors="replace")


def instrument_rsids(path: Path) -> set[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if "rsid" not in (reader.fieldnames or []):
            raise SystemExit(f"instrument file missing rsid: {reader.fieldnames}")
        return {r["rsid"].strip() for r in reader if r.get("rsid", "").strip()}


def split_line(line: str, delimiter: str | None) -> list[str]:
    if delimiter is None:
        return re.split(r"\s+", line.strip())
    return line.rstrip("\n\r").split(delimiter)


def value(row: dict[str, str], col: str | None) -> str:
    return "" if col is None else row.get(col, "")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--schema", choices=sorted(SCHEMAS), required=True)
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--instruments", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--summary", type=Path, required=True)
    ap.add_argument("--phenotype", required=True)
    ap.add_argument("--ancestry", required=True)
    args = ap.parse_args()

    schema = SCHEMAS[args.schema]
    targets = instrument_rsids(args.instruments)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    seen: set[str] = set()
    matched = 0
    duplicate_rows = 0

    with open_text(args.input) as src, gzip.open(args.output, "wt", encoding="utf-8", newline="") as dst:
        header = split_line(src.readline(), schema["delimiter"])
        missing = [schema[k] for k in ("rsid","ea","oa","beta","se","p") if schema[k] and schema[k] not in header]
        if missing:
            raise SystemExit(f"{args.schema}: missing columns {missing}; header={header}")
        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "phenotype","ancestry","rsid","effect_allele","other_allele",
            "beta","se","p_value","eaf","n","chr","pos","effect_scale","source_schema"
        ])
        total = 0
        for line in src:
            if not line.strip():
                continue
            vals = split_line(line, schema["delimiter"])
            if len(vals) != len(header):
                # Published files occasionally contain malformed trailing records; skip audibly via summary.
                continue
            total += 1
            row = dict(zip(header, vals))
            rsid = value(row, schema["rsid"]).strip()
            if rsid not in targets:
                continue
            if rsid in seen:
                duplicate_rows += 1
            seen.add(rsid)
            writer.writerow([
                args.phenotype,args.ancestry,rsid,
                value(row,schema["ea"]),value(row,schema["oa"]),
                value(row,schema["beta"]),value(row,schema["se"]),value(row,schema["p"]),
                value(row,schema["eaf"]),value(row,schema["n"]),
                value(row,schema["chr"]),value(row,schema["pos"]),
                schema["effect_scale"],args.schema,
            ])
            matched += 1

    summary = {
        "schema": args.schema,
        "phenotype": args.phenotype,
        "ancestry": args.ancestry,
        "instrument_unique_rsids": len(targets),
        "source_rows_parsed": total,
        "matched_rows": matched,
        "matched_unique_rsids": len(seen),
        "duplicate_matched_rows": duplicate_rows,
        "match_fraction": len(seen) / len(targets) if targets else 0,
        "note": "alleles are source-published; orientation harmonisation not yet applied",
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    if matched == 0:
        raise SystemExit("no instrument rsIDs matched outcome")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
