#!/usr/bin/env python3
"""Extract Sun et al. 2023 UKB-PPP Supplementary Table 16 cis signals.

The workbook has a two-row header. This script preserves the published
conditional association statistics and adds parsed variant coordinates/alleles
for downstream harmonisation. No additional LD pruning is performed here:
ST16 already represents statistically independent SuSiE/conditional signals.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def parse_variant(variant_id: str) -> tuple[str, str, str, str]:
    parts = str(variant_id).split(":")
    if len(parts) < 4:
        return "", "", "", ""
    return parts[0], parts[1], parts[2], parts[3]


def gene_symbol(protein_id: str) -> str:
    return str(protein_id).split(":", 1)[0].strip()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--summary", type=Path, required=True)
    args = ap.parse_args()

    try:
        import openpyxl
    except ImportError as exc:
        raise SystemExit("openpyxl is required: pip install openpyxl") from exc

    wb = openpyxl.load_workbook(args.xlsx, read_only=True, data_only=True)
    if "ST16" not in wb.sheetnames:
        raise SystemExit(f"ST16 not found; sheets={wb.sheetnames}")
    ws = wb["ST16"]
    rows = list(ws.iter_rows(values_only=True))
    if len(rows) < 5:
        raise SystemExit("ST16 has too few rows")

    # Published workbook layout:
    # row 1 title, row 2 note, row 3 grouped headers, row 4 leaf headers.
    header = ["" if v is None else str(v).strip() for v in rows[4]]
    header[0] = "UKBPPP ProteinID"
    index = {name: i for i, name in enumerate(header)}
    required = [
        "UKBPPP ProteinID", "Variant ID", "rsID", "PIP", "log10(BF)",
        "Beta(cond)", "SE(cond)", "-log10(P)(cond)", "ALT freq",
        "Cis/trans", "Test region hg19", "CS size", "CS variant IDs",
        "CS rsIDs", "CS PIPs",
    ]
    missing = [c for c in required if c not in index]
    if missing:
        raise SystemExit(f"ST16 missing columns: {missing}; header={header}")

    out_header = [
        "protein_id", "gene_symbol", "variant_id", "rsid",
        "chrom_hg19", "pos_hg19", "ref_allele", "alt_allele",
        "pip", "log10_bf", "beta_cond", "se_cond", "minus_log10p_cond",
        "alt_freq", "cis_trans", "test_region_hg19",
        "ld_with_other_signal_top_variants_r2", "credible_set_size",
        "credible_set_variant_ids", "credible_set_rsids", "credible_set_pips",
    ]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if args.output.suffix == ".gz" else open
    cis_rows = 0
    proteins: set[str] = set()
    variants: set[str] = set()
    rsids: set[str] = set()

    with opener(args.output, "wt", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(out_header)
        for raw in rows[5:]:
            vals = ["" if v is None else str(v).strip() for v in raw]
            if len(vals) < len(header):
                vals += [""] * (len(header) - len(vals))
            if vals[index["Cis/trans"]].lower() != "cis":
                continue
            protein = vals[index["UKBPPP ProteinID"]]
            variant = vals[index["Variant ID"]]
            rsid = vals[index["rsID"]]
            chrom, pos, ref, alt = parse_variant(variant)
            writer.writerow([
                protein,
                gene_symbol(protein),
                variant,
                rsid,
                chrom,
                pos,
                ref,
                alt,
                vals[index["PIP"]],
                vals[index["log10(BF)"]],
                vals[index["Beta(cond)"]],
                vals[index["SE(cond)"]],
                vals[index["-log10(P)(cond)"]],
                vals[index["ALT freq"]],
                vals[index["Cis/trans"]],
                vals[index["Test region hg19"]],
                vals[index["LD with top variants for other signals in the region (r2)"]],
                vals[index["CS size"]],
                vals[index["CS variant IDs"]],
                vals[index["CS rsIDs"]],
                vals[index["CS PIPs"]],
            ])
            cis_rows += 1
            proteins.add(protein)
            variants.add(variant)
            if rsid:
                rsids.add(rsid)

    if cis_rows != 10750:
        raise SystemExit(f"Expected 10,750 published cis signals, observed {cis_rows}")

    summary = {
        "source": "Sun et al. Nature 2023, Supplementary Table 16",
        "doi": "10.1038/s41586-023-06592-6",
        "sheet": "ST16",
        "filter": "Cis/trans == cis",
        "n_cis_signals": cis_rows,
        "n_unique_protein_ids": len(proteins),
        "n_unique_variant_ids": len(variants),
        "n_unique_rsids": len(rsids),
        "coordinate_build": "GRCh37/hg19 as published in ST16",
        "statistics": "published conditional beta/se and SuSiE PIP; no extra LD pruning",
        "source_workbook": str(args.xlsx),
        "source_workbook_sha256": sha256(args.xlsx),
        "output": str(args.output),
        "output_sha256": sha256(args.output),
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
