#!/usr/bin/env python3
"""Use signed external EAS LD to audit rs671 proxy coverage in EUR AIS GWAS.

Does not estimate EUR LD, effect sharing, imputed association or causality.
"""
import argparse
import csv
import hashlib
import json
import math
import struct
from pathlib import Path

RS671 = "12:112241766:G:A"
THRESHOLDS = (0.95, 0.90, 0.80, 0.50, 0.20, 0.10)
MATCHED = "EXACT_ALLELE_PAIR_MATCHED"


def read_ld_row(path, n, index):
    p = Path(path)
    if p.stat().st_size != n * n * 8:
        raise ValueError(f"LD dimensions inconsistent: {p.stat().st_size}, n={n}")
    with p.open("rb") as f:
        f.seek(index * n * 8)
        payload = f.read(n * 8)
    r = struct.unpack(f"<{n}d", payload)
    if not math.isfinite(r[index]) or abs(r[index] - 1) > 1e-6:
        raise ValueError("LD index/order/diagonal sanity check failed")
    if any(not math.isfinite(x) or abs(x) > 1.00001 for x in r):
        raise ValueError("invalid LD correlation value")
    return r


def audit(variants, coverage, ld, lead=RS671):
    with Path(variants).open(newline="") as f:
        eas = list(csv.DictReader(f, delimiter="\t"))
    if not eas or any(x["status"] != "ALT_SIGNED_CONFIRMED" for x in eas):
        raise ValueError("EAS variants missing or unsigned")
    vids = [r["variant_id"] for r in eas]
    if len(set(vids)) != len(vids) or lead not in vids:
        raise ValueError("missing lead or duplicate EAS variants")
    with Path(coverage).open(newline="") as f:
        eur_rows = list(csv.DictReader(f, delimiter="\t"))
    eur = {x["variant_id"]: x for x in eur_rows}
    if len(eur) != len(eur_rows) or set(eur) != set(vids):
        raise ValueError("EUR coverage is not same exact EAS variant universe")
    idx = vids.index(lead)
    ldrow = read_ld_row(ld, len(eas), idx)
    evidence = []
    for i, (variant, r) in enumerate(zip(vids, ldrow)):
        row = eur[variant]
        evidence.append({
            "variant_id": variant,
            "r_eas_to_rs671": format(r, ".9g"),
            "r2_eas_to_rs671": format(r * r, ".9g"),
            "eas_ref_alt_af": eas[i]["reference_alt_af"],
            "eur_ais_source_coverage": row["eur_source_coverage"],
            "eur_ais_alt_eaf": row["eur_gwas_alt_eaf"],
            "eur_ais_p": row["eur_p"],
        })
    thresholds = {}
    for threshold in THRESHOLDS:
        subset = [v for v in evidence if float(v["r2_eas_to_rs671"]) >= threshold - 1e-12]
        present = [v for v in subset if v["eur_ais_source_coverage"] == MATCHED]
        thresholds[str(threshold)] = {
            "eas_ld_variants": len(subset),
            "eur_source_exact_matched": len(present),
            "eur_source_missing_or_unmatched": len(subset) - len(present),
            "eur_source_coverage_fraction": len(present) / len(subset) if subset else None,
        }
    return evidence, {
        "audit": "IS_G0022_RS671_EAS_LD_PROXY_EUR_SOURCE_COVERAGE_V1",
        "scientific_state": "REFERENCE_BASED_PROXY_SOURCE_COVERAGE_ONLY",
        "anchor": lead,
        "n_eas_reference_variants": len(vids),
        "external_ld_reference": "1000_Genomes_EAS_504",
        "ld_rank_upper_bound": 503,
        "lead_ld_self": ldrow[idx],
        "thresholds": thresholds,
        "interpretation": "The EAS-LD proxy set has reduced EUR AIS summary coverage; this does NOT estimate EUR LD or demonstrate that these variants are monomorphic in EUR.",
        "multi_ancestry_causal_gain_demonstrated": False,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--variants", type=Path, required=True)
    ap.add_argument("--eur-coverage", type=Path, required=True)
    ap.add_argument("--eas-ld-matrix", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    rows, summary = audit(args.variants, args.eur_coverage, args.eas_ld_matrix)
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    with (out / "G0022_RS671_EAS_LD_PROXY_EUR_COVERAGE.tsv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    summary["source_files"] = {name: str(path.resolve()) for name, path in (
        ("eas_variants", args.variants), ("eur_gwas_source_coverage", args.eur_coverage),
        ("eas_reference_signed_ld", args.eas_ld_matrix))}
    (out / "G0022_RS671_EAS_LD_PROXY_SUMMARY.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
