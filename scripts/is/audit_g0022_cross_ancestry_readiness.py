#!/usr/bin/env python3
"""Fail-closed G0022 EAS vs EUR AIS source coverage audit; no inferred variants."""
import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


def open_table(path):
    path = Path(path)
    return gzip.open(path, "rt", newline="") if path.suffix == ".gz" else path.open("r", newline="")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for buf in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(buf)
    return h.hexdigest()


def load_eas(variants_path, summary_path):
    meta = json.loads(Path(summary_path).read_text())
    cs = {v for entries in meta["credible_set_variant_ids"].values() for v in entries}
    variants = {}
    with open_table(variants_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            variant = row["variant_id"]
            if variant in variants:
                raise ValueError(f"duplicate EAS variant: {variant}")
            bits = variant.split(":")
            if len(bits) != 4 or (row["ref"].upper(), row["alt"].upper()) != (bits[2].upper(), bits[3].upper()):
                raise ValueError(f"EAS reference/ALT mismatch: {variant}")
            if row["status"] != "ALT_SIGNED_CONFIRMED":
                raise ValueError(f"unverified EAS orientation: {variant}")
            variants[variant] = row
    if not variants or not cs or not cs.issubset(variants):
        raise ValueError("missing EAS input variants or credible-set SNPs")
    chromosomes = {x["chr"].replace("chr", "") for x in variants.values()}
    if len(chromosomes) != 1:
        raise ValueError("mixed chromosomes in a locus")
    return variants, cs, meta


def collect(eas, cs, eur_path, eur_ld_path=None):
    chrom = next(iter(eas.values()))["chr"].replace("chr", "")
    position_map = defaultdict(list)
    for variant, row in eas.items():
        position_map[int(row["pos"])].append(variant)
    min_pos = min(position_map)
    max_pos = max(position_map)
    matches = defaultdict(list)
    positions_observed = set()
    total_rows, chrom_rows, region_rows = 0, 0, 0
    with open_table(eur_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        expected = {"chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "p_value"}
        if not expected.issubset(set(reader.fieldnames or ())):
            raise ValueError(f"EUR source columns missing: {expected - set(reader.fieldnames or ())}")
        for row in reader:
            total_rows += 1
            if row["chromosome"].removeprefix("chr") != chrom:
                continue
            chrom_rows += 1
            try:
                pos = int(row["base_pair_location"])
            except (ValueError, TypeError):
                continue
            if min_pos <= pos <= max_pos:
                region_rows += 1
            if pos not in position_map:
                continue
            positions_observed.add(pos)
            eur_effect = row["effect_allele"].upper()
            eur_other = row["other_allele"].upper()
            for var in position_map[pos]:
                e = eas[var]
                ref, alt = e["ref"].upper(), e["alt"].upper()
                if (eur_effect, eur_other) == (alt, ref):
                    flip = 1
                elif (eur_effect, eur_other) == (ref, alt):
                    flip = -1
                else:
                    continue
                try:
                    freq = float(row["effect_allele_frequency"])
                    signed_beta = float(row["beta"]) * flip
                    se = float(row["standard_error"])
                    p = float(row["p_value"])
                except (TypeError, ValueError):
                    continue
                if not (0 <= freq <= 1 and se > 0 and 0 <= p <= 1):
                    continue
                matches[var].append((freq if flip == 1 else 1 - freq, signed_beta, se, p, flip))
    outrows = []
    statuses = Counter()
    for var, e in sorted(eas.items(), key=lambda z: int(z[1]["pos"])):
        found = matches[var]
        pos = int(e["pos"])
        if len(found) == 1:
            status = "EXACT_ALLELE_PAIR_MATCHED"
        elif len(found) > 1:
            status = "AMBIGUOUS_MULTIPLE_MATCHES"
        elif pos in positions_observed:
            status = "POSITION_PRESENT_ALLELES_UNMATCHED"
        else:
            status = "POSITION_ABSENT_EUR_SOURCE"
        statuses[status] += 1
        good = found[0] if len(found) == 1 else None
        outrows.append({
            "variant_id": var,
            "is_eas_credible_set": int(var in cs),
            "eas_reference_alt_af": e["reference_alt_af"],
            "eas_gwas_alt_eaf": e["gwas_alt_eaf"],
            "eur_source_coverage": status,
            "eur_gwas_alt_eaf": format(good[0], ".8g") if good else "",
            "eur_alt_beta": format(good[1], ".8g") if good else "",
            "eur_se": format(good[2], ".8g") if good else "",
            "eur_p": format(good[3], ".8g") if good else "",
            "eur_effect_allele_sign": good[4] if good else "",
        })
    cs_rows = [r for r in outrows if r["is_eas_credible_set"]]
    cs_matched = sum(r["eur_source_coverage"] == "EXACT_ALLELE_PAIR_MATCHED" for r in cs_rows)
    total_matched = statuses["EXACT_ALLELE_PAIR_MATCHED"]
    ld_ok = bool(eur_ld_path and Path(eur_ld_path).is_file())
    return outrows, {
        "audit": "IS_G0022_EAS_EUR_GWAS_SOURCE_COVERAGE_V1",
        "scientific_status": "SOURCE_COVERAGE_AUDIT_NOT_FINE_MAPPING",
        "eur_source_gwas_total_rows": total_rows,
        "eur_source_chromosome_rows": chrom_rows,
        "eur_source_g0022_region_rows": region_rows,
        "eas_reference_variant_count": len(eas),
        "eas_credible_set_variant_count": len(cs),
        "eur_exact_matched_eas_variants": total_matched,
        "eur_exact_matched_eas_fraction": total_matched / len(eas),
        "eur_position_only_or_mismatched": statuses["POSITION_PRESENT_ALLELES_UNMATCHED"],
        "eur_missing_positions": statuses["POSITION_ABSENT_EUR_SOURCE"],
        "eur_ambiguous_multiple_matches": statuses["AMBIGUOUS_MULTIPLE_MATCHES"],
        "eur_matched_eas_credible_set_variants": cs_matched,
        "eur_missing_or_unmatched_cs_variants": len(cs) - cs_matched,
        "eur_ancestry_matched_signed_ld_provided": ld_ok,
        "multi_ancestry_direct_cs_ready": cs_matched == len(cs) and total_matched >= 100 and ld_ok,
        "ready_for_multi_ancestry_fine_mapping": False,
        "readiness_reason": "No verified ancestry-specific EUR LD, per-variant N/cohort overlap, and EAS/EUR common variant QC; CS coverage also required",
        "warning": "No imputation or invented EUR effect. Missing source variant is not monomorphic genotype proof. Position and allele matching does not establish study independence.",
    }


def build(args):
    variants, cs, model = load_eas(args.eas_variants, args.eas_summary)
    detail, summary = collect(variants, cs, args.eur_original, args.eur_signed_ld)
    summary["sources"] = {
        "eas_variants_path": str(Path(args.eas_variants).resolve()),
        "eas_summary_path": str(Path(args.eas_summary).resolve()),
        "eur_original_path": str(Path(args.eur_original).resolve()),
        "eas_variants_sha256": digest(args.eas_variants),
        "eas_summary_sha256": digest(args.eas_summary),
        "eur_original_sha256": digest(args.eur_original),
    }
    summary["source_model"] = {
        "model": "EAS AIS full-locus SuSiE-RSS",
        "n_snps_reported": model.get("n_snps"),
        "n_credible_sets": model.get("n_credible_sets"),
        "approximate_n_eff": model.get("approximate_n_eff"),
        "warning": "1000G EAS n=504 external LD, rank<=503; GWAS N_eff is not verified per variant",
    }
    if args.outdir:
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        tsv = outdir / "G0022_EAS_EUR_SOURCE_VARIANT_COVERAGE.tsv"
        with tsv.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(detail[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(detail)
        (outdir / "G0022_EAS_EUR_READINESS_SUMMARY.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--eas-variants", type=Path, required=True)
    p.add_argument("--eas-summary", type=Path, required=True)
    p.add_argument("--eur-original", type=Path, required=True)
    p.add_argument("--eur-signed-ld", type=Path, default=None, help="Validated EUR signed LD path (not inferred from EUR GWAS)")
    p.add_argument("--outdir", type=Path, required=True)
    args = p.parse_args()
    print(json.dumps(build(args), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
