#!/usr/bin/env python3
"""Read-only four-locus legacy GWAS/eQTL input QC.

Read the original user-controlled Drive stats.tsv only when --source-root is
supplied. GitHub CI has no Google Drive credentials and validates the frozen
aggregate provenance contract only. Neither mode reruns MR/coloc/SuSiE.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
from pathlib import Path

GENES = {
    "FGF5": (4, 1639, 174013, "65135b9e984c085e945a4b00990e47343ce5c6199d02375ccb029184321ebe35", 175, 293, 438, 21490568),
    "ALDH2": (12, 1796, 195639, "e4c5dae71eea1f9e59aa89f465990d82f3c4b663d977513d73c794c4429f46a6", 584, 233, 454, 25804928),
    "SH3PXD2A": (10, 1557, 169470, "bf02be38e1bf4efbbf6405bdec17c55fe5141b02bf35f6547006de7a52ba01d1", 584, 212, 106, 19393992),
    "COL4A2": (13, 2977, 322971, "c074094bfe123f6421cbe8672b45b20f33e3eb8c57fd324ece8bbead0505dafe", 170, 380, 329, 70900232),
}
COLUMNS = ["snp", "gwas_beta", "gwas_se", "gwas_maf",
           "eqtl_beta", "eqtl_se", "eqtl_maf", "qtl_n"]
REQUIRED_FALSE = (
    "per_variant_effect_n_independently_verified",
    "source_genome_build_independently_verified",
    "historical_posteriors_recomputed_here",
    "MR_instrument_independence_verified",
    "approved_formal_coloc_reanalysis",
)


def guard(condition: bool, msg: str):
    if not condition:
        raise ValueError(msg)


def profile(raw: bytes, chrom: int) -> dict:
    """Stat-only QC; no variant rows leave this function."""
    guard(len(raw) <= 1_000_000, "unexpected source size")
    text = raw.decode("utf-8-sig")
    rd = csv.DictReader(io.StringIO(text), delimiter="\t")
    guard(rd.fieldnames == COLUMNS, "source table column drift")
    count = 0
    seen = set()
    palindrome = 0
    gap = 0
    source_n = set()
    for row in rd:
        count += 1
        guard(count <= 10000 and None not in row and None not in row.values(),
              "unbounded or malformed input")
        parts = row["snp"].split(":")
        guard(len(parts) == 4 and parts[0] == str(chrom)
              and parts[1].isdigit() and int(parts[1]) > 0
              and all(a in ("A", "C", "G", "T") for a in parts[2:])
              and parts[2] != parts[3], "SNP coordinate/allele invalid")
        guard(row["snp"] not in seen, "duplicate SNP")
        seen.add(row["snp"])
        try:
            v = {key: float(row[key]) for key in COLUMNS[1:]}
        except (ValueError, TypeError, OverflowError) as e:
            raise ValueError("non-numeric source statistic") from e
        guard(all(math.isfinite(x) for x in v.values()), "non-finite source statistic")
        guard(v["gwas_se"] > 0 and v["eqtl_se"] > 0,
              "nonpositive original standard error")
        guard(0 < v["gwas_maf"] <= 0.5 and 0 < v["eqtl_maf"] <= 0.5,
              "original MAF out of bounds")
        guard(v["qtl_n"] > 0, "invalid source qtl_n")
        source_n.add(v["qtl_n"])
        if frozenset((parts[2], parts[3])) in (frozenset(("A", "T")), frozenset(("C", "G"))):
            palindrome += 1
        if abs(v["gwas_maf"] - v["eqtl_maf"]) > 0.2:
            gap += 1
    guard(count > 0, "source data empty")
    return {
        "source_rows": count, "unique_variant_keys": len(seen),
        "source_chromosome": chrom,
        "palindromic_variant_count": palindrome,
        "abs_gwas_eqtl_maf_gap_gt_0_2_count": gap,
        "distinct_qtl_n_values": sorted(source_n),
        "source_size_bytes": len(raw),
        "source_sha256": hashlib.sha256(raw).hexdigest(),
    }


def validate_snapshot(data: dict) -> dict:
    guard(isinstance(data, dict)
          and data.get("schema") == "is-legacy-four-locus-private-source-input-qc-20261122-v1",
          "missing original source audit schema")
    guard(data.get("source_data_read_in_github_ci") is False
          and data.get("original_source_per_variant_effects_committed") is False
          and data.get("scientific_status") == "HOLD"
          and data.get("mr_effect_estimates_independently_recomputed") is False
          and data.get("coloc_hypothesis_posteriors_independently_recomputed") is False
          and data.get("strict_per_trait_qtl_N_cohort_provenance_approved") is False
          and data.get("independently_matched_qtl_ld_approved") is False,
          "source CI or scientific causal claims have been incorrectly promoted")
    items = data.get("genes")
    guard(isinstance(items, list) and [x.get("gene") for x in items] == list(GENES),
          "wrong four-gene order or missing locus")
    totals = {"rows": 0, "pal": 0, "gap": 0}
    for item in items:
        gene = item["gene"]
        chrom, rows, size, sha, n, pal, gap, ldsize = GENES[gene]
        expected = {
            "source_chromosome": chrom,
            "source_rows": rows,
            "source_size_bytes": size,
            "source_sha256": sha,
            "unique_variant_keys": rows,
            "per_row_n_constant": n,
            "palindromic_variant_count": pal,
            "abs_gwas_eqtl_maf_gap_gt_0_2_count": gap,
            "invalid_numeric_or_nonpositive_se_rows": 0,
            "observed_maf_outside_0_0p5_rows": 0,
        }
        guard(all(item.get(k) == v for k, v in expected.items()),
              "source row/hash/statistics drift: " + gene)
        guard(item.get("source_columns") == COLUMNS, "source columns drift: " + gene)
        guard(all(item.get(k) is False for k in REQUIRED_FALSE),
              "unjustified MR, coloc or LD approval: " + gene)
        ld = item.get("reference_ld_file_size_only", {})
        guard(ld.get("gtx_eur_ld_f64_bytes") == ldsize
              and ld.get("bbj_eas_ld_f64_bytes") == ldsize
              and ld.get("matrix_n_double") == rows
              and ld.get("expected_n_squared_times_8_bytes") == 8*rows*rows
              and ld.get("binary_content_digest_and_snp_order_verified") is False
              and ld.get("reference_population_matches_qtl_donors") is False
              and ld.get("bbj_cohort_matched_ld_verified") is False,
              "LD size-only comparison misrepresented: " + gene)
        totals["rows"] += rows
        totals["pal"] += pal
        totals["gap"] += gap
    guard(data.get("gene_count") == 4
          and data.get("total_rows") == totals["rows"] == 7969
          and data.get("total_palindromic_rows") == totals["pal"] == 1118
          and data.get("total_abs_gwas_eqtl_maf_gap_gt_0_2") == totals["gap"] == 1327,
          "original four-gene denominators did not reconcile")
    exploratory = data.get("separate_source_model", {})
    guard(exploratory.get("gene") == "FGF5"
          and exploratory.get("exploratory_EBI_QTD000161_BBJ_ABF_baseline_H4_reported") == 0.9132
          and exploratory.get("historical_GTEx_model_independent_replication") is False,
          "exploratory alternative QTL source misleadingly called replication")
    return {
        "schema": "is-four-locus-archive-snapshot-contract-check-v1",
        "technical_status": "PASS_FROZEN_SOURCE_AGGREGATE_CONTRACT_ONLY",
        "source_data_re_downloaded_in_ci": False,
        "source_rows_read_in_ci": 0,
        "historical_GWAS_eQTL_input_rows_verified_in_prior_local_read": totals["rows"],
        "cohort_matched_LD_verified": False,
        "replicated_MR_estimates": 0,
        "independently_recomputed_coloc_results": 0,
        "scientific_status": "HOLD",
        "next_needed": [
            "Preserve exact original source accession/build/alleles and compare raw effects",
            "Approve donor QTL effective N and actual signed SNP-matched study LD",
            "Recompute existing ABF/SuSiE with explicit prior sensitivity before causality language",
            "Separate MR instruments, harmonization, pleiotropy and independent outcome replication",
        ],
    }


def compare_source_dir(source_root: Path, snapshot: dict) -> None:
    """Optional direct private-Drive source re-read; never upload those files."""
    for item in snapshot["genes"]:
        gene = item["gene"]
        path = source_root / (gene + "_stats.tsv")
        if not path.is_file():
            path = source_root / gene / "stats.tsv"
        guard(path.is_file() and not path.is_symlink(),
              "source not found in permitted local staging: " + gene)
        observed = profile(path.read_bytes(), GENES[gene][0])
        for key in ("source_rows", "unique_variant_keys", "source_chromosome",
                    "palindromic_variant_count", "abs_gwas_eqtl_maf_gap_gt_0_2_count",
                    "source_size_bytes", "source_sha256"):
            guard(observed[key] == item[key], "actual private Drive source drift: " + gene + ":" + key)
        guard(observed["distinct_qtl_n_values"] == [float(item["per_row_n_constant"])],
              "source qtl_n not constant as documented: " + gene)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--snapshot", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--source-root", type=Path)
    a = p.parse_args()
    try:
        snapshot = json.loads(a.snapshot.read_text(encoding="utf-8"))
        result = validate_snapshot(snapshot)
        if a.source_root is not None:
            compare_source_dir(a.source_root, snapshot)
            result["source_data_re_read_locally_in_this_run"] = True
        else:
            result["source_data_re_read_locally_in_this_run"] = False
    except (OSError, ValueError, TypeError, KeyError) as exc:
        result = {"schema": "is-four-locus-archive-snapshot-contract-check-v1",
                  "technical_status": "HOLD_SOURCE_SNAPSHOT_OR_SOURCE_FILES",
                  "scientific_status": "HOLD", "source_data_re_downloaded_in_ci": False,
                  "error": str(exc)[:240]}
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False))
    return 0 if result["technical_status"] == "PASS_FROZEN_SOURCE_AGGREGATE_CONTRACT_ONLY" else 1


if __name__ == "__main__":
    raise SystemExit(main())
