#!/usr/bin/env python3
"""Reconcile 36 actual P3A QTL assays against pinned original-source CI evidence.

This is a read-only, fail-closed *input readiness* audit, not a coloc run.
Archived cis counts are gene-assay row counts, not unique independent SNPs.
LD is not a prerequisite for coloc.abf; ancestry/allele-matched signed LD
is required to evaluate coloc.susie.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

KEYS = ("region_id", "QTL_dataset", "gene_id")
EXPECTED_QTL = {"QTD000131": 10, "QTD000171": 11, "QTD000434": 15}
PINNED = {
    "recover": 45019987189,
    "gwas": 45020605713,
    "eas24": 45022609773,
}
REQUIRED_ABF = (
    "QTL_ORIGINAL_FULL_CIS_SIGNED_BETA_SE_MAF",
    "QTL_DONOR_N_ANCESTRY_AND_TISSUE",
    "GWAS_ORIGINAL_FULL_CIS_SIGNED_BETA_SE_MAF",
    "GWAS_SOURCE_EXACT_CASE_CONTROL_N_AND_PHENOTYPE",
    "GENOME_BUILD_REF_ALT_AND_EFFECT_ALLELE_HARMONIZATION",
    "MATCHED_UNBIASED_FULL_CIS_VARIANT_SET",
    "ARCHIVED_SOURCE_STAT_MEMBERSHIP_AND_PROVENANCE",
)
REQUIRED_SUSIE_EXTRA = (
    "ANCESTRY_AND_ALLELE_MATCHED_SIGNED_LD_MATRIX",
    "CONSISTENT_SNP_ORDER_AND_FINE_MAPPING_QC",
)


def require(condition: bool, reason: str) -> None:
    if not condition:
        raise ValueError(reason)


def group_key(item: dict) -> tuple:
    return tuple(item.get(k) for k in KEYS)


def load(path: Path) -> dict:
    result = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(result, dict), "expected JSON object: " + str(path))
    return result


def audit(manifest: dict, recovery: dict, original: dict, eas24: dict) -> dict:
    groups = manifest.get("groups", [])
    require(manifest.get("schema") == "is-p3a-36-real-group-keys-from-legacy-393-v1"
            and manifest.get("scientific_status") == "HOLD"
            and manifest.get("origin", {}).get("commit_sha")
                == "0bf646fd6f29322d17da878253f37def24c01c74",
            "source manifest identity not pinned")
    require(len(groups) == 36 and len({group_key(x) for x in groups}) == 36
            and len({x.get("gene_id") for x in groups}) == 20
            and Counter(x.get("QTL_dataset") for x in groups) == EXPECTED_QTL
            and all(x.get("region_id") == "IS_XDATA_G0015"
                    and x.get("approved_for_new_coloc") is False
                    and x.get("original_P3A_source_stat_positive_proven") is False
                    for x in groups), "36 actual assay identities or science state drift")
    require(manifest.get("denominators", {}).get(
                "per_group_original_server_source_stats_positive_verified") == 0,
            "legacy source-stat 25/11 membership was incorrectly assigned")
    by_manifest = {group_key(x): x for x in groups}

    recovered = recovery.get("groups", [])
    require(recovery.get("schema") == "is-p3a-restored-36-pinned-legacy-source-audit-v1"
            and recovery.get("technical_status", "").startswith("PASS_")
            and recovery.get("scientific_status") == "HOLD"
            and recovery.get("restored_group_key_count") == 36
            and recovery.get("distinct_gene_id_count") == 20
            and recovery.get("public_BH05_eGene_assay_count") == 25
            and recovery.get("fully_approved_new_P3A_coloc_groups") == 0
            and recovery.get("new_P3A_coloc_computations") == 0
            and len(recovered) == 36
            and {group_key(x) for x in recovered} == set(by_manifest),
            "pinned original QTL source readback missing/modified")
    by_recovery = {group_key(x): x for x in recovered}
    coverage = recovery.get("source_only_25_group_fullcis_coverage", {})
    require(coverage.get("sum_tested_fullcis_rows") == 142196
            and coverage.get("sum_BBJ_EAS_nonpal_exact_gene_assay_rows") == 59090
            and coverage.get("sum_nonpal_QTL_EAS_AF_delta_le_0p1_gene_assay_rows") == 31313,
            "old assay-variant source coverage drift")

    original_rows = original.get("per_group_status", [])
    require(original.get("schema") == "is-p3a-chr10-original-gwas-combined-36-qtl-groups-v1"
            and original.get("technical_status")
                == "PASS_BOTH_PUBLISHER_ORIGINAL_GWAS_CHR10_CONTEXT_SCAN"
            and original.get("scientific_status") == "HOLD"
            and original.get("publisher_accessions")
                == {"EUR": "GCST90104540", "EAS": "GCST90104545"}
            and original.get("EUR_chr10_context_valid_original_variant_rows") == 4397
            and original.get("EAS_chr10_context_valid_original_variant_rows") == 3913
            and original.get("real_group_keys") == 36
            and original.get("real_gene_ids") == 20
            and original.get("distinct_analysed_GWAS_region_count") == 1
            and original.get("fully_approved_per_group_gwas_qtl_source_matches") == 0
            and original.get("case_control_N_approved") is False
            and original.get("coloc_abf_executed") == 0
            and original.get("coloc_susie_executed") == 0
            and len(original_rows) == 36
            and {group_key(x) for x in original_rows} == set(by_manifest)
            and all(x.get("QTL_GWAS_per_variant_same_source_effect_verified") is False
                    for x in original_rows),
            "official original chr10 source evidence unavailable or incorrectly promoted")

    require(eas24.get("schema") == "is-gws-eas-original-24-variant-cross-region-readback-v1"
            and eas24.get("source_audit_status")
                == "PASS_ORIGINAL_ROWS_RECONCILED_TO_PINNED_SOURCE_JSON"
            and eas24.get("scientific_status") == "HOLD"
            and eas24.get("original_source") == "GCST90104545"
            and eas24.get("nonpalindromic_leads_tested") == 24
            and eas24.get("signed_beta_p_eaf_reproduced") == 24
            and eas24.get("concordant_chromosome_count") == 16
            and eas24.get("P3A_QTL_gene_assay_count") == 36
            and eas24.get("P3A_exact_QTL_GWAS_source_rows_paired") == 0
            and eas24.get("P3A_coloc_abf") == 0
            and eas24.get("P3A_coloc_susie") == 0
            and eas24.get("source_GWAS_case_fraction_N_approved") is False
            and {x.get("region_id") for x in
                 eas24.get("held_palindromic_leads", [])}.issuperset({"IS_XDATA_G0015"}),
            "24 cross-chrom signed lead audit not an eligible chr10 assay provenance")
    results = []
    for key in sorted(by_manifest):
        m, r = by_manifest[key], by_recovery[key]
        has_bh = m["official_studywide_eGene_BH05_assay"]
        require((r.get("source_BH05_eGene_assay") == "YES") is has_bh,
                "public eGene membership differs across immutable sources: " + repr(key))
        require(r.get("complete_cis_tested_SNP_from_BH05_source")
                == m.get("official_full_cis_tested_variant_count"),
                "source fullcis coverage not equal to original record: " + repr(key))
        require(r.get("source_stat_group_membership")
                == "UNVERIFIED_25_OF_36_AGGREGATE_ONLY",
                "old source-positive 25 falsely assigned to particular gene")
        blockers = list(REQUIRED_ABF)
        results.append({
            "region_id": key[0],
            "QTL_dataset": key[1],
            "gene_id": key[2],
            "gene_symbol": m["gene_symbol"],
            "public_BH05_eGene_assay": has_bh,
            "public_BH05_QTL_q": m.get("official_QTL_studywide_BH_q"),
            "archived_fullcis_tested_QTL_assay_variant_rows":
                m.get("official_full_cis_tested_variant_count"),
            "archived_nonpal_BBJ_EAS_overlap_assay_variant_rows":
                m.get("archived_nonpalindromic_BBJ_EAS_SNP_overlap"),
            "publisher_chr10_original_EAS_variant_rows": 3913,
            "publisher_chr10_original_EUR_variant_rows": 4397,
            "individual_old_server_source_stat_membership": "UNVERIFIED",
            "QTL_variant_level_original_effect_SE_N_proven": False,
            "GWAS_variant_level_original_effect_SE_matched_QTL": False,
            "GWAS_exact_Ncase_Ncontrol_and_case_fraction_proven": False,
            "fasta_ref_alt_strand_harmonization_proven": False,
            "matched_complete_cis_SNP_universe_proven": False,
            "signed_matched_LD_proven": False,
            "coloc_abf_missing_requirements": blockers,
            "coloc_susie_additional_requirements": list(REQUIRED_SUSIE_EXTRA),
            "coloc_abf_eligible": False,
            "coloc_susie_eligible": False,
            "coloc_abf_executed": False,
            "coloc_susie_executed": False,
            "status": "HOLD_UNVALIDATED_PER_ASSAY_SOURCE_STATISTICS",
        })
    require(len(results) == 36
            and sum(x["public_BH05_eGene_assay"] for x in results) == 25,
            "new 36-group readiness denominator incorrect")
    by_gene = {}
    for x in results:
        gene = by_gene.setdefault(x["gene_id"], {
            "gene_id": x["gene_id"], "gene_symbol": x["gene_symbol"],
            "assay_count": 0, "public_BH05_eGene_assays": 0,
            "QTL_datasets": [], "coloc_eligible_assays": 0,
        })
        gene["assay_count"] += 1
        gene["public_BH05_eGene_assays"] += int(x["public_BH05_eGene_assay"])
        gene["QTL_datasets"].append(x["QTL_dataset"])
    return {
        "schema": "is-p3a-36-coloc-input-readiness-matrix-v1",
        "technical_status": "PASS_PINNED_REAL_SOURCE_CROSS_AUDIT",
        "scientific_status": "HOLD",
        "pinned_ci_run_ids": PINNED,
        "source_loci": 1,
        "assay_groups": 36,
        "candidate_gene_ids": 20,
        "public_BH05_eGene_assays": 25,
        "source_non_BH05_assays": 11,
        "legacy_per_group_source_stat_positive_membership_known": False,
        "source_coverage_unit": "GENE_ASSAY_VARIANT_ROWS_NOT_UNIQUE_SNP",
        "P3A_36_formal_coloc_abf_eligible": 0,
        "P3A_36_formal_coloc_abf_executed": 0,
        "P3A_36_coloc_susie_eligible": 0,
        "P3A_36_coloc_susie_executed": 0,
        "cross_chrom_24_leads_NOT_same_as_chr10_fullcis": True,
        "ABF_does_NOT_require_LD": True,
        "SuSiE_requires_matched_LD": True,
        "groups": results,
        "genes": sorted(by_gene.values(), key=lambda x: x["gene_id"]),
    }


def write_files(result: dict, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "readiness.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    cols = ("region_id", "QTL_dataset", "gene_id", "gene_symbol",
            "public_BH05_eGene_assay", "public_BH05_QTL_q",
            "archived_fullcis_tested_QTL_assay_variant_rows",
            "individual_old_server_source_stat_membership",
            "coloc_abf_eligible", "coloc_susie_eligible", "status")
    with (out_dir / "readiness.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=cols, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(result.get("groups", []))
    summary = [
        "# IS P3A — all 36 QTL assays scientific input readiness",
        "",
        "## Outcome",
        "",
        "- Technical cross-audit: **" + result["technical_status"] + "**.",
        "- Formal coloc.abf eligible/executed: **0/36 / 0**.",
        "- Formal coloc.susie eligible/executed: **0/36 / 0**.",
        "- 36 gene-assay groups, 20 distinct genes, ONE chr10 GWAS locus.",
        "- 25/36 public BH05 eGene source matches; 11 without this public record.",
        "- Publisher EUR chr10 rows 4,397; EAS chr10 rows 3,913.",
        "- Independently audited 24 EAS nonpal lead SNPs across 16 chromosomes are "
        "not 36-group per-variant cis harmonization.",
        "- Legacy server 25 source-stat-positive vs 11 absent: individual "
        "membership **not recovered**, irrespective of coincident BH05 25/11.",
        "",
        "## Statistical gates",
        "",
        "- coloc.abf requires harmonized complete cis QTL+GWAS effect/SE/MAF, "
        "true sample N and binary GWAS case fraction, consistent alleles/build, "
        "prior sensitivity and single-causal-variant interpretation. "
        "No LD matrix required for the standard coloc.abf model.",
        "- coloc.susie additionally requires signed, allele-ordered "
        "ancestry-appropriate LD and fine-mapping diagnostics.",
        "- Neither result nor causal attribution can be produced from source "
        "file integrity and old QTL association counts.",
        "",
        "## Source evidence",
        "",
        "- Original 36-group recovery CI: #45019987189",
        "- Official EUR/EAS chr10 publisher scan CI: #45020605713",
        "- 24 nonpal EAS original source independent readback CI: #45022609773",
        "",
        "## Next",
        "",
        "Provide genuine per-assay full cis QTL effect/SE/MAF and donor N, "
        "original GWAS exact phenotype and sample N, reference FASTA allele "
        "normalization, and complete SNP intersection. Re-run audit before "
        "any eligible assay is passed to coloc.abf.",
        "",
    ]
    (out_dir / "summary.md").write_text("\n".join(summary), encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    for name in ("manifest", "recovery", "original", "eas24"):
        ap.add_argument("--" + name, required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()
    try:
        result = audit(*(load(getattr(args, k)) for k in
                         ("manifest", "recovery", "original", "eas24")))
    except (ValueError, KeyError, TypeError, OSError, json.JSONDecodeError) as exc:
        result = {
            "schema": "is-p3a-36-coloc-input-readiness-matrix-v1",
            "technical_status": "HOLD_MISSING_OR_CONFLICTING_SOURCE",
            "scientific_status": "HOLD",
            "P3A_36_formal_coloc_abf_eligible": 0,
            "P3A_36_formal_coloc_abf_executed": 0,
            "P3A_36_coloc_susie_executed": 0,
            "groups": [],
            "error": str(exc)[:300],
        }
        args.out_dir.mkdir(parents=True, exist_ok=True)
        (args.out_dir / "readiness.json").write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(result))
        return 1
    write_files(result, args.out_dir)
    print(json.dumps({k: result[k] for k in (
        "technical_status", "scientific_status", "assay_groups",
        "candidate_gene_ids", "P3A_36_formal_coloc_abf_eligible")},
        ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
