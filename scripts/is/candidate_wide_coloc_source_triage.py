#!/usr/bin/env python3
"""Candidate-wide IS coloc source triage: 2225 positional genes, 393 assays, 36 ABF flags.

Two modes:
  (a) baseline: validates cited aggregate from previous real source audits;
  (b) --archived-audit: additionally reads an original, already-produced
      is-archived-fasta-abf-audit-v1 JSON including its 36 per-group records.

Neither mode performs MR, coloc, SuSiE, invents metadata, or approves causality.
Optional detailed groups may be independently audited from a trusted original
science/archived_fasta_abf_audit.json; no raw clinical GWAS rows are uploaded.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

BASELINE_SCHEMA = "is-candidate-wide-coloc-baseline-v1"
GROUP_SCHEMA = "is-archived-fasta-abf-audit-v1"
OUTPUT_SCHEMA = "is-candidate-wide-coloc-source-triage-v1"
GENES = ("FGF5", "ALDH2", "SH3PXD2A", "COL4A2")
COLUMNS = (
    "region_id", "QTL_dataset", "gene_id", "stable_group_id",
    "fullcis_archived_rows", "source_flagged_stat_rows",
    "source_stat_fraction", "source_stat_screen",
    "coloc_abf_eligible", "coloc_abf_executed",
    "independent_replication_verified", "scientific_status",
    "blockers",
)

def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError(reason)

def integer(value, what: str) -> int:
    require(type(value) is int and value >= 0, what + " must be nonnegative integer")
    return value

def validate_baseline(b: dict, four: dict) -> dict:
    require(isinstance(b, dict) and b.get("schema") == BASELINE_SCHEMA, "wrong baseline schema")
    require(b.get("scientific_status") == "HOLD" and b.get("do_not_run_unapproved_coloc") is True,
            "baseline scientific flag was promoted")
    ev = b.get("prior_evidence", {})
    require(ev.get("independent_reexecution_in_this_ci") is False
            and ev.get("source_repo") == "kimtk94/IS_Analysis_V3",
            "source audit provenance lost")
    c = b.get("candidate_universe", {})
    require(c.get("all_positional_gene_ids") == 2225
            and c.get("all_2225_are_coloc_eligible") is False,
            "positional gene universe drift")
    components = ("public_qtl_nominal_available", "regional_pilot_unexecuted",
                  "targeted_slices_no_nominal_rows", "not_in_eas_30_regions")
    nums = [integer(c.get(x), x) for x in components]
    require(nums == [149,629,91,1356] and sum(nums) == 2225, "public QTL gene count partition invalid")
    a = b.get("assay_group_universe", {})
    require(a.get("readiness_region_gene_tissue_groups") == 393
            and a.get("archived_source_gene_dataset_groups") == 98
            and a.get("archived_source_gene_assay_snp_rows") == 321538
            and a.get("abf_followup_source_flagged_groups") == 36
            and a.get("abf_followup_groups_with_some_valid_source_stats") == 25
            and a.get("abf_followup_groups_without_valid_source_stats") == 11
            and a.get("qualifying_archived_stat_rows_summed_across_groups") == 68612,
            "readiness assay-count or source-row-count drift")
    for field in ("source_allele_orientations_independently_approved",
                  "full_cis_coverage_verified",
                  "gwas_genuine_beta_se_and_case_control_metadata_verified",
                  "qtl_individual_sample_n_verified"):
        require(a.get(field) is False, "invalid full-source study approval: " + field)
    for field in ("coloc_abf_eligible_groups",
                  "coloc_abf_executed_new_stage_groups",
                  "qtl_cohort_ancestry_matched_ld_groups",
                  "independent_stroke_replication_groups"):
        require(a.get(field) == 0, "unsupported biological group approval: " + field)
    h = b.get("archived_four_locus_panel", {})
    require(h.get("genes") == list(GENES)
            and h.get("archived_source_variant_rows") == 7969
            and h.get("historical_coloc_not_equal_to_current_p3a") is True
            and h.get("current_run_mr_recomputed") == 0
            and h.get("current_run_coloc_recomputed") == 0,
            "legacy four-locus studies accidentally joined as new coloc")
    require(isinstance(four, dict)
            and four.get("schema") == "is-legacy-four-locus-private-source-input-qc-20261122-v1"
            and four.get("scientific_status") == "HOLD"
            and four.get("total_rows") == h["archived_source_variant_rows"]
            and [x.get("gene") for x in four.get("genes", [])] == list(GENES)
            and all(x.get("approved_formal_coloc_reanalysis") is False for x in four["genes"]),
            "original four-locus aggregate was changed or misclassified")
    require(b.get("dataset_roles", {}).get("gene_count_vs_group_count_different") is True,
            "different denominators wrongly merged")
    return a

def validate_per_group(raw: dict, a: dict) -> list[dict]:
    require(isinstance(raw, dict) and raw.get("schema") == GROUP_SCHEMA
            and raw.get("integrity") == "PASS", "per-group original scientific report is invalid")
    require(raw.get("scientific_status") == "ARCHIVED_SOURCE_FLAGS_AND_ABF_READINESS_ONLY",
            "per-group audit scientific status promoted")
    readiness = raw.get("readiness", {})
    require(readiness.get("groups") == 393 and readiness.get("abf_candidate_rows") == 36,
            "not same 393-group source universe")
    gate = raw.get("abf_gate_summary", {})
    require(gate.get("eligible_for_coloc_abf") == 0
            and gate.get("groups_with_source_stat_rows") == 25
            and gate.get("groups_with_no_source_stat_rows") == 11
            and gate.get("source_stat_rows_total") == 68612,
            "per-group P3A source evidence drift")
    screen = raw.get("abf_screen", {})
    require(screen.get("screened") == 36
            and screen.get("coloc_abf_executed") == 0,
            "per-group source report incorrectly claims actual coloc runs")
    source = raw.get("p3a_gwas_source_identity", {})
    require(source.get("accepted_sources") == []
            and source.get("registry_status") == "MISSING_APPROVED_REGISTRY",
            "official GWAS source matching not independently reverified")
    header = raw.get("p3a_source_header_audit", {}).get("field_presence", {})
    for key in ("bbj_effect_se", "eas_effect_se", "bbj_sample_size",
                "eas_sample_size", "stroke_case_fraction"):
        require(header.get(key) is False, "source lacks required GWAS SE/N and case fractions: " + key)

    groups = raw.get("abf_flagged_groups")
    require(isinstance(groups, list) and len(groups) == 36, "36 actual source-keyed groups unavailable")
    ids = set()
    out = []
    for item in groups:
        require(isinstance(item, dict), "invalid group")
        keys = [item.get(k) for k in ("region_id", "QTL_dataset", "gene_id")]
        require(all(isinstance(x, str) and x and len(x) <= 150 for x in keys),
                "missing original region/QTL/gene key")
        tup = tuple(keys)
        require(tup not in ids, "duplicate region/gene/QTL assay")
        ids.add(tup)
        total = integer(item.get("fullcis_archived_rows"), "fullcis_archived_rows")
        stat = integer(item.get("source_flagged_stat_rows"), "source_flagged_stat_rows")
        require(stat <= total and item.get("coloc_abf_eligible") is False
                and item.get("scientific_status") == "ABF_FOLLOWUP_FLAG_ONLY_NOT_COLOC_VALIDATION",
                "candidate falsely called eligible or source rows invalid")
        blockers = item.get("blocking_reasons")
        require(isinstance(blockers, list)
                and "QTL_GWAS_SAMPLE_SIZES_UNVERIFIED" in blockers
                and "INDEPENDENT_ALLELE_REFERENCE_PENDING" in blockers
                and "GWAS_EFFECT_VARIANCE_UNVERIFIED" in blockers,
                "required validation blockers missing")
        stable = hashlib.sha256(json.dumps(tup, ensure_ascii=False).encode()).hexdigest()[:16]
        out.append({
            "region_id":keys[0], "QTL_dataset":keys[1], "gene_id":keys[2],
            "stable_group_id":stable,
            "fullcis_archived_rows":total,
            "source_flagged_stat_rows":stat,
            "source_stat_fraction":round(stat/total,6) if total else None,
            "source_stat_screen":"BASIC_SOURCE_STATS_PRESENT" if stat else "NO_VALID_SOURCE_STATS",
            "coloc_abf_eligible":False,
            "coloc_abf_executed":False,
            "independent_replication_verified":False,
            "scientific_status":"HOLD",
            "blockers":";".join(sorted(set(blockers))),
        })
    require(sum(x["source_flagged_stat_rows"] > 0 for x in out) == 25
            and sum(x["source_flagged_stat_rows"] for x in out) == 68612,
            "per-group original 25 / 68,612 source-stat count mismatch")
    return sorted(out,key=lambda x:(-x["source_flagged_stat_rows"],x["region_id"],x["QTL_dataset"],x["gene_id"]))

def generate(b: dict, four: dict, raw: dict | None) -> tuple[dict, list[dict]]:
    a = validate_baseline(b,four)
    groups = validate_per_group(raw,a) if raw is not None else []
    unique_genes = len({x["gene_id"] for x in groups}) if raw is not None else None
    return {
        "schema":OUTPUT_SCHEMA,
        "technical_status":"PASS_SOURCE_AGGREGATES_AND_PER_GROUP_COHORT"
          if raw is not None else "PASS_SOURCE_AGGREGATE_ONLY_NO_GROUP_IDS",
        "scientific_status":"HOLD",
        "candidate_total_positional_genes":2225,
        "readiness_assay_groups":393,
        "p3a_source_flagged_assay_groups":36,
        "p3a_groups_with_basic_stats":25,
        "p3a_groups_missing_stats":11,
        "p3a_science_approved_groups":0,
        "p3a_coloc_executed_in_this_run":0,
        "unique_genes_in_36_assay_groups":unique_genes,
        "all_36_original_ids_loaded":raw is not None,
        "per_group_rows_output":len(groups),
        "historical_four_locus_panel":list(GENES),
        "historical_four_locus_panel_separate":True,
        "independent_replication_approved":False,
        "candidate_prioritization_is_causality":False,
        "matched_qtl_ld_verified":False,
        "next_gate":"Verify official GWAS/QTL effect alleles, full cis coverage, signed effect variances, cohort-specific N and prior sensitivity before running any new coloc",
    },groups

def main() -> int:
    p=argparse.ArgumentParser()
    p.add_argument("--baseline",type=Path,required=True)
    p.add_argument("--four-locus",type=Path,required=True)
    p.add_argument("--archived-audit",type=Path)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    try:
        baseline=json.loads(a.baseline.read_text(encoding="utf-8"))
        four=json.loads(a.four_locus.read_text(encoding="utf-8"))
        report=json.loads(a.archived_audit.read_text(encoding="utf-8")) if a.archived_audit else None
        summary,groups=generate(baseline,four,report)
    except (KeyError,OSError,TypeError,ValueError) as e:
        summary={"schema":OUTPUT_SCHEMA,"technical_status":"HOLD_INPUT_NOT_VERIFIED",
                 "scientific_status":"HOLD","error":str(e)[:240]}
        groups=[]
    a.out_dir.mkdir(parents=True,exist_ok=True)
    (a.out_dir/"candidate_wide_summary.json").write_text(
        json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    if groups:
        with (a.out_dir/"candidate_group_triage.tsv").open("w",encoding="utf-8",newline="") as f:
            w=csv.DictWriter(f,fieldnames=COLUMNS,delimiter="\t")
            w.writeheader()
            w.writerows(groups)
    print(json.dumps(summary,ensure_ascii=False))
    return 0 if summary["technical_status"].startswith("PASS_") else 1

if __name__=="__main__":
    raise SystemExit(main())
