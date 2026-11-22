#!/usr/bin/env python3
"""Read-only cross-model TPMI 433.21 sample-size provenance gate.

Consumes pinned, previously audited GitHub Actions aggregate JSONs for:
  1. original TPMI PLINK2 full-GWAS QC,
  2. original TPMI SAIGE full-GWAS QC,
  3. original PheWeb links and 2025 primary-paper methods audit.

A technical PASS is NOT an LDSC N/N_eff approval. Never compute an
effective sample size, use one model's N for the other, or write SNP data.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

PLINK_SHA = "ce14daef23f9fadf921030d9ed2429383d927f6da391c620d9abf870045b2e62"
SAIGE_SHA = "9f235bfec6265442b7699accecbbd2a69c5acb369f7ee1bcc9e79ec94285a33e"
PLINK_RUN = 44872892525
SAIGE_RUN = 44907119885
METHODS_RUN = 44906623458
PLINK_OBS = 204099
SAIGE_CASES = 9249
SAIGE_CONTROLS = 304660
PLINK_ROWS = 1739048
SAIGE_ROWS = 15588353
SAIGE_PASS = 8042797
SAIGE_FAIL = 7545556


def require(condition: bool, reason: str) -> None:
    if not condition:
        raise ValueError(reason)


def validate(plink: dict, saige: dict, methods: dict) -> dict:
    require(all(isinstance(x, dict) for x in (plink, saige, methods)),
            "aggregate evidence must contain three JSON objects")
    require(
        plink.get("schema") == "TPMI43321-original-gwas-complete-row-QC-v1"
        and plink.get("technical_status") ==
            "PASS_FULL_SOURCE_SHA256_AND_ALL_ROWS_PARSED_WITH_QC_FLAGS"
        and plink.get("source") == "TPMI_433.21_PheWeb_GRCh38_logistic_PLINK2"
        and plink.get("full_original_compressed_sha256") == PLINK_SHA
        and plink.get("pinned_prior_original_sha256_match") is True
        and plink.get("rows_total_in_original") == PLINK_ROWS,
        "PLINK2 original source/hash/model/QC drift",
    )
    extremes = plink.get("value_min_max", {})
    require(
        extremes.get("min_OBS_CT") == PLINK_OBS
        and extremes.get("max_OBS_CT") == PLINK_OBS
        and plink.get("OBS_CT_is_not_SNP_case_control_effective_N") is True
        and plink.get("scientific_status") == "HOLD"
        and plink.get("ready_for_ldsc") is False
        and plink.get("MR_coloc_LDSC_executed") is False,
        "PLINK2 OBS_CT no longer exact or incorrectly promoted to approved N_eff",
    )
    require(
        saige.get("schema") == "TPMI43321-SAIGE-full-original-casecontrol-QC-v1"
        and saige.get("technical_status") == "PASS_FULL_GZIP_AND_SOURCE_CASE_CONTROL_SCAN"
        and saige.get("source") == "https://pheweb.ibms.sinica.edu.tw/download/saige_gwas/433.21"
        and saige.get("computed_original_sha256") == SAIGE_SHA
        and saige.get("records") == SAIGE_ROWS
        and saige.get("full_gzip_crc_verified") is True,
        "SAIGE original source/hash/model/QC drift",
    )
    pairs = saige.get("most_frequent_case_control_pairs", [])
    require(
        saige.get("distinct_case_control_pairs") == 1
        and len(pairs) == 1
        and pairs[0] == {
            "N_case": SAIGE_CASES,
            "N_ctrl": SAIGE_CONTROLS,
            "count": SAIGE_ROWS,
            "sum": SAIGE_CASES + SAIGE_CONTROLS,
        }
        and saige.get("all_N_pairs_identical_to_phenotype_page") is True
        and saige.get("QC_counts", {}).get("source_explicit_QC_PASS") == SAIGE_PASS
        and saige.get("QC_counts", {}).get("source_explicit_QC_FAIL") == SAIGE_FAIL
        and saige.get("PLINK2_model_OBS_CT") == PLINK_OBS
        and saige.get("SAIGE_case_control_N_are_NOT_PLINK2_model_case_control_N") is True
        and saige.get("SAIGE_LDSC_neff_approved") is False
        and saige.get("scientific_status") == "HOLD",
        "SAIGE per-source N/QC or cross-model separation drift",
    )
    require(
        methods.get("schema") == "TPMI43321-PLINK2-SAIGE-sample-basis-audit-v1"
        and methods.get("TPMI_pheweb_cases") == SAIGE_CASES
        and methods.get("TPMI_pheweb_controls") == SAIGE_CONTROLS
        and methods.get("TPMI_pheweb_total") == SAIGE_CASES + SAIGE_CONTROLS
        and methods.get("source_PLINK2_OBS_CT") == PLINK_OBS
        and methods.get("N_discrepancy") == SAIGE_CASES + SAIGE_CONTROLS - PLINK_OBS
        and methods.get("regression_case_count") is None
        and methods.get("regression_control_count") is None
        and methods.get("per_SNP_case_control_effective_N") is None
        and methods.get("scientific_status") == "HOLD_NO_SOURCE_SUPPORTED_CASE_CONTROL_BREAKDOWN",
        "methodology evidence or unresolved N semantics changed",
    )
    links = methods.get("pheweb", {}).get("data_links", [])
    require(
        methods.get("pheweb", {}).get("status") == "PASS_PHEWEB_LINKS"
        and any(x.get("href") == "/download/gwas_plink2/433.21" for x in links)
        and any(x.get("href") == "/download/saige_gwas/433.21" for x in links),
        "distinct original PheWeb model downloads not verified",
    )
    paper = methods.get("primary_papers", {}).get("Chen_TPMI_PRS_2025", {})
    snippets = paper.get("method_evidence", {}).get("SAIGE", {}).get("excerpts", [])
    require(paper.get("status") == "PASS_PRIMARY_XML" and snippets,
            "primary-paper original methods XML unavailable")
    normalized = re.sub(r"\\s+", " ", " ".join(snippets)).lower().replace("\\u2009", " ")
    require(
        "saige" in normalized and "plink2" in normalized
        and "248,754" in normalized and "unrelated subset" in normalized
        and "1:10" in normalized and "matching" in normalized,
        "published SAIGE/PLINK2 subset and imbalance-method evidence missing",
    )
    return {
        "schema": "TPMI43321-crossmodel-N-LDSC-safety-gate-v1",
        "technical_status": "PASS_PINNED_ORIGINAL_MODELS_AND_PUBLISHED_METHODS",
        "original_runs": {
            "PLINK2": PLINK_RUN, "SAIGE": SAIGE_RUN, "primary_methods": METHODS_RUN
        },
        "primary_methods_source": paper["source"],
        "PLINK2": {
            "model": "unrelated-subset_generalized-linear-model",
            "source_constant_OBS_CT": PLINK_OBS,
            "phenotype_specific_case_control_counts": None,
            "eligible_for_source_neff": False,
        },
        "SAIGE": {
            "model": "relatedness-adjusted_mixed-model",
            "source_N_case": SAIGE_CASES,
            "source_N_control": SAIGE_CONTROLS,
            "source_N_total": SAIGE_CASES + SAIGE_CONTROLS,
            "source_SNP_rows": SAIGE_ROWS,
            "source_explicit_QC_PASS": SAIGE_PASS,
            "source_explicit_QC_FAIL": SAIGE_FAIL,
            "eligible_for_direct_transfer_to_PLINK2": False,
        },
        "difference": {
            "source_total_minus_PLINK2_OBS_CT": SAIGE_CASES + SAIGE_CONTROLS - PLINK_OBS,
            "model_design_difference_supported_by_paper": True,
            "exact_43321_selection_and_subset_composition_verified": False,
            "do_not_assume_exact_cause_of_all_excluded_participants": True,
        },
        "LDSC": {
            "N_eff_approved": False,
            "allow_PLINK2_OBS_CT_as_binary_Neff": False,
            "allow_SAIGE_total_as_PLINK2_Neff": False,
            "allow_unapproved_case_control_formula": False,
            "reference_panel_regression_SNPs_verified": False,
            "h2_executed": False,
            "ready_for_ldsc": False,
        },
        "scientific_status": "HOLD",
        "raw_GWAS_or_sample_data_exported": False,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--plink", required=True, type=Path)
    ap.add_argument("--saige", required=True, type=Path)
    ap.add_argument("--methods", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    try:
        result = validate(
            json.loads(a.plink.read_text(encoding="utf-8")),
            json.loads(a.saige.read_text(encoding="utf-8")),
            json.loads(a.methods.read_text(encoding="utf-8")),
        )
    except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError) as exc:
        result = {
            "schema": "TPMI43321-crossmodel-N-LDSC-safety-gate-v1",
            "technical_status": "HOLD_SOURCE_MODEL_OR_METHODS",
            "error_type": type(exc).__name__,
            "error": str(exc)[:250],
            "scientific_status": "HOLD",
            "ready_for_ldsc": False,
            "N_eff_approved": False,
            "h2_executed": False,
        }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\\n", encoding="utf-8")
    print(json.dumps({
        "technical_status": result["technical_status"],
        "PLINK2_OBS_CT": result.get("PLINK2", {}).get("source_constant_OBS_CT"),
        "SAIGE_source_N": result.get("SAIGE", {}).get("source_N_total"),
        "N_difference": result.get("difference", {}).get("source_total_minus_PLINK2_OBS_CT"),
        "scientific_status": result["scientific_status"],
        "LDSC_ready": result.get("LDSC", {}).get("ready_for_ldsc", False),
        "error": result.get("error"),
    }, ensure_ascii=False))
    return 0 if result["technical_status"] == "PASS_PINNED_ORIGINAL_MODELS_AND_PUBLISHED_METHODS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
