#!/usr/bin/env python3
"""Source-artifact-only LDSC readiness matrix. Never run LDSC or invent N_eff.

Input artifacts must be pinned real GitHub Actions aggregate reports, not GWAS
rows, study participant information, or ad hoc sample-size approximations.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

TRAITS = {
    "AIS": ("GCST90104540", 1_006_932),
    "LAA": ("GCST90104542", 956_929),
    "CE": ("GCST90104541", 996_662),
    "SVD": ("GCST90104543", 967_356),
}
BMI_SOURCE_SNP = 1_007_899
LOCKE_SOURCE_SNP = 1_070_753
PROVENANCE_SHA = "ed8f713827cd3281bea5f88a9cb42a3eabfecb1ee67c912ca7e86f151e4399bf"
REFERENCE_SHA = "9537f00eb0d163a935aaa2cf04b358b7cf21852279b9c7925802526f6060b069"


def check(predicate: bool, message: str) -> None:
    if not predicate:
        raise ValueError(message)


def reject_approval(j: dict, *paths: tuple[str, ...]) -> None:
    for path in paths:
        part = j
        for key in path:
            if not isinstance(part, dict) or key not in part:
                raise ValueError("missing fail-closed flag " + ".".join(path))
            part = part[key]
        if part is not False:
            raise ValueError("unexpected LDSC approval " + ".".join(path))


def validate(ref: dict, prereg: dict, bmi_provenance: dict, tpmi: dict) -> dict:
    check(all(isinstance(x, dict) for x in (ref, prereg, bmi_provenance, tpmi)),
          "expected four original aggregate JSON objects")
    check(ref.get("schema") == "ckm-EUR-LDSC-reference-weight-and-outcome-N-source-v1"
          and ref.get("technical_status") == "PASS_22CHR_LD_REFERENCE_AND_SOURCE_NEFF_HOLD_AUDITED"
          and ref.get("reference_BMI_candidate_count") == BMI_SOURCE_SNP
          and ref.get("reference_BMI_candidate_22chr_merge_count") == BMI_SOURCE_SNP
          and ref.get("reference_BMI_candidate_also_in_HM3") == BMI_SOURCE_SNP,
          "source 22chr reference SNP candidate evidence invalid")
    check(ref.get("scientific_status") == "HOLD"
          and ref.get("GC_intercept_and_h2_NOT_RUN") is True
          and ref.get("rg_MR_coloc_NOT_RUN") is True
          and ref.get("case_control_N_eff_per_variant_known") is False
          and ref.get("outcome_post_munge_SNP_selection_executed") is False
          and ref.get("same_archive_weight_reference_not_independent_weight_validation") is True
          and ref.get("reference_must_NOT_be_called_final_regression_SNP_count") is True,
          "source outcome N or regression-eligibility flags falsely promoted")
    reject_approval(ref, ("ready_for_ldsc",))

    ld = ref.get("LD_score_source", {})
    check(ld.get("publisher_archive_md5_verified") is True
          and ld.get("publisher_archive_sha256_verified") is True
          and ld.get("reference_LD_score_chromosomes") == 22
          and ld.get("standard_M_5_50_chromosomes") == 22
          and ld.get("non_partitioned_weights_share_LD_score_file") is True,
          "22chr original LD score, normalizer or weight provenance invalid")
    chroms = ref.get("chr_reference", [])
    check(isinstance(chroms, list) and len(chroms) == 22
          and [x.get("chr") for x in chroms] == list(range(1, 23)),
          "22chr reference order/count invalid")
    check(all(type(x.get("pinned_reference_BMI_candidate_SNPs")) is int
              and x["pinned_reference_BMI_candidate_SNPs"] > 0
              and x.get("BMI_candidate_also_in_HM3_and_LD_score")
                  == x["pinned_reference_BMI_candidate_SNPs"] for x in chroms)
          and sum(x["pinned_reference_BMI_candidate_SNPs"] for x in chroms) == BMI_SOURCE_SNP,
          "original chromosome reference count mismatch")

    check(bmi_provenance.get("technical_status")
          == "PASS_CONSERVATIVE_CANDIDATE_CREATED_NOT_LDSC_APPROVED"
          and bmi_provenance.get("candidate_plaintext_sha256") == PROVENANCE_SHA
          and bmi_provenance.get("counts", {}).get("candidate_rows") == BMI_SOURCE_SNP
          and bmi_provenance.get("sources", {}).get("zenodo_archive_sha256") == REFERENCE_SHA
          and bmi_provenance.get("scientific_status") == "HOLD",
          "full-source BMI 2018 or reference checksum provenance invalid")
    reject_approval(bmi_provenance, ("ready_for_ldsc",))

    check(prereg.get("schema") == "ckm-2015-2018-BMI-and-stroke-h2-preregression-preflight-v1"
          and prereg.get("technical_status")
              == "PASS_INPUT_DESCRIPTIVE_BMI_CHI2_OUTCOME_N_GATE_ONLY"
          and prereg.get("scientific_status") == "HOLD"
          and prereg.get("LDSC_h2_or_intercept_run") is False
          and prereg.get("LDSC_rg_MR_coloc_run") is False
          and prereg.get("cross_trait_intercept_zero_fixed") is False,
          "original source descriptive chi2 or no-regression gate invalid")
    reject_approval(prereg, ("ready_for_ldsc",))
    bm = prereg.get("BMI", [])
    check(isinstance(bm, list) and [x.get("name") for x in bm] == ["locke_2015", "yengo_2018"]
          and [x.get("SNPs") for x in bm] == [LOCKE_SOURCE_SNP, BMI_SOURCE_SNP],
          "BMI source universe or distinct 2015/2018 models invalid")
    check(bm[1].get("canonical_plaintext_sha256") == PROVENANCE_SHA
          and bm[0].get("reported_N", {}).get("fractional_N_rows") == 11_506
          and bm[1].get("reported_N", {}).get("fractional_N_rows") == 0
          and all(x.get("descriptive_only", {}).get("not_an_LDSC_h2_or_intercept") is True
                  and x.get("ready_for_ldsc") is False for x in bm),
          "BMI original N rounding, hashes or descriptive-only status changed")

    outcomes = ref.get("outcomes_pre_munge_counts", [])
    nrefs = ref.get("cohort_meta_analysis_neff_illustrations", [])
    preg_out = prereg.get("stroke_outcomes", [])
    check(all(isinstance(x, list) and len(x) == 4 for x in (outcomes, nrefs, preg_out))
          and [x.get("phenotype") for x in outcomes] == list(TRAITS)
          and [x.get("phenotype") for x in nrefs] == list(TRAITS)
          and [x.get("phenotype") for x in preg_out] == list(TRAITS),
          "stroke 4 subtype missing or wrongly reordered")
    resolved = []
    for x, n, p in zip(outcomes, nrefs, preg_out):
        trait = x["phenotype"]
        accession, count = TRAITS[trait]
        check(type(x.get("valid_pre_munge_SNP_count")) is int
              and x["valid_pre_munge_SNP_count"] == count and count <= BMI_SOURCE_SNP,
              "source pre-munge SNP count drift " + trait)
        check(n.get("accession") == accession
              and p.get("accession") == accession
              and n.get("N_SNP_specific_observed") is False
              and n.get("scientific_status") == "HOLD"
              and p.get("outcome_per_variant_N_in_source") is False
              and p.get("LDSC_h2_intercept_rg_NOT_RUN") is True
              and p.get("total_first_Neff_not_approved") is True
              and p.get("cohortwise_Neff_not_approved") is True
              and p.get("individual_sample_overlap_unknown") is True
              and p.get("allow_constrain_cross_trait_intercept_to_zero") is False
              and p.get("input_gate") == "HOLD_CASE_CONTROL_META_ANALYSIS_EFFECTIVE_N_NOT_VERIFIED",
              "meta-analysis outcome N, overlap or effect gate drift " + trait)
        resolved.append({
            "phenotype": trait,
            "accession": accession,
            "pre_munge_allele_aligned_valid_SNPs": count,
            "post_munge_regression_SNPs": None,
            "effective_N": None,
            "h2": None,
            "intercept": None,
            "input_gate": "HOLD_N_EFF_AND_FINAL_REGRESSION_SNP_SET",
        })

    check(tpmi.get("schema") == "TPMI43321-crossmodel-N-LDSC-safety-gate-v1"
          and tpmi.get("technical_status")
              == "PASS_PINNED_ORIGINAL_MODELS_AND_PUBLISHED_METHODS"
          and tpmi.get("PLINK2", {}).get("source_constant_OBS_CT") == 204_099
          and tpmi.get("SAIGE", {}).get("source_N_total") == 313_909
          and tpmi.get("difference", {}).get("source_total_minus_PLINK2_OBS_CT") == 109_810
          and tpmi.get("scientific_status") == "HOLD"
          and tpmi.get("LDSC", {}).get("h2_executed") is False,
          "unverified EAS cross-model sample-size evidence changed")
    reject_approval(tpmi, ("LDSC", "ready_for_ldsc"), ("LDSC", "N_eff_approved"))

    return {
        "schema": "IS-LDSC-original-aggregates-joined-readiness-v1",
        "technical_status": "PASS_PINNED_ORIGINAL_AGGREGATE_CONSISTENCY",
        "scientific_status": "HOLD",
        "pinned_source_runs": {
            "eur_reference": 44866196133,
            "bmi_preregression": 44855617332,
            "tpmi_n_model": 45013758427,
        },
        "reference_EUR_22chr_M_5_50_integrity": "PASS_TECHNICAL",
        "reference_BMI_2018_pre_munge_SNPs": BMI_SOURCE_SNP,
        "BMI_2015_mean_chi2_descriptive": bm[0]["descriptive_only"]["mean_chi2"],
        "BMI_2018_mean_chi2_descriptive": bm[1]["descriptive_only"]["mean_chi2"],
        "EUR_stroke": resolved,
        "EAS_TPMI": {
            "SAIGE_phenotype_total": 313_909,
            "PLINK2_OBS_CT": 204_099,
            "both_model_N_not_interchangeable": True,
            "scientific_status": "HOLD",
        },
        "ldsc_h2_executed": False,
        "ldsc_rg_executed": False,
        "outcome_effective_N_approved": False,
        "post_munge_regression_SNP_count_verified": False,
        "ready_for_ldsc": False,
        "no_GWAS_rows_exported": True,
    }


def main() -> int:
    p = argparse.ArgumentParser()
    for name in ("reference", "preregression", "bmi_provenance", "tpmi", "out"):
        p.add_argument("--" + name.replace("_", "-"), type=Path, required=True)
    a = p.parse_args()
    try:
        result = validate(*[json.loads(getattr(a, name).read_text(encoding="utf-8"))
                            for name in ("reference", "preregression", "bmi_provenance", "tpmi")])
    except (OSError, KeyError, TypeError, ValueError) as exc:
        result = {
            "schema": "IS-LDSC-original-aggregates-joined-readiness-v1",
            "technical_status": "HOLD_UNVERIFIED_SOURCE_AGGREGATE",
            "scientific_status": "HOLD", "ready_for_ldsc": False,
            "ldsc_h2_executed": False, "error": str(exc)[:300],
        }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + chr(10),
                     encoding="utf-8")
    print(json.dumps({"technical_status": result["technical_status"],
                      "scientific_status": result["scientific_status"],
                      "ready_for_ldsc": result["ready_for_ldsc"],
                      "EUR_stroke": result.get("EUR_stroke"),
                      "error": result.get("error")}, ensure_ascii=False))
    return 0 if result["technical_status"] == "PASS_PINNED_ORIGINAL_AGGREGATE_CONSISTENCY" else 1


if __name__ == "__main__":
    raise SystemExit(main())
