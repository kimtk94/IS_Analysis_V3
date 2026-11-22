"""Offline negative cases: source-pinned LDSC preflight must never promote h2."""
import copy
import unittest

from scripts.ckm_stroke.ldsc_original_artifact_readiness import (
    BMI_SOURCE_SNP, LOCKE_SOURCE_SNP, PROVENANCE_SHA, REFERENCE_SHA,
    TRAITS, validate,
)


def originals():
    refs = [{"chr": i,
             "pinned_reference_BMI_candidate_SNPs": BMI_SOURCE_SNP - 21 if i == 1 else 1,
             "BMI_candidate_also_in_HM3_and_LD_score": BMI_SOURCE_SNP - 21 if i == 1 else 1}
            for i in range(1, 23)]
    ref = {
        "schema": "ckm-EUR-LDSC-reference-weight-and-outcome-N-source-v1",
        "technical_status": "PASS_22CHR_LD_REFERENCE_AND_SOURCE_NEFF_HOLD_AUDITED",
        "reference_BMI_candidate_count": BMI_SOURCE_SNP,
        "reference_BMI_candidate_22chr_merge_count": BMI_SOURCE_SNP,
        "reference_BMI_candidate_also_in_HM3": BMI_SOURCE_SNP,
        "scientific_status": "HOLD",
        "GC_intercept_and_h2_NOT_RUN": True,
        "rg_MR_coloc_NOT_RUN": True,
        "case_control_N_eff_per_variant_known": False,
        "outcome_post_munge_SNP_selection_executed": False,
        "same_archive_weight_reference_not_independent_weight_validation": True,
        "reference_must_NOT_be_called_final_regression_SNP_count": True,
        "ready_for_ldsc": False,
        "LD_score_source": {
            "publisher_archive_md5_verified": True,
            "publisher_archive_sha256_verified": True,
            "reference_LD_score_chromosomes": 22,
            "standard_M_5_50_chromosomes": 22,
            "non_partitioned_weights_share_LD_score_file": True,
        },
        "chr_reference": refs,
        "outcomes_pre_munge_counts": [
            {"phenotype": x, "valid_pre_munge_SNP_count": num}
            for x, (_, num) in TRAITS.items()],
        "cohort_meta_analysis_neff_illustrations": [
            {"phenotype": x, "accession": acc, "N_SNP_specific_observed": False,
             "scientific_status": "HOLD"}
            for x, (acc, _) in TRAITS.items()],
    }
    bmi_source = {
        "technical_status": "PASS_CONSERVATIVE_CANDIDATE_CREATED_NOT_LDSC_APPROVED",
        "candidate_plaintext_sha256": PROVENANCE_SHA,
        "counts": {"candidate_rows": BMI_SOURCE_SNP},
        "sources": {"zenodo_archive_sha256": REFERENCE_SHA},
        "scientific_status": "HOLD",
        "ready_for_ldsc": False,
    }
    pre = {
        "schema": "ckm-2015-2018-BMI-and-stroke-h2-preregression-preflight-v1",
        "technical_status": "PASS_INPUT_DESCRIPTIVE_BMI_CHI2_OUTCOME_N_GATE_ONLY",
        "scientific_status": "HOLD",
        "LDSC_h2_or_intercept_run": False,
        "LDSC_rg_MR_coloc_run": False,
        "cross_trait_intercept_zero_fixed": False,
        "ready_for_ldsc": False,
        "BMI": [
            {"name": "locke_2015", "SNPs": LOCKE_SOURCE_SNP,
             "reported_N": {"fractional_N_rows": 11506},
             "descriptive_only": {"mean_chi2": 1.25082047, "not_an_LDSC_h2_or_intercept": True},
             "ready_for_ldsc": False},
            {"name": "yengo_2018", "SNPs": BMI_SOURCE_SNP,
             "canonical_plaintext_sha256": PROVENANCE_SHA,
             "reported_N": {"fractional_N_rows": 0},
             "descriptive_only": {"mean_chi2": 3.98119481, "not_an_LDSC_h2_or_intercept": True},
             "ready_for_ldsc": False},
        ],
        "stroke_outcomes": [
            {"phenotype": x, "accession": acc,
             "outcome_per_variant_N_in_source": False,
             "LDSC_h2_intercept_rg_NOT_RUN": True,
             "total_first_Neff_not_approved": True,
             "cohortwise_Neff_not_approved": True,
             "individual_sample_overlap_unknown": True,
             "allow_constrain_cross_trait_intercept_to_zero": False,
             "input_gate": "HOLD_CASE_CONTROL_META_ANALYSIS_EFFECTIVE_N_NOT_VERIFIED"}
            for x, (acc, _) in TRAITS.items()],
    }
    tpmi = {
        "schema": "TPMI43321-crossmodel-N-LDSC-safety-gate-v1",
        "technical_status": "PASS_PINNED_ORIGINAL_MODELS_AND_PUBLISHED_METHODS",
        "PLINK2": {"source_constant_OBS_CT": 204099},
        "SAIGE": {"source_N_total": 313909},
        "difference": {"source_total_minus_PLINK2_OBS_CT": 109810},
        "scientific_status": "HOLD",
        "LDSC": {"h2_executed": False, "ready_for_ldsc": False, "N_eff_approved": False},
    }
    return ref, pre, bmi_source, tpmi


class JoinedArtifactGateTest(unittest.TestCase):
    def test_original_joined_aggregate_is_technical_pass_scientific_hold(self):
        r = validate(*originals())
        self.assertEqual(r["reference_BMI_2018_pre_munge_SNPs"], 1007899)
        self.assertEqual(len(r["EUR_stroke"]), 4)
        self.assertEqual(r["scientific_status"], "HOLD")
        self.assertFalse(r["ldsc_h2_executed"])
        self.assertFalse(r["ready_for_ldsc"])
        self.assertIsNone(r["EUR_stroke"][0]["effective_N"])

    def test_wrong_ce_laa_accession_rejected(self):
        r, p, b, t = originals()
        r["cohort_meta_analysis_neff_illustrations"][1]["accession"] = "GCST90104541"
        with self.assertRaisesRegex(ValueError, "outcome N"):
            validate(r, p, b, t)

    def test_one_chr_missing_is_rejected(self):
        r, p, b, t = originals()
        r["chr_reference"].pop()
        with self.assertRaisesRegex(ValueError, "22chr"):
            validate(r, p, b, t)

    def test_reference_checksum_claim_rejected(self):
        r, p, b, t = originals()
        b["sources"]["zenodo_archive_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "checksum"):
            validate(r, p, b, t)

    def test_post_munge_claim_rejected(self):
        r, p, b, t = originals()
        r["outcome_post_munge_SNP_selection_executed"] = True
        with self.assertRaisesRegex(ValueError, "falsely promoted"):
            validate(r, p, b, t)

    def test_rounded_2015_fractional_N_rejected(self):
        r, p, b, t = originals()
        p["BMI"][0]["reported_N"]["fractional_N_rows"] = 0
        with self.assertRaisesRegex(ValueError, "rounding"):
            validate(r, p, b, t)

    def test_tpmi_neff_elevated_rejected(self):
        r, p, b, t = originals()
        t["LDSC"]["N_eff_approved"] = True
        with self.assertRaisesRegex(ValueError, "unexpected LDSC approval"):
            validate(r, p, b, t)

    def test_h2_implicitly_run_rejected(self):
        r, p, b, t = originals()
        r["GC_intercept_and_h2_NOT_RUN"] = False
        with self.assertRaisesRegex(ValueError, "falsely promoted"):
            validate(r, p, b, t)

    def test_neff_sensitivity_not_per_variant_rejected(self):
        r, p, b, t = originals()
        r["cohort_meta_analysis_neff_illustrations"][0]["N_SNP_specific_observed"] = True
        with self.assertRaisesRegex(ValueError, "outcome N"):
            validate(r, p, b, t)

    def test_snp_reference_counts_do_not_define_final_LDSC(self):
        r, p, b, t = originals()
        r["reference_must_NOT_be_called_final_regression_SNP_count"] = False
        with self.assertRaisesRegex(ValueError, "falsely promoted"):
            validate(r, p, b, t)


if __name__ == "__main__":
    unittest.main()
