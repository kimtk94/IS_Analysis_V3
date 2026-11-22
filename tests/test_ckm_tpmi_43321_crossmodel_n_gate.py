"""Fail-closed source model separation: TPMI Phecode 433.21, no LDSC inference."""
import copy
import unittest

from scripts.ckm_stroke import reconcile_tpmi_43321_crossmodel_n as m


def evidence():
    plink = {
        "schema": "TPMI43321-original-gwas-complete-row-QC-v1",
        "technical_status": "PASS_FULL_SOURCE_SHA256_AND_ALL_ROWS_PARSED_WITH_QC_FLAGS",
        "source": "TPMI_433.21_PheWeb_GRCh38_logistic_PLINK2",
        "full_original_compressed_sha256": m.PLINK_SHA,
        "pinned_prior_original_sha256_match": True,
        "rows_total_in_original": m.PLINK_ROWS,
        "value_min_max": {"min_OBS_CT": 204099, "max_OBS_CT": 204099},
        "OBS_CT_is_not_SNP_case_control_effective_N": True,
        "scientific_status": "HOLD",
        "ready_for_ldsc": False,
        "MR_coloc_LDSC_executed": False,
    }
    saige = {
        "schema": "TPMI43321-SAIGE-full-original-casecontrol-QC-v1",
        "technical_status": "PASS_FULL_GZIP_AND_SOURCE_CASE_CONTROL_SCAN",
        "source": "https://pheweb.ibms.sinica.edu.tw/download/saige_gwas/433.21",
        "computed_original_sha256": m.SAIGE_SHA,
        "records": m.SAIGE_ROWS,
        "full_gzip_crc_verified": True,
        "distinct_case_control_pairs": 1,
        "most_frequent_case_control_pairs": [{
            "N_case": 9249, "N_ctrl": 304660,
            "count": m.SAIGE_ROWS, "sum": 313909
        }],
        "all_N_pairs_identical_to_phenotype_page": True,
        "QC_counts": {
            "source_explicit_QC_PASS": 8042797,
            "source_explicit_QC_FAIL": 7545556,
        },
        "PLINK2_model_OBS_CT": 204099,
        "SAIGE_case_control_N_are_NOT_PLINK2_model_case_control_N": True,
        "SAIGE_LDSC_neff_approved": False,
        "scientific_status": "HOLD",
    }
    methods = {
        "schema": "TPMI43321-PLINK2-SAIGE-sample-basis-audit-v1",
        "TPMI_pheweb_cases": 9249,
        "TPMI_pheweb_controls": 304660,
        "TPMI_pheweb_total": 313909,
        "source_PLINK2_OBS_CT": 204099,
        "N_discrepancy": 109810,
        "regression_case_count": None,
        "regression_control_count": None,
        "per_SNP_case_control_effective_N": None,
        "scientific_status": "HOLD_NO_SOURCE_SUPPORTED_CASE_CONTROL_BREAKDOWN",
        "pheweb": {
            "status": "PASS_PHEWEB_LINKS",
            "data_links": [
                {"href": "/download/gwas_plink2/433.21"},
                {"href": "/download/saige_gwas/433.21"}
            ]
        },
        "primary_papers": {
            "Chen_TPMI_PRS_2025": {
                "status": "PASS_PRIMARY_XML",
                "source": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12675292/fullTextXML",
                "method_evidence": {"SAIGE": {"excerpts": [
                    "The mixed-effect model used SAIGE. In the GWAS set, "
                    "we selected an unrelated subset (n = 248,754) to perform GWAS "
                    "using a generalized linear model with PLINK2, and conducted "
                    "1:10 age, sex-matching for the imbalanced traits."
                ]}}
            }
        }
    }
    return plink, saige, methods


class SourceModelGateTests(unittest.TestCase):
    def test_two_original_source_models_pass_technical_not_scientific(self):
        out = m.validate(*evidence())
        self.assertEqual(out["technical_status"],
                         "PASS_PINNED_ORIGINAL_MODELS_AND_PUBLISHED_METHODS")
        self.assertEqual(out["difference"]["source_total_minus_PLINK2_OBS_CT"], 109810)
        self.assertTrue(out["difference"]["model_design_difference_supported_by_paper"])
        self.assertFalse(out["difference"]["exact_43321_selection_and_subset_composition_verified"])
        self.assertFalse(out["LDSC"]["ready_for_ldsc"])
        self.assertFalse(out["LDSC"]["N_eff_approved"])
        self.assertEqual(out["scientific_status"], "HOLD")

    def test_saige_and_plink_model_swap_rejected(self):
        p, s, d = evidence()
        p["source"] = s["source"]
        with self.assertRaisesRegex(ValueError, "PLINK2"):
            m.validate(p, s, d)

    def test_saige_count_drift_rejected(self):
        p, s, d = evidence()
        s["most_frequent_case_control_pairs"][0]["N_case"] += 1
        with self.assertRaisesRegex(ValueError, "SAIGE"):
            m.validate(p, s, d)

    def test_saige_qc_fail_cannot_be_treated_as_pass(self):
        p, s, d = evidence()
        s["QC_counts"]["source_explicit_QC_PASS"] = m.SAIGE_ROWS
        with self.assertRaisesRegex(ValueError, "SAIGE"):
            m.validate(p, s, d)

    def test_source_sha_drift_rejected(self):
        p, s, d = evidence()
        s["computed_original_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "SAIGE"):
            m.validate(p, s, d)

    def test_original_phenotype_specific_plink_case_counts_unknown(self):
        p, s, d = evidence()
        d["regression_case_count"] = 10000
        with self.assertRaisesRegex(ValueError, "unresolved N"):
            m.validate(p, s, d)

    def test_methods_source_missing_rejected(self):
        p, s, d = evidence()
        d["primary_papers"]["Chen_TPMI_PRS_2025"]["method_evidence"]["SAIGE"]["excerpts"] = []
        with self.assertRaisesRegex(ValueError, "primary-paper"):
            m.validate(p, s, d)

    def test_plink_obs_not_approved_as_binary_neff(self):
        p, s, d = evidence()
        p["OBS_CT_is_not_SNP_case_control_effective_N"] = False
        with self.assertRaisesRegex(ValueError, "N_eff"):
            m.validate(p, s, d)

    def test_pheweb_routes_must_remain_distinct(self):
        p, s, d = evidence()
        d["pheweb"]["data_links"].pop()
        with self.assertRaisesRegex(ValueError, "PheWeb"):
            m.validate(p, s, d)

    def test_published_ratio_not_proof_of_exact_case_selection(self):
        p, s, d = evidence()
        d = copy.deepcopy(d)
        snippets = d["primary_papers"]["Chen_TPMI_PRS_2025"]["method_evidence"]["SAIGE"]["excerpts"]
        snippets[0] = snippets[0].replace("1:10", "2:10")
        with self.assertRaisesRegex(ValueError, "published"):
            m.validate(p, s, d)


if __name__ == "__main__":
    unittest.main()
