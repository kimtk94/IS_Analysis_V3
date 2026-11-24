#!/usr/bin/env python3
"""Synthetic mutation tests for optional MAGMA/LDSC-SEG design, not biological analyses."""
import copy
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "is"))
from validate_is_optional_enrichment import validate, summarize  # noqa: E402


class OptionalMethodContractTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ext = json.loads((ROOT / "config/research/is_optional_magma_ldscseg_v1.json").read_text(encoding="utf-8"))
        cls.core = json.loads((ROOT / "config/research/is_eur_eas_regulatory_v1.json").read_text(encoding="utf-8"))

    def test_design_valid_but_all_optional_methods_hold(self):
        self.assertEqual(validate(self.ext, self.core), [])
        statuses = summarize(self.ext, self.core)
        self.assertTrue(all(value == "HOLD_SOURCE" for item in statuses.values() for value in item.values()))

    def test_cannot_replace_primary_endpoint(self):
        e = copy.deepcopy(self.ext)
        e["primary_endpoint_unchanged"] = False
        self.assertTrue(validate(e, self.core))

    def test_cannot_claim_new_primary_hypothesis(self):
        e = copy.deepcopy(self.ext)
        e["new_primary_hypothesis"] = True
        self.assertTrue(validate(e, self.core))

    def test_wrong_EAS_gwas_accession_rejected(self):
        e = copy.deepcopy(self.ext)
        e["sources"]["GWAS_EAS_IS"]["accession"] = "GCST90104544"
        self.assertTrue(validate(e, self.core))

    def test_cross_ancestry_reference_swap_rejected(self):
        e = copy.deepcopy(self.ext)
        e["sources"]["LDSC_SEG_EAS_ANNOTATION_LDSCORES"]["ancestry"] = "EUR"
        self.assertTrue(validate(e, self.core))

    def test_missing_eas_reference_rejected(self):
        e = copy.deepcopy(self.ext)
        del e["sources"]["LDSC_SEG_EAS_ANNOTATION_LDSCORES"]
        self.assertTrue(validate(e, self.core))

    def test_no_fake_pass_without_real_analysis(self):
        e = copy.deepcopy(self.ext)
        e["analyses"][0]["eas_status"] = "PASS"
        self.assertTrue(validate(e, self.core))

    def test_evidence_must_precede_verified_sources(self):
        e = copy.deepcopy(self.ext)
        e["sources"]["LD_EUR"]["source_qc_verified"] = True
        self.assertTrue(validate(e, self.core))

    def test_magma_requires_ancestry_matched_ld(self):
        e = copy.deepcopy(self.ext)
        e["analyses"][0]["depends_on"].remove("LD_EUR")
        self.assertTrue(validate(e, self.core))

    def test_ldsc_seg_needs_baseline_model_documentation(self):
        e = copy.deepcopy(self.ext)
        e["analyses"][2]["requirements"] = ["EAS sample size", "EUR reference panel"]
        self.assertTrue(validate(e, self.core))

    def test_tissue_preselection_cannot_duplicate(self):
        e = copy.deepcopy(self.ext)
        e["prespecified_tissue_panel"]["proposed"].append(e["prespecified_tissue_panel"]["proposed"][0])
        self.assertTrue(validate(e, self.core))

    def test_h2_cutoff_not_immutable_truth(self):
        e = copy.deepcopy(self.ext)
        e["analyses"][2]["readiness_heuristic"]["interpretation"] = "This is an absolute validity criterion"
        self.assertTrue(validate(e, self.core))


if __name__ == "__main__":
    unittest.main()
