#!/usr/bin/env python3
"""Unit tests of protocol invariants; synthetic mutations only, no patient data."""
import copy
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "is"))
from validate_is_eas_study_design import validate  # noqa: E402


class StudyDesignContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.design = json.loads(
            (ROOT / "config/research/is_eur_eas_regulatory_v1.json").read_text(encoding="utf-8")
        )
        cls.registry = json.loads(
            (ROOT / "data/metadata/stroke_eas_phenotype_register_2026-11-11.json").read_text(encoding="utf-8")
        )

    def test_current_design_schema_passes_but_still_holds(self):
        self.assertEqual(validate(self.design, self.registry), [])
        self.assertTrue(all(x["status"] != "PASS" for x in self.design["decision_gates"]))
        self.assertFalse(self.design["comparison_design"]["external_validation_claim"])

    def test_any_stroke_cannot_replace_ischemic_primary(self):
        d = copy.deepcopy(self.design)
        d["primary_outcome"]["eas_gwas_accession"] = "GCST90104544"
        self.assertTrue(validate(d, self.registry))

    def test_independent_validation_claim_rejected(self):
        d = copy.deepcopy(self.design)
        d["comparison_design"]["external_validation_claim"] = True
        self.assertTrue(validate(d, self.registry))

    def test_missing_eas_eqtl_contract_rejected(self):
        d = copy.deepcopy(self.design)
        d["data_contracts"] = [x for x in d["data_contracts"] if x["id"] != "EQTL_EAS_BLOOD"]
        self.assertTrue(validate(d, self.registry))

    def test_premature_gate_pass_rejected_without_evidence(self):
        d = copy.deepcopy(self.design)
        d["decision_gates"][0]["status"] = "PASS"
        self.assertTrue(validate(d, self.registry))

    def test_treating_uncertain_case_definitions_as_verified_rejected(self):
        register = copy.deepcopy(self.registry)
        register["eas_studies"][1]["case_ascertainment"] = "CLINICAL_IMAGING_VERIFIED"
        self.assertTrue(validate(self.design, register))

    def test_non_evaluable_cannot_be_treated_as_discordant(self):
        d = copy.deepcopy(self.design)
        d["primary_endpoint"]["missingness_rule"] = "Missing EAS coverage means ancestry discordance."
        self.assertTrue(validate(d, self.registry))

    def test_no_selection_on_eas_check(self):
        d = copy.deepcopy(self.design)
        d["comparison_design"]["discovery"] = "EAS_IS_GWAS"
        self.assertTrue(validate(d, self.registry))

    def test_mismatched_phenotype_registry_rejected(self):
        register = copy.deepcopy(self.registry)
        register["eas_studies"][1]["phenotype"] = "any_stroke"
        self.assertTrue(validate(self.design, register))


if __name__ == "__main__":
    unittest.main()
