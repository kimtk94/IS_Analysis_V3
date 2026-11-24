import copy
import importlib.util
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/is/audit_bbj_2021_source_freeze.py"
EVIDENCE = ROOT / "research/is/evidence/bbj2021_eas5_source_comparison_log_20261018.json"
spec = importlib.util.spec_from_file_location("bbj_freeze", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class HistoricalComparisonEvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original = json.loads(EVIDENCE.read_text(encoding="utf-8"))

    def test_historical_subset_and_ambiguity(self):
        result = module.audit_evidence(copy.deepcopy(self.original))
        self.assertEqual(result["historical_numeric_concordance_screen"], "PASS")
        self.assertEqual(result["compared_variant_keys"], 25441)
        self.assertEqual(result["effect_A2_only_under_tolerances"], 25421)
        self.assertEqual(result["effect_A1_and_A2_both_compatible_under_tolerances"], 20)
        self.assertFalse(result["new_original_data_comparison_performed"])
        self.assertFalse(result["sample_size_adopted"])
        self.assertFalse(result["coloc_executed"])
        self.assertEqual(result["scientific_status"], "HOLD")

    def test_reject_count_tamper(self):
        x = copy.deepcopy(self.original)
        x["original_log_payload"]["summary"]["same_p_with_log_tolerance"] -= 1
        with self.assertRaises(ValueError):
            module.audit_evidence(x)

    def test_reject_unverified_source_promotion(self):
        x = copy.deepcopy(self.original)
        x["original_log_payload"]["gwas_source_accepted"] = True
        with self.assertRaises(ValueError):
            module.audit_evidence(x)

    def test_reject_fake_archive_hash(self):
        x = copy.deepcopy(self.original)
        x["archive_sha256"] = "f" * 64
        with self.assertRaises(ValueError):
            module.audit_evidence(x)

    def test_reject_missing_signed_ambiguity(self):
        x = copy.deepcopy(self.original)
        x["original_log_payload"]["summary"]["hypothesis_effect_A1_agree"] = 0
        with self.assertRaises(ValueError):
            module.audit_evidence(x)

    def test_reject_nonhistorical_provenance(self):
        x = copy.deepcopy(self.original)
        x["original_run"] = "https://github.com/other/repository/actions/runs/1"
        with self.assertRaises(ValueError):
            module.audit_evidence(x)


if __name__ == "__main__":
    unittest.main()
