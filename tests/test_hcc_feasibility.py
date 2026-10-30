"""Regression tests for HCC public-metadata-only feasibility gate."""
import copy
import json
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "hcc"))
from feasibility_audit import decide, render, validate  # noqa: E402

DATA = ROOT / "research" / "hcc" / "feasibility" / "datasets.json"


class HCCFeasibilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original = json.loads(DATA.read_text(encoding="utf-8"))

    def test_real_registry_holds_without_verified_file_join(self):
        result = decide(self.original)
        self.assertEqual(result["overall_decision"], "HOLD")
        self.assertEqual(result["mechanistic_decision"], "HOLD")
        self.assertEqual(result["recurrence_prediction_decision"], "HOLD")
        self.assertEqual(len(self.original["datasets"]), 9)
        self.assertTrue(all(not d["verified_pairing"] for d in self.original["datasets"]))

    def test_publication_count_cannot_pass_file_gate(self):
        r = copy.deepcopy(self.original)
        d = r["datasets"][0]
        self.assertEqual(d["publication_paired_rna_protein_specimens"], 51)
        self.assertFalse(decide(r)["gate_results"]["file_level_rna_protein_pairs"])

    def test_verified_claim_requires_evidence(self):
        r = copy.deepcopy(self.original)
        d = r["datasets"][0]
        d["verified_pairing"] = True
        d["verified_pair_n"] = 51
        with self.assertRaises(ValueError):
            validate(r)

    def test_reject_duplicate_ids(self):
        r = copy.deepcopy(self.original)
        r["datasets"].append(copy.deepcopy(r["datasets"][0]))
        with self.assertRaises(ValueError):
            validate(r)

    def test_go_requires_proven_pairs_etiology_and_external_endpoint(self):
        r = copy.deepcopy(self.original)
        d = r["datasets"][0]
        d.update({"verified_pairing": True, "verified_pair_n": 51,
                  "verified_etiology_counts": {"alcohol": 31, "HCV": 13},
                  "evidence_artifact": "results/hcc/secure_aggregate_join.tsv"})
        r["datasets"][2].update({"verified_recurrence_events": 36,
                                 "evidence_artifact": "results/hcc/tcga_event_audit.json"})
        r["datasets"][3].update({"verified_external_recurrence_events": 24,
                                 "evidence_artifact": "results/hcc/gse14520_event_audit.json"})
        result = decide(r)
        self.assertEqual(result["overall_decision"], "GO")
        self.assertIn("alcohol", result["eligible_etiology_groups"])
        self.assertIn("HCV", result["eligible_etiology_groups"])

    def test_same_dataset_cannot_validate_itself(self):
        r = copy.deepcopy(self.original)
        d = r["datasets"][0]
        d.update({"verified_pairing": True, "verified_pair_n": 51,
                  "verified_etiology_counts": {"alcohol": 31, "HCV": 13},
                  "verified_recurrence_events": 36,
                  "verified_external_recurrence_events": 24,
                  "evidence_artifact": "results/hcc/audit.json"})
        self.assertEqual(decide(r)["recurrence_prediction_decision"], "HOLD")

    def test_render_marks_unknown_not_zero_patients(self):
        output = render(self.original, decide(self.original), [])
        self.assertIn("PUBLICATION ONLY", output)
        self.assertIn("488 microarrays", output)
        self.assertIn("Overall: HOLD", output)


if __name__ == "__main__":
    unittest.main()
