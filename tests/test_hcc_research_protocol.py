"""Tests for evidence-protected HCC workflow and claim significance gates."""
import copy
import json
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "hcc"))
from validate_research_protocol import summarize, validate

PATH = ROOT / "research" / "hcc" / "validation" / "checkpoints.json"


class HCCResearchProtocolTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original = json.loads(PATH.read_text(encoding="utf-8"))

    def test_initial_status_not_overclaimed(self):
        result = summarize(self.original)
        self.assertFalse(result["primary_claim_allowed"])
        self.assertFalse(result["clinical_claim_allowed"])
        self.assertEqual(result["counts"]["PASS"], 0)
        self.assertEqual(len(self.original["checkpoints"]), 16)

    def test_duplicate_checkpoint_rejected(self):
        reg = copy.deepcopy(self.original)
        reg["checkpoints"][1]["id"] = "V00"
        with self.assertRaises(ValueError):
            validate(reg)

    def test_pass_without_review_rejected(self):
        reg = copy.deepcopy(self.original)
        reg["checkpoints"][0]["status"] = "PASS"
        with self.assertRaises(ValueError):
            validate(reg)

    def test_pass_without_real_artifact_rejected(self):
        reg = copy.deepcopy(self.original)
        reg["checkpoints"][0]["status"] = "PASS"
        reg["checkpoints"][0]["reviewer_approved"] = True
        with self.assertRaises(ValueError):
            validate(reg)

    def test_approved_artifact_can_pass_nonclinical_gate(self):
        reg = copy.deepcopy(self.original)
        reg["checkpoints"][0].update({
            "status": "PASS", "reviewer_approved": True,
            "analysis_evidence": [{
                "path": "results/hcc/audit/source_integrity.json",
                "source_sha256": "ab" * 32,
                "code_sha": "cd" * 20,
                "checked_at": "2026-10-30"
            }]
        })
        self.assertTrue(validate(reg))

    def test_paired_claim_requires_unique_patients(self):
        reg = copy.deepcopy(self.original)
        gate = reg["checkpoints"][2]
        gate.update({"status": "PASS", "reviewer_approved": True,
                     "analysis_evidence": [{
                         "path": "results/hcc/biopsy_join.json",
                         "source_sha256": "ab" * 32,
                         "code_sha": "cd" * 20,
                         "checked_at": "2026-10-30"
                     }]})
        with self.assertRaises(ValueError):
            validate(reg)

    def test_clinical_claim_requires_actual_events(self):
        reg = copy.deepcopy(self.original)
        gate = reg["checkpoints"][12]
        gate.update({"status": "PASS", "reviewer_approved": True,
                     "analysis_evidence": [{
                         "path": "results/hcc/recurrence_codebook.json",
                         "source_sha256": "ab" * 32,
                         "code_sha": "cd" * 20,
                         "checked_at": "2026-10-30"
                     }]})
        with self.assertRaises(ValueError):
            validate(reg)

    def test_fail_requires_evidence(self):
        reg = copy.deepcopy(self.original)
        reg["checkpoints"][8]["status"] = "FAIL"
        with self.assertRaises(ValueError):
            validate(reg)

    def test_clinical_gate_cannot_bypass_biological_validation(self):
        reg = copy.deepcopy(self.original)
        for i in (12, 13):
            gate = reg["checkpoints"][i]
            gate.update({
                "status": "PASS", "reviewer_approved": True,
                "independent_event_n": 34,
                "analysis_evidence": [{
                    "path": "results/hcc/clinical/event_summary.json",
                    "source_sha256": "ab" * 32,
                    "code_sha": "cd" * 20,
                    "checked_at": "2026-10-30"
                }],
            })
        summary = summarize(reg)
        self.assertFalse(summary["clinical_claim_allowed"])
        self.assertFalse(summary["primary_claim_allowed"])

    def test_invalid_analysis_hash_rejected(self):
        reg = copy.deepcopy(self.original)
        gate = reg["checkpoints"][0]
        gate.update({
            "status": "PASS", "reviewer_approved": True,
            "analysis_evidence": [{
                "path": "results/hcc/checks.json",
                "source_sha256": "unverified",
                "code_sha": "cd" * 20,
                "checked_at": "2026-10-30",
            }],
        })
        with self.assertRaises(ValueError):
            validate(reg)

    def test_unsafe_analysis_artifact_location_rejected(self):
        reg = copy.deepcopy(self.original)
        gate = reg["checkpoints"][0]
        gate.update({
            "status": "PASS", "reviewer_approved": True,
            "analysis_evidence": [{
                "path": "../tmp/private.json",
                "source_sha256": "ab" * 32,
                "code_sha": "cd" * 20,
                "checked_at": "2026-10-30",
            }],
        })
        with self.assertRaises(ValueError):
            validate(reg)

    def test_statistical_non_significance_is_not_workflow_fail(self):
        reg = copy.deepcopy(self.original)
        # Low-power/no-result is an interpretation, not decisive FAIL.
        gate = reg["checkpoints"][8]
        gate["status"] = "FAIL"
        gate["observation"] = "p=0.2 with wide CI"
        with self.assertRaises(ValueError):
            validate(reg)



if __name__ == "__main__":
    unittest.main()
