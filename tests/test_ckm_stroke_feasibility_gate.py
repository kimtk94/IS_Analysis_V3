"""Offline metadata-only regression tests for CKM→IS feasibility gate."""
import unittest

from scripts.ckm_stroke.feasibility_gate import assess


class FeasibilityGateTests(unittest.TestCase):
    def test_conflicting_subtype_accession_is_hold(self):
        p = [
            {"phenotype": "LAA", "ancestry": "EUR", "role": "outcome",
             "accession": "GCST90104541", "build": "GRCh37", "data_file": "REVIEW_REQUIRED"},
            {"phenotype": "CE", "ancestry": "EUR", "role": "outcome",
             "accession": "GCST90104542", "build": "GRCh37", "data_file": "REVIEW_REQUIRED"},
        ]
        observed = assess(p, [], {})
        self.assertEqual(observed["scientific_status"], "HOLD")
        self.assertEqual({r["phenotype"] for r in observed["stroke_eur_accession_conflicts"]},
                         {"LAA", "CE"})
        self.assertFalse(observed["can_run_ldsc"])

    def test_correct_official_labels_do_not_trigger_false_conflict(self):
        p = [
            {"phenotype": n, "ancestry": "EUR", "role": "outcome",
             "accession": a, "build": "GRCh37", "data_file": "REVIEW_REQUIRED"}
            for n,a in [("LAA","GCST90104542"),("CE","GCST90104541"),("SVD","GCST90104543")]
        ]
        result=assess(p,[],{})
        self.assertEqual(result["stroke_eur_accession_conflicts"],[])
        self.assertEqual(result["stroke_eur_unverified_builds"],[])
        self.assertFalse(result["can_run_ldsc"])
        self.assertEqual(result["scientific_status"],"HOLD")

    def test_url_is_not_full_gwas_qc(self):
        p = [{"phenotype": "LDL_C", "ancestry": "EUR", "role": "exposure",
              "accession": "REVIEW_REQUIRED", "build": "REVIEW_REQUIRED", "data_file": "REVIEW_REQUIRED"}]
        s = [{"trait": "LDL_C", "ancestry": "EUR", "source": "GLGC",
              "kind": "file", "url": "https://example.org/source.gz"}]
        observed = assess(p, s, {})
        self.assertIn("LDL_C", observed["eur_source_url_present_traits"])
        self.assertIn("LDL_C", observed["eur_missing_approved_normalized_input_paths"])
        self.assertFalse(observed["can_run_ldsc"])

    def test_duplicate_keys_rejected(self):
        p = [{"phenotype": "LAA", "ancestry": "EUR", "role": "outcome",
              "accession": "X", "build": "X", "data_file": "X"}]
        with self.assertRaises(ValueError):
            assess(p + p, [], {})


if __name__ == "__main__":
    unittest.main()
