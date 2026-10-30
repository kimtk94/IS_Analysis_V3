"""Synthetic privacy-preserving Xing 2023 layout diagnostic regression tests."""
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "hcc"))
from probe_xing2023_clinical_layout import (
    paired_numeric_header_keys, score_sheet, transpose_rows
)


class LayoutTests(unittest.TestCase):
    PROTEIN = ["001T", "001N", "2T", "2N"]

    def test_canonical_numeric_tn_pairs(self):
        self.assertEqual(paired_numeric_header_keys(self.PROTEIN), {"1", "2"})
        self.assertEqual(paired_numeric_header_keys(["1T", "2N", "x3T"]), set())

    def test_vertical_patient_column_candidate(self):
        rows = [(1, {1: "study"}), (2, {1: "Patient"}), (3, {1: "ID"}),
                (4, {1: "001", 2: "64"}), (5, {1: "2", 2: "46"})]
        out = score_sheet(self.PROTEIN, rows)
        self.assertEqual(out["vertical_best_numeric_overlap"], 2)
        self.assertFalse(out["clinical_field_identity_confirmed"])
        self.assertFalse(out["patient_level_join_confirmed"])
        self.assertNotIn("001", str(out))

    def test_horizontal_orientation_candidate(self):
        rows = [(1, {1: "kind", 2: "A", 3: "B"}),
                (2, {1: "ID", 2: "1", 3: "2"}),
                (3, {1: "Age", 2: "64", 3: "46"})]
        out = score_sheet(self.PROTEIN, rows)
        self.assertEqual(out["horizontal_best_numeric_overlap"], 2)
        self.assertFalse(out["patient_level_join_confirmed"])
        self.assertEqual(len(transpose_rows(rows)), 3)

    def test_matching_protein_labels_do_not_verify_clinical_join(self):
        rows = [(1, {1: "1T", 2: "1N"}), (2, {1: "2T", 2: "2N"})]
        out = score_sheet(self.PROTEIN, rows)
        self.assertEqual(out["matching_numeric_tn_header_keys_in_cells"], 2)
        self.assertFalse(out["patient_level_join_confirmed"])

    def test_empty_and_unrelated_rows(self):
        out = score_sheet(self.PROTEIN, [(4, {1: "UNRELATED", 2: "60"})])
        self.assertEqual(out["vertical_best_numeric_overlap"], 0)
        self.assertEqual(out["horizontal_best_numeric_overlap"], 0)
        self.assertEqual(out["matching_numeric_tn_header_keys_in_cells"], 0)
        self.assertFalse(out["patient_level_join_confirmed"])


if __name__ == "__main__":
    unittest.main()
