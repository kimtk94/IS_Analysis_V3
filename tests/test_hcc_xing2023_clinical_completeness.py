"""Offline synthetic regression for matched clinical completeness (no real patient data)."""
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "hcc"))
from probe_xing2023_clinical_completeness import (
    classify_missing, finite_numeric, safe_label, strict_rows, inspect_availability,
    candidate_endpoint_qc
)

def mock_data():
    head = [x for k in range(1, 153) for x in (f"{k}T", f"{k}N")]
    rows = [(1, {3: "Age", 4: "Recurrence time (months)", 5: "BCLC Stage"}),
            (2, {3: "", 4: "", 5: ""})]
    for k in range(1, 153):
        row = 2 + k
        rows.append((row, {1: f"{k}T", 2: f"{k}N",
                           3: str(39 + k % 30),
                           4: "NA" if k == 1 else str(k % 25 + 1),
                           5: "A" if k % 2 else "B"}))
    return head, rows


class Tests(unittest.TestCase):
    def test_conservative_header_family_does_not_emit_raw(self):
        self.assertIn("RFS", safe_label("Recurrence time (months)"))
        self.assertIn("survival_duration", safe_label("Recurrence time (months)"))
        self.assertEqual(safe_label("T_1234"), [])
        self.assertEqual(safe_label("patient Z-12345"), [])
        self.assertTrue(classify_missing("N/A"))
        self.assertTrue(classify_missing("."))
        self.assertFalse(classify_missing("0"))
        self.assertFalse(finite_numeric("NaN"))
        self.assertTrue(finite_numeric("10.0"))

    def test_exact_152_pairs_and_column_missingness(self):
        protein, rows = mock_data()
        rep = inspect_availability(protein, rows)
        self.assertEqual(rep["matched_sample_label_pairs"], 152)
        self.assertEqual(rep["matched_rows_with_any_covar_value"], 152)
        self.assertFalse(rep["clinical_codebook_verified"])
        self.assertIsNone(rep["endpoint_event_counts"])
        self.assertFalse(rep["recurrence_or_OS_model_ready"])
        self.assertFalse(rep["rna_protein_quantitative_join"])
        age = rep["variables"][0]
        recurrence = rep["variables"][1]
        self.assertEqual((age["n_nonmissing"], age["n_numeric"]), (152, 152))
        self.assertEqual((recurrence["n_missing"], recurrence["n_numeric"]), (1, 151))
        self.assertNotIn("39", str(age))
        self.assertFalse(rep["raw_patient_identifiers_or_phenotypes_logged"])

    def test_uninterpreted_binary_event_and_duration(self):
        matched = {"a": {29: "0", 31: "3.5"},
                   "b": {29: "1", 31: "0"},
                   "c": {29: "0", 31: "NA"}}
        rep = candidate_endpoint_qc(matched, 29, 31)
        self.assertEqual(rep["event_observed_code_0_n"], 2)
        self.assertEqual(rep["event_observed_code_1_n"], 1)
        self.assertTrue(rep["event_code_set_is_binary"])
        self.assertEqual(rep["both_event_and_duration_numeric"], 2)
        self.assertEqual(rep["duration_zero_n"], 1)
        self.assertFalse(rep["actual_recurrence_events_claimed"])
        self.assertFalse(rep["survival_model_eligible"])

    def test_nonbinary_event_blocks_any_event_inference(self):
        matched = {"a": {29: "2", 31: "-2"},
                   "b": {29: "1", 31: "5"}}
        rep = candidate_endpoint_qc(matched, 29, 31)
        self.assertFalse(rep["event_code_set_is_binary"])
        self.assertEqual(rep["event_other_numeric_code_n"], 1)
        self.assertEqual(rep["duration_negative_n"], 1)
        self.assertFalse(rep["time_origin_and_censoring_confirmed"])

    def test_tn_keys_must_be_in_same_row(self):
        protein, rows = mock_data()
        altered = [(r, {**f}) for r, f in rows]
        altered[3][1][2] = "3N"
        with self.assertRaises(ValueError):
            inspect_availability(protein, altered)

    def test_duplicate_clinical_key_rejected(self):
        protein, rows = mock_data()
        with self.assertRaises(ValueError):
            strict_rows(protein, rows + [(200, {1: "1T", 2: "1N", 3: "30"})])

    def test_extra_key_fails_closed(self):
        protein, rows = mock_data()
        rows = rows[:-1]
        with self.assertRaises(ValueError):
            strict_rows(protein, rows)

if __name__ == "__main__":
    unittest.main()
