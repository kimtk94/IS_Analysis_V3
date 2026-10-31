"""Synthetic Xing etiology markers; never asserts infection from laboratory codes."""
import pathlib,sys,unittest
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"/"hcc"))
from probe_xing2023_etiology_codebook_preflight import (
    numeric_shape, classify_literal_text, summarize_etiology
)

def synthetic():
    rows={}
    for k in range(152):
        rows[f"dummy-{k}"]={
            8:"1" if k%3==0 else "0",
            13:"0" if k%4 else "1",
            16:"HBV" if k<70 else "HCV" if k<78 else
               "NAFLD" if k<80 else "Unclear"
        }
        if k<7:rows[f"dummy-{k}"][13]="?"
    return rows

class Tests(unittest.TestCase):
    def test_numeric_binary_has_no_automatic_diagnosis(self):
        v=numeric_shape([{8:"0"},{8:"1"},{8:"NA"},{8:"2"}],8)
        self.assertFalse(v["all_finite_values_binary_0_1"])
        self.assertEqual(v["numeric_zero_count"],1)
        self.assertEqual(v["numeric_one_count"],1)
        self.assertFalse(v["clinical_positive_negative_classification_verified"])

    def test_literal_labels_are_conservative(self):
        self.assertEqual(classify_literal_text("HBV"),"HBV_LITERAL")
        self.assertEqual(classify_literal_text("HBV, HCV"),"MULTIPLE_LITERAL")
        self.assertEqual(classify_literal_text("NASH"),"METABOLIC_LITERAL")
        self.assertEqual(classify_literal_text("Unknown"),"UNCLASSIFIED_OR_OTHER")
        self.assertEqual(classify_literal_text("HBsAg 12.1"),"UNCLASSIFIED_OR_OTHER")
        self.assertEqual(classify_literal_text("."),"MISSING")

    def test_152_rows_and_small_cell_suppression(self):
        r=summarize_etiology(synthetic(),{8:"HBV infection",13:"HCV infection",16:"Liver Disease"})
        self.assertEqual(r["matched_protein_to_clinical_rows"],152)
        self.assertEqual(r["joint_two_marker_numeric_coverage_n"],145)
        c=r["underlying_liver_disease_text_field"]["category_count_with_small_cells_suppressed"]
        self.assertEqual(c["HBV_LITERAL"],70)
        self.assertEqual(c["HCV_LITERAL"],8)
        self.assertIsNone(c["METABOLIC_LITERAL"])
        self.assertFalse(r["etiology_patient_group_assignment_confirmed"])
        self.assertFalse(r["etiology_interaction_power_calculable"])

    def test_no_unvalidated_source_header_accepted(self):
        with self.assertRaises(ValueError):
            summarize_etiology(synthetic(),{8:"Another",13:"HCV",16:"?"})

if __name__=="__main__":unittest.main()
