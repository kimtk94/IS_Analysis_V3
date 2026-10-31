"""Synthetic source-alias join: never expose raw patient/biopsy identifiers."""
import pathlib,sys,unittest
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"/"hcc"))
from probe_ng_suter_clinical_sample_join import compare
COLS={"tumorbiopsyid":0,"patientid":1,"rnaseq":2,"proteomics":3}
def fixture():
    labels=[f"syntheticcase{k}" for k in range(50)]
    rows=[{0:f"SyntheticCase-{k}",1:f"Patient{k if k<48 else k-48}",
           2:"Y",3:"Y"} for k in range(50)]
    return labels,rows
class Tests(unittest.TestCase):
    def test_exact_fifty_biopsies_48_distinct_patients(self):
        a,b=fixture()
        r=compare(a,b,COLS)
        self.assertTrue(r["strict_same_biopsy_source_label_match"])
        self.assertEqual(r["exact_normalized_biopsy_label_matched_n"],50)
        self.assertEqual(r["distinct_nonempty_matched_patient_labels_n"],48)
        self.assertFalse(r["scientific_etiology_coding_confirmed"])
        self.assertNotIn("patient0",str(r).lower())

    def test_unknown_source_alias_remains_hold(self):
        a,b=fixture()
        b[0][0]="not_matching"
        r=compare(a,b,COLS)
        self.assertEqual(r["exact_normalized_biopsy_label_matched_n"],49)
        self.assertFalse(r["strict_same_biopsy_source_label_match"])

    def test_assay_flag_missing_not_paired(self):
        a,b=fixture()
        b[0][3]=""
        r=compare(a,b,COLS)
        self.assertEqual(r["matched_biopsies_with_both_assay_flags_n"],49)
        self.assertFalse(r["strict_same_biopsy_source_label_match"])

    def test_alias_collision_fails(self):
        a,b=fixture()
        b.append({0:"SyntheticCase-0",1:"ConflictingPatient",2:"Y",3:"Y"})
        with self.assertRaises(ValueError):compare(a,b,COLS)
if __name__=="__main__":unittest.main()
