"""PDC specimen-table/matrix column correspondence tests (offline)."""
import pathlib
import sys
import unittest
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"/"hcc"))
from probe_pdc_sample_mapping import schema_summary

class TestPDCMapping(unittest.TestCase):
    def test_overlap_counts_only(self):
        matrix=[["GeneSymbol","SAMPLE-ABC01","SAMPLE-DEF02"],["APOE","1","2"]]
        meta=[["Sample Name","Patient ID","Plex"],["SAMPLE-ABC01","Patient42","Plex1"],
              ["SAMPLE-DEF02","Patient43","Plex1"]]
        result=schema_summary(matrix,meta)
        self.assertTrue(result["any_exact_normalized_overlap"])
        self.assertEqual(result["top_5_metadata_fields_by_identifier_overlap"][0]["exact_normalized_header_matches"],2)
        self.assertNotIn("SAMPLE-ABC01",str(result))
        self.assertFalse(result["rna_protein_same_patient_join_verified"])

    def test_unrelated_metadata_not_false_positive(self):
        matrix=[["Gene","SAMPLE-ABC01"],["GATA3","1"]]
        meta=[["Patient"],["PersonOne"],["PersonTwo"]]
        self.assertFalse(schema_summary(matrix,meta)["any_exact_normalized_overlap"])

    def test_empty_is_unverified(self):
        self.assertEqual(schema_summary([],[]),{"status":"EMPTY_TABLE"})
if __name__=="__main__":
    unittest.main()
