"""Xing S2-proteome to S1-clinical sample-ID in-memory join tests."""
import pathlib,sys,unittest
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"/"hcc"))
from probe_xing2023_clinical_id_join import canon_id,match_clinical_ids

class Tests(unittest.TestCase):
    def test_numeric_id_canonicalization(self):
        self.assertEqual(canon_id("001"),"1")
        self.assertEqual(canon_id("152.0"),"152")
        self.assertIsNone(canon_id("AGE-45"))
        self.assertIsNone(canon_id("UP123"))
    def test_candidate_join_is_not_clinical_confirmation(self):
        prot=["1T","1N","2T","2N"]
        clinical=[(1,{1:"Header"}),(2,{1:"Patient"}),(3,{1:"ID"}),(4,{1:"1",2:"56"}),(5,{1:"2",2:"64"})]
        r=match_clinical_ids(prot,clinical)
        self.assertEqual(r["proteome_header_paired_numeric_case_keys"],2)
        self.assertEqual(r["best_clinical_column_match_count"],2)
        self.assertTrue(r["all_protein_pairs_match_clinical_candidate"])
        self.assertFalse(r["clinical_patient_column_annotated_and_confirmed"])
        self.assertFalse(r["protein_rna_quantitative_join"])
        self.assertNotIn("56",str({k:v for k,v in r.items() if k not in ("best_clinical_column_numeric_cells",)}))
    def test_false_overlap_does_not_pass(self):
        r=match_clinical_ids(["1T","1N","2T","2N"],[(4,{1:"80",2:"A1"}),(5,{1:"90",2:"A2"})])
        self.assertFalse(r["all_protein_pairs_match_clinical_candidate"])
        self.assertEqual(r["best_clinical_column_match_count"],0)

if __name__=="__main__":unittest.main()
