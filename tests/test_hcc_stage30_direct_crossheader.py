"""Direct PDC quant array source cross-header overlap, privacy and inference gates."""
import pathlib,sys,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from probe_stage30_gao_direct_quant_crossheader import compare,review

class TestCrossHeader(unittest.TestCase):
    def test_exact_common_alias_counts(self):
        d=compare(["PRIVATE_A","PRIVATE_B","PRIVATE_C"],
                  ["PRIVATE_C","PRIVATE_B","PRIVATE_D"])
        self.assertEqual(d["cross_study_common_source_header_token_n"],2)
        self.assertEqual(d["identical_source_header_order_position_n"],1)
        self.assertEqual(d["proteome_only_header_token_n"],1)
        self.assertNotIn("PRIVATE_",str(d))
    def test_one_study_no_proof(self):
        x={"schema":"hcc-stage30-gao-direct-quant-crossheader-v1",
           "studies":{pid:{"status":"SOURCE_HEADER_ACCESS_HOLD",
                 "aliquot_header_labels_exported":False} for pid in ("PDC000198","PDC000199")},
           "comparison":{"status":"CROSS_HEADER_NOT_VERIFIED",
                 "patient_header_alias_identifiers_exported":False},
           "official_patient_level_ID_pairing_proven":False,
           "matched_RNA_protein_matrix_available":False,
           "independent_biological_replication_completed":False}
        self.assertEqual(review(x)["G30_D"],"STOP_NO_VALIDATED_CASE_RNA_PROTEOME_JOIN")
    def test_no_export(self):
        d=compare(["A","B"],["B","A"])
        self.assertFalse(d["patient_header_alias_identifiers_exported"])
        with self.assertRaises(ValueError):compare(["A","A"],["A","B"])
if __name__=="__main__":unittest.main()
