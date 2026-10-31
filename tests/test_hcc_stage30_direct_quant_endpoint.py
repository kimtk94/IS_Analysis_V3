"""Stage30 PDC direct quantDataMatrix read limits and privacy guards."""
import pathlib,sys,unittest
from unittest.mock import patch
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from probe_stage30_gao_direct_quant_endpoint import review,direct,MAX_BYTES,inspect_oversized_matrix_prefix

def fixture():
    return {"schema":"hcc-stage30-gao-direct-quant-api-bounded-v1",
      "studies":[{"pdc_study_id":pid,
        "status":"RESPONSE_OVER_BYTE_LIMIT",
        "full_matrix_fetched_and_verified":False,
        "sample_or_gene_names_exported":False,
        "cell_values_exported":False,
        "download_prefix_bytes_capped":MAX_BYTES}
       for pid in ("PDC000198","PDC000199")],
      "case_and_sample_matching_performed":False}

class QuantTests(unittest.TestCase):
    def test_oversize_is_hold(self):
        x=review(fixture())
        self.assertEqual(x["studies_with_over_byte_budget_n"],2)
        self.assertTrue(x["G30_B"].startswith("HOLD"))
    def test_sample_names_cannot_be_exported(self):
        x=fixture();x["studies"][1]["sample_or_gene_names_exported"]=True
        with self.assertRaises(ValueError):review(x)
    def test_numerical_matrix_is_not_claimed(self):
        x=fixture();x["studies"][0]["full_matrix_fetched_and_verified"]=True
        with self.assertRaises(ValueError):review(x)
    def test_missing_study_rejected(self):
        x=fixture();x["studies"][0]["pdc_study_id"]="PDC000000"
        with self.assertRaises(ValueError):review(x)
    def test_oversized_real_schema_prefix(self):
        sample=(b'{"data":{"quantDataMatrix":[{"Gene/Aliquot":"PRIVATE_GENE",'
                b'"PRIVATE_SAMPLE_1":"1.23","PRIVATE_SAMPLE_2":-0.31},'
                b'{"Gene/Aliquot":"PRIVATE_GENE_2",'
                b'"PRIVATE_SAMPLE_1":0.12,"PRIVATE_SAMPLE_2":"0.4"},'
                b'{"Gene/Aliquot":"PRIVATE_GENE_3",'
                b'"PRIVATE_SAMPLE_1":"NA","PRIVATE_SAMPLE_2":1.15},')
        x=inspect_oversized_matrix_prefix(sample)
        self.assertEqual(x["first_complete_gene_records_seen_n"],3)
        self.assertEqual(x["first_record_field_count_n"],3)
        self.assertEqual(x["numeric_cells_in_complete_first_records_n"],5)
        self.assertTrue(x["sample_gene_field_present_in_partial_rows"])
        self.assertNotIn("PRIVATE_",str(x))
    def test_missing_or_malformed_prefix_has_no_leakage(self):
        x=inspect_oversized_matrix_prefix(b'{"other":[{"ID":"PRIVATE"}]}')
        self.assertFalse(x["array_prefix_marker_observed"])
        self.assertNotIn("PRIVATE",str(x))
    def test_no_unknown_study_call(self):
        with self.assertRaises(ValueError):direct("PDC000000")
    def test_source_access_timeout_is_hold(self):
        with patch("probe_stage30_gao_direct_quant_endpoint.urllib.request.urlopen",
                   side_effect=TimeoutError()):
            x=direct("PDC000198")
        self.assertEqual(x["status"],"ACCESS_HOLD")
        self.assertFalse(x["full_matrix_fetched_and_verified"])

if __name__=="__main__":unittest.main()
