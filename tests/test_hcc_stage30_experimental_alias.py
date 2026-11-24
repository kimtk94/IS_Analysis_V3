"""PDC experimentalMetadata official alias-to-aliquot, no patient ID leak."""
import pathlib,sys,unittest
from unittest.mock import patch
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from probe_stage30_gao_experimental_alias import metadata,source,review

class AliasTests(unittest.TestCase):
    def fake(self,query):
        if "experimentalMetadata" in query:
            return ({"experimentalMetadata":[{
                "study_run_metadata":[{
                    "study_run_metadata_submitter_id":"PRIVATE_PLEX",
                    "aliquot_run_metadata":[
                        {"aliquot_id":"PRIVATE_A1","aliquot_submitter_id":"PRIVATE_S1",
                         "label":"TMT_127N","alias":"PRIVATE_ALIAS_1"},
                        {"aliquot_id":"PRIVATE_A2","aliquot_submitter_id":"PRIVATE_S2",
                         "label":"TMT_128N","alias":"PRIVATE_ALIAS_2"}]}]}]},"a"*64)
        if "filesPerStudy" in query:
            return ({"filesPerStudy":[{"file_name":"protein_quant.tsv",
                "file_format":"tsv","data_category":"Protein Assembly",
                "signedUrl":[{"url":"https://foo.s3.amazonaws.com/x?signature=SECRET"}]}]},"b"*64)
        if "study(" in query:
            return ({"study":{"pdc_study_id":"PDC000198",
                "study_id":"01234567-abcd-abcd-abcd-123456789abc"}},"c"*64)
        raise ValueError("unexpected source")
    def test_metadata_no_raw_alias_print(self):
        with patch("probe_stage30_gao_experimental_alias.pdc_query",side_effect=self.fake):
            a,metrics=metadata("PDC000198")
        self.assertEqual(metrics["distinct_aliquot_run_alias_n"],2)
        self.assertEqual(metrics["aliquot_run_rows_n"],2)
        self.assertNotIn("PRIVATE_",str(metrics))
        self.assertEqual(a["PRIVATE_ALIAS_1"],{"PRIVATE_A1"})
    def test_alias_match_is_not_gene_by_aliquot_proof(self):
        with patch("probe_stage30_gao_experimental_alias.pdc_query",side_effect=self.fake),\
             patch("probe_stage30_gao_experimental_alias.read_header",
                   return_value=({"status":"BOUNDED_TABULAR_HEADER"},
                                 {"PRIVATE_ALIAS_1","OTHER"})):
            x=source("PDC000198")
        self.assertEqual(x["status"],"OFFICIAL_EXPERIMENTAL_METADATA_ALIAS_AND_HEADER_CHECKED")
        self.assertEqual(x["bounded_header_alias_match"][0]["exact_experimental_metadata_alias_token_match_n"],1)
        self.assertEqual(x["bounded_header_alias_match"][0]["aliquot_UUIDs_resolved_by_exact_alias_n"],1)
        self.assertNotIn("PRIVATE_",str(x))
        self.assertNotIn("SECRET",str(x))
        self.assertFalse(x["complete_quant_numeric_matrix_verified"])
    def test_source_transport_hold_not_zero_aliquot(self):
        with patch("probe_stage30_gao_experimental_alias.pdc_query",side_effect=TimeoutError()):
            x=source("PDC000198")
        self.assertEqual(x["status"],"ACCESS_HOLD")
    def test_claim_guard(self):
        data={"schema":"hcc-stage30-experimental-alias-quant-header-v1",
              "studies":[{"pdc_study_id":p,"status":"ACCESS_HOLD",
              "source_alias_header_aliquot_identifiers_exported":False,
              "matched_RNA_protein_patients_verified":False,
              "complete_quant_numeric_matrix_verified":False}
              for p in ("PDC000198","PDC000199")]}
        self.assertEqual(review(data)["G30_D"],"STOP_NO_SAME_PATIENT_RNA_PROTEIN_JOIN")
        data["studies"][0]["matched_RNA_protein_patients_verified"]=True
        with self.assertRaises(ValueError):review(data)

if __name__=="__main__":unittest.main()
