"""Stage30 PDC exact file-header label vs biospecimen aliquot counter."""
import pathlib,sys,unittest
from unittest.mock import patch
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from probe_stage30_gao_quant_biospecimen_join import study,review,get_biospecimen

class OfficialPDCHarness(unittest.TestCase):
    def source(self,q):
        if "biospecimenPerStudy" in q:
            return ({"biospecimenPerStudy":[{"case_submitter_id":"CASE_PRIVATE",
                 "sample_submitter_id":"SAMPLE_PRIVATE",
                 "aliquot_submitter_id":"ALIQUOT_PRIVATE",
                 "sample_type":"Primary Tumor"}]},"a"*64)
        if "filesPerStudy" in q:
            return ({"filesPerStudy":[{"file_name":"protein_summary.tsv",
                 "data_category":"Protein Assembly","file_format":"tsv",
                 "signedUrl":[{"url":"https://foo.s3.amazonaws.com/report?token=private"}]}]},"b"*64)
        return ({"study":{"pdc_study_id":"PDC000198",
                 "study_id":"01234567-abcd-abcd-abcd-0123456789ab"}},"c"*64)
    def test_id_within_runner_only(self):
        with patch("probe_stage30_gao_quant_biospecimen_join.pdc_query",side_effect=self.source):
            with patch("probe_stage30_gao_quant_biospecimen_join.sample_header",
                       return_value=({"ALIQUOT_PRIVATE","G1"},"HEADER_READ")):
                x=study("PDC000198")
        self.assertEqual(x["status"],"SOURCE_ASSEMBLY_HEADER_VS_PDC_ALIQUOT_METADATA_CHECKED")
        self.assertEqual(x["per_candidate_aliquot_ID_overlap_counts"][0]["exact_matched_PDC_aliquot_submitter_IDs_n"],1)
        self.assertNotIn("ALIQUOT_PRIVATE",str(x))
        self.assertNotIn("token=private",str(x))
        self.assertFalse(x["full_matrix_numeric_values_and_gene_annotations_verified"])
    def test_no_exact_match_not_biological_absence(self):
        with patch("probe_stage30_gao_quant_biospecimen_join.pdc_query",side_effect=self.source):
            with patch("probe_stage30_gao_quant_biospecimen_join.sample_header",
                       return_value=({"OTHER"},"HEADER_READ")):
                x=study("PDC000198")
        self.assertEqual(x["per_candidate_aliquot_ID_overlap_counts"][0]["exact_matched_PDC_aliquot_submitter_IDs_n"],0)
        self.assertFalse(x["gene_aliquot_full_matrix_verified"])
    def test_source_status_and_no_RNA_claim(self):
        x={"schema":"hcc-stage30-quant-header-vs-pdc-aliquot-source-v1",
           "studies":[{"pdc_study_id":p,"status":"SOURCE_ASSEMBLY_HEADER_VS_PDC_ALIQUOT_METADATA_CHECKED",
               "PDC_source_unique_aliquot_n":331,"num_assembly_headers_probed_n":0,
               "per_candidate_aliquot_ID_overlap_counts":[],
               "gene_aliquot_full_matrix_verified":False,
               "signed_urls_and_raw_identifiers_exported":False,
               "RNA_case_identity_verified":False,
               "full_matrix_numeric_values_and_gene_annotations_verified":False}
                 for p in ("PDC000198","PDC000199")]}
        self.assertEqual(review(x)["G30_D"],"STOP_NO_INDEPENDENT_MATCHED_RNA_AND_PROTEIN")
        x["studies"][0]["RNA_case_identity_verified"]=True
        with self.assertRaises(ValueError):review(x)
if __name__=="__main__":unittest.main()
