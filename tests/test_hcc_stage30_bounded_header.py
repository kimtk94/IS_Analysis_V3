"""Pure synthetic bounded sampling and no clinical/biological overclaim."""
import pathlib,sys,unittest
from unittest.mock import patch
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from probe_stage30_gao_quant_header import safe_public_https,inspect_bytes,study,review

class Gate(unittest.TestCase):
    def test_url_ssrf_guard(self):
        self.assertFalse(safe_public_https("http://pdc.cancer.gov/big.csv"))
        self.assertFalse(safe_public_https("https://pdc.cancer.gov.evil.net/data"))
        self.assertFalse(safe_public_https("https://127.0.0.1"))
        self.assertTrue(safe_public_https("https://foo.s3.amazonaws.com/abc?token=hidden"))
    def test_partial_quantified_table_is_not_patient_matrix(self):
        x=inspect_bytes(b'Gene\\tA\\tB\\nABC1\\t1.2\\t3.4\\nDEF2\\t5.6\\t7.8\\n'.replace(b"\\t",b"\t").replace(b"\\n",b"\n"))
        self.assertEqual(x["status"],"TEXT_TABLE_PREVIEWED")
        self.assertEqual(x["finite_numeric_cells_observed_n"],4)
        self.assertNotIn("ABC1",str(x))
        self.assertFalse(x["sample_by_gene_matrix_verified"])
    def test_access_api_and_header_no_patient_claim(self):
        def fake_pdc(q):
            if "filesPerStudy" in q:return ({"filesPerStudy":[{
                "file_name":"protein_summary.tsv",
                "data_category":"Protein Assembly","file_format":"tsv",
                "signedUrl":[{"url":"https://foo.s3.amazonaws.com/report.tsv?token=secret"}]}]}, "e"*64)
            return ({"study":{"pdc_study_id":"PDC000198",
                "study_id":"01234567-abcd-abcd-abcd-0123456789ab"}},"f"*64)
        with patch("probe_stage30_gao_quant_header.pdc_query",side_effect=fake_pdc):
            with patch("probe_stage30_gao_quant_header.browse_one",return_value={
                "status":"TEXT_TABLE_PREVIEWED","approx_header_columns":9,
                "finite_numeric_cells_observed_n":12,
                "sample_by_gene_matrix_verified":False,
                "header_or_sample_identifiers_exported":False}):
                x=study("PDC000198")
        self.assertEqual(x["status"],"READ_ONLY_ASSEMBLY_PREVIEW_COMPLETED")
        self.assertNotIn("secret",str(x))
        self.assertTrue(x["has_tabular_numeric_header_evidence"])
        self.assertFalse(x["matched_RNA_protein_cohort_verified"])
    def test_source_unauthorized_biological_claim_rejected(self):
        rows=[{"pdc_study_id":p,"status":"READ_ONLY_ASSEMBLY_PREVIEW_COMPLETED",
            "signed_urls_or_filename_exported":False,
            "patient_gene_or_aliquot_identifiers_exported":False,
            "full_file_downloaded":False,
            "matched_RNA_protein_cohort_verified":False,
            "gene_by_aliquot_matrix_and_tumor_normal_verified":False,
            "candidate_header_tests_n":0,
            "candidate_header_findings":[]}
            for p in ("PDC000198","PDC000199")]
        obj={"schema":"hcc-stage30-pdc-protein-assembly-bounded-header-v1",
             "studies":rows}
        self.assertEqual(review(obj)["G30_D"],"STOP_NO_SAME_PATIENT_TUMOR_RNA_PROTEOME_JOIN")
        rows[0]["matched_RNA_protein_cohort_verified"]=True
        with self.assertRaises(ValueError):review(obj)

if __name__=="__main__":unittest.main()
