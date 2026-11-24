"""Assembly category cannot substitute quantitative matrix or patient join."""
import pathlib,sys,unittest
from unittest.mock import patch
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from probe_stage30_gao_assembly_files import classify,norm_suffix,read,review

class TestAssembly(unittest.TestCase):
    def test_file_lexical_not_patient_matrix(self):
        self.assertIn("protein_quant_keyword",classify("protein_quant.tsv"))
        self.assertEqual(norm_suffix("protein_quant.tsv.gz"),".tsv.gz")
    def test_pdc_assembly_query_and_no_data_download(self):
        def api(q):
            if "filesPerStudy" in q:
                return ({"filesPerStudy":[
                    {"file_name":"study_protein_summary.tsv",
                     "data_category":"Protein Assembly","file_type":"Text",
                     "file_format":"tsv","file_size":1234,"md5sum":"a"*32}]},
                     "a"*64)
            return ({"study":{"pdc_study_id":"PDC000198",
                "study_id":"01234567-abcd-abcd-abcd-0123456789ab",
                "cases_count":171,"aliquots_count":331}}, "b"*64)
        with patch("probe_stage30_gao_assembly_files.pdc_query",side_effect=api):
            x=read("PDC000198")
        self.assertEqual(x["status"],"VERIFIED_ASSEMBLY_FILE_METADATA")
        self.assertFalse(x["quantitative_gene_by_sample_content_verified"])
        self.assertNotIn("study_protein_summary.tsv",str(x))
    def test_wrong_category_detected(self):
        def api(q):
            if "filesPerStudy" in q:
                return ({"filesPerStudy":[{"file_name":"X.raw","data_category":"Raw Mass Spectra"}]},"x"*64)
            return ({"study":{"pdc_study_id":"PDC000198",
                 "study_id":"01234567-abcd-abcd-abcd-0123456789ab"}},"x"*64)
        with patch("probe_stage30_gao_assembly_files.pdc_query",side_effect=api):
            x=read("PDC000198")
        self.assertEqual(x["status"],"ACCESS_HOLD")
    def test_assembly_content_is_not_claimed(self):
        a=[{"pdc_study_id":p,"status":"VERIFIED_ASSEMBLY_FILE_METADATA",
            "source_file_names_and_identifiers_not_exported":True,
            "signed_download_urls_not_exported":True,
            "quantitative_matrix_downloaded":False,
            "gene_aliquot_mapping_validated":False,
            "quantitative_gene_by_sample_content_verified":False,
            "assembly_file_count_in_response":8,
            "assembly_files_with_recorded_md5_n":8} for p in ("PDC000198","PDC000199")]
        x={"schema":"hcc-stage30-gao-assembly-metadata-v1","studies":a}
        self.assertEqual(review(x)["G30_D"],"STOP_NO_REAL_CASE_SAMPLE_ALIQUOT_JOIN")
        a[1]["gene_aliquot_mapping_validated"]=True
        with self.assertRaises(ValueError):review(x)

if __name__=="__main__":unittest.main()
