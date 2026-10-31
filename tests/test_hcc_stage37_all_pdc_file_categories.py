import pathlib
import sys
import unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage37_all_pdc_file_categories import source, query_page, filename_flags, verify

class Stage37Tests(unittest.TestCase):
    def test_no_category_filter_in_query_and_no_source_names_in_output(self):
        queries=[]
        def f(q):
            queries.append(q)
            return {"filesPerStudy":[
                {"file_id":"u1","file_name":"author_quant.gct.gz",
                 "data_category":"Alternate Processing Pipeline",
                 "file_format":"GCT","md5sum":"f"*32},
                {"file_id":"u2","file_name":"A.sample.txt",
                 "data_category":"Protein Assembly","file_format":"txt",
                 "md5sum":"f"*32}
            ]},"source_sha"
        x=source("PDC000198",fetch=f)
        self.assertEqual(x["status"],"PDC_ALL_FILE_CATEGORIES_ENUMERATED")
        self.assertEqual(x["alternate_processing_pipeline_files_n"],1)
        self.assertEqual(x["original_PDC_GCT_file_catalog_candidate_n"],1)
        self.assertTrue(x["all_pages_exhausted"])
        self.assertNotIn("data_category:",queries[0])
        self.assertNotIn("author_quant",str(x))
        self.assertFalse(x["quant_330_aliquot_join_certified"])

    def test_empty_first_page_is_not_proof_of_absence(self):
        x=source("PDC000199",lambda q:({"filesPerStudy":[]},"sha"))
        self.assertEqual(x["status"],"ACCESS_HOLD")

    def test_duplicate_file_identity_is_access_hold(self):
        def f(q):
            return {"filesPerStudy":[{"file_id":"a","file_name":"a.txt","data_category":"X"},
                 {"file_id":"a","file_name":"b.txt","data_category":"X"}]},"sha"
        self.assertEqual(source("PDC000198",f)["status"],"ACCESS_HOLD")

    def test_truncated_max_pages_never_claims_exhaustive(self):
        def f(q):
            return {"filesPerStudy":[{"file_id":f"f{q.split('offset:')[1].split(',')[0]}_{i}",
                "file_name":"file.tsv","data_category":"X"} for i in range(75)]},"sha"
        x=source("PDC000198",f)
        self.assertEqual(x["status"],"TRUNCATED_MAX_PAGES_REVIEW")
        self.assertFalse(x["all_pages_exhausted"])

    def test_filename_classification_is_only_candidate(self):
        r=filename_flags("study_sample_map.tsv")
        self.assertTrue(r["author_design"])
        self.assertFalse(r["gct"])

    def test_scientific_gate_cannot_be_bypassed(self):
        v={"schema":"hcc-stage37-all-PDC-source-file-categories-v1",
          "studies":[{"pdc_study_id":p,"status":"ACCESS_HOLD",
             "patient_aliquot_or_filename_exported":False,
             "signed_url_or_raw_quant_exported":False,
             "quant_330_aliquot_join_certified":False,
             "RNA_protein_pairing_certified":False} for p in ("PDC000198","PDC000199")]}
        self.assertEqual(verify(v)["G30_D"],"STOP_NO_INDEPENDENT_SAME_BIOPSY_RNA_PROTEIN")
        v["studies"][0]["quant_330_aliquot_join_certified"]=True
        with self.assertRaises(ValueError):verify(v)

if __name__=="__main__":unittest.main()
