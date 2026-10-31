import pathlib,sys,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage38_targeted_source_categories import query_category,declared_category_counts,study,review,CATEGORIES

class Stage38Tests(unittest.TestCase):
    def test_exact_category_filter_and_anonymous_counts(self):
        seen=[]
        def get(q):
            seen.append(q)
            return {"filesPerStudy":[{"file_id":"one","file_name":"author-quant.tsv","file_format":"tsv",
                "data_category":"Alternate Processing Pipeline","md5sum":"f"*32}]},"sha"
        x=query_category("PDC000198","Alternate Processing Pipeline",get)
        self.assertEqual(x["status"],"CATEGORY_SOURCE_LIST_CHECKED")
        self.assertEqual(x["candidate_source_files_n"],1)
        self.assertIn('data_category:"Alternate Processing Pipeline"',seen[0])
        self.assertNotIn("author-quant.tsv",str(x))

    def test_category_mismatch_holds(self):
        x=query_category("PDC000198","Alternate Processing Pipeline",
            lambda q:({"filesPerStudy":[{"file_id":"a","data_category":"Protein Assembly"}]},"sha"))
        self.assertEqual(x["status"],"ACCESS_HOLD")
        self.assertEqual(x["error_class"],"ValueError")

    def test_official_source_counts_not_summed_without_ints(self):
        def get(q):
            return {"study":{"pdc_study_id":"PDC000199","filesCount":[
                {"data_category":"Protein Assembly","file_type":"Text","files_count":8},
                {"data_category":"Protein Assembly","file_type":"TSV","files_count":2}]}},"digest"
        counts,sha=declared_category_counts("PDC000199",get)
        self.assertEqual(counts["Protein Assembly"],10)

    def test_missing_source_count_keeps_unconfirmed(self):
        def fake(q):
            if "{study(" in q:raise ConnectionError("sensitive URL")
            return {"filesPerStudy":[]},"sha"
        d=study("PDC000198",fake)
        self.assertEqual(d["official_study_counts_status"],"ACCESS_HOLD")
        self.assertNotIn("sensitive",str(d))
        self.assertEqual(len(d["categories"]),3)

    def test_review_cannot_claim_patient_join(self):
        d={"schema":"hcc-stage38-targeted-PDC-author-category-v1","studies":[{
            "pdc_study_id":p,"patient_or_raw_source_names_exported":False,
            "data_join_complete":False,"official_study_counts_status":"ACCESS_HOLD",
            "categories":[{"category":c,"status":"ACCESS_HOLD",
              "original_file_names_urls_or_ids_exported":False,
              "same_biopsy_RNA_protein_join_verified":False} for c in CATEGORIES]
        } for p in ("PDC000198","PDC000199")]}
        self.assertEqual(review(d)["G30_D"],"STOP_NO_SAME_BIOPSY_TRANSCRIPT_PROTEIN_JOIN")
        d["studies"][0]["data_join_complete"]=True
        with self.assertRaises(ValueError):review(d)

if __name__=="__main__":unittest.main()
