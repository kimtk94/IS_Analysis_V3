"""Stage40 strict source names; source names do NOT establish patient identity."""
import pathlib,sys,unittest,urllib.error
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage40_cdap_crossassay_source_labels import cross,retry,review

class Stage40Tests(unittest.TestCase):
    @staticmethod
    def labels(prefix="S",n=330):
        return {"biological_labels":{f"{prefix}{i:04}" for i in range(n)},
                "analytical_samples":{f"plex_{j:02}" for j in range(33)}}
    def test_full_exact_overlap_is_not_patient_join(self):
        a,b=self.labels(),self.labels()
        d=cross(a,b)
        self.assertTrue(d["source_same_330_label_universe_candidate"])
        self.assertEqual(d["source_exact_matching_biological_names_n"],330)
        self.assertFalse(d["official_PDC_aliquot_or_case_join_certified"])
        self.assertFalse(d["RNA_matching_proven"])
        self.assertNotIn("S0001",str(d))
    def test_disjoint_labels_does_not_infer_patient(self):
        d=cross(self.labels("P"),self.labels("Q"))
        self.assertEqual(d["source_exact_matching_biological_names_n"],0)
        self.assertFalse(d["source_same_330_label_universe_candidate"])
    def test_partial_overlap_is_counted_without_names(self):
        b=self.labels();b["biological_labels"].remove("S0001");b["biological_labels"].add("DIFFERENT")
        d=cross(self.labels(),b)
        self.assertEqual(d["source_exact_matching_biological_names_n"],329)
        self.assertNotIn("DIFFERENT",str(d))
    def test_casefold_is_reported_not_upgraded_to_exact(self):
        b=self.labels()
        b["biological_labels"]={x.lower() for x in b["biological_labels"]}
        d=cross(self.labels(),b)
        self.assertEqual(d["source_exact_matching_biological_names_n"],0)
        self.assertEqual(d["source_casefold_matching_biological_names_n"],330)
    def test_missing_source_avoids_false_empty_cohort(self):
        self.assertEqual(cross(self.labels(),None)["status"],"SOURCE_CROSS_ASSAY_HOLD")
    def test_retry_limited_transient_http_503(self):
        calls=[]
        def fn():
            calls.append(1)
            if len(calls)==1:raise urllib.error.HTTPError("https://pdc.cancer.gov",503,"Busy",{},None)
            return "data"
        self.assertEqual(retry(fn,sleep=lambda x:None),"data")
        self.assertEqual(len(calls),2)
    def test_retry_does_not_mask_404(self):
        def fn():raise urllib.error.HTTPError("https://pdc.cancer.gov",404,"Missing",{},None)
        with self.assertRaises(urllib.error.HTTPError):retry(fn,sleep=lambda x:None)
    def test_safety_review_rejects_patient_claim(self):
        r={"schema":"hcc-stage40-Gao-original-CDAP-crossassay-names-v1",
           "studies":[{"pdc_study_id":x,"status":"ACCESS_HOLD",
             "source_names_or_urls_exported":False,"cross_assay_patient_identity_certified":False,
             "PDC_quant_matrix_to_aliquot_certified":False,"external_same_patient_RNA_protein_proven":False}
             for x in ("PDC000198","PDC000199")],
           "cross_assay":cross(None,None)}
        self.assertEqual(review(r)["G30_D"],"STOP_NO_SOURCE_AUTHORIZED_RNA_PROTEIN_SAME_PATIENT_JOIN")
        r["studies"][0]["cross_assay_patient_identity_certified"]=True
        with self.assertRaises(ValueError):review(r)

if __name__=="__main__":
    unittest.main()
