"""Stage30 two PDC study headers synthetic shared/discordant sample groups."""
import pathlib,sys,unittest
from unittest.mock import patch
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from probe_stage30_gao_crossassembly_headers import compare,review,sample
class Stage30(unittest.TestCase):
    def test_overlapping_headers_not_join(self):
        x=compare([{"A","B","C"}],[{"A","B","Z"}])
        self.assertEqual(x["max_common_header_token_n"],2)
        self.assertAlmostEqual(x["max_header_token_jaccard"],.5)
        self.assertTrue(x["this_is_not_a_verified_same_patient_join"])
    def test_disjoint_empty_does_not_mean_no_patients(self):
        x=compare([{"A"}],[{"B"}])
        self.assertEqual(x["max_common_header_token_n"],0)
        self.assertTrue(x["this_is_not_a_verified_same_patient_join"])
    def test_safe_quant_preview_no_names(self):
        raw=b"Gene\tBiopsy1\tBiopsy2\nABO\t0.2\t0.4\nPROCR\t1.2\t3.4\n"
        class Resp:
            status=206
            def __enter__(self):return self
            def __exit__(self,*a):return None
            def read(self,n):return raw[:n]
        with patch("probe_stage30_gao_crossassembly_headers.urllib.request.urlopen",return_value=Resp()):
            x,keys=sample("https://foo.s3.amazonaws.com/file.tsv?token=secret")
        self.assertEqual(x["status"],"BOUNDED_TABULAR_HEADER")
        self.assertEqual(len(keys),2)
        self.assertNotIn("Biopsy1",str(x))
        self.assertFalse(x["sample_to_case_mapping_verified"])
    def test_strict_raw_ID_export_claim_blocked(self):
        rows=[{"pdc_study_id":p,"status":"CANDIDATE_HEADER_AUDIT_EXECUTED",
               "actual_header_preview_count_n":0,
               "source_file_or_patient_identifiers_exported":False,
               "header_previews":[]} for p in ("PDC000198","PDC000199")]
        cross=compare([],[])
        obj={"schema":"hcc-stage30-crossassembly-header-overlap-v1",
             "studies":rows,"header_set_overlap":cross}
        self.assertEqual(review(obj)["G30_D"],"STOP_NO_SAME_PATIENT_TUMOR_TISSUE_JOIN")
        cross["this_is_not_a_verified_same_patient_join"]=False
        with self.assertRaises(ValueError):review(obj)
if __name__=="__main__":unittest.main()
