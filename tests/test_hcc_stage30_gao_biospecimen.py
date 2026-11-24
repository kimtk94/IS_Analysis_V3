"""Privacy-preserving PDC case→sample→aliquot and tumor/normal source tests."""
import pathlib,sys,unittest
from unittest.mock import patch
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from probe_stage30_gao_biospecimen_overlap import records,overlap,verify,normalize_kind

class PDCBiospecimenTests(unittest.TestCase):
    def rows(self):
        return {"biospecimenPerStudy":[
            {"case_submitter_id":"SECRET_A","sample_submitter_id":"SECRET_T",
             "aliquot_submitter_id":"SECRET_T1","sample_type":"Primary Tumor"},
            {"case_submitter_id":"SECRET_A","sample_submitter_id":"SECRET_N",
             "aliquot_submitter_id":"SECRET_N1","sample_type":"Solid Tissue Normal"},
            {"case_submitter_id":"SECRET_B","sample_submitter_id":"SECRET_B_N",
             "aliquot_submitter_id":"SECRET_B1","sample_type":"Normal"}]}
    def test_direct_biospecimen_aggregation_and_ID_privacy(self):
        with patch("probe_stage30_gao_biospecimen_overlap.pdc_query",
                   return_value=(self.rows(),"a"*64)):
            summary,internal=records("PDC000198")
        self.assertEqual(summary["unique_case_n"],2)
        self.assertEqual(summary["cases_with_both_tumor_and_normal_metadata_n"],1)
        self.assertNotIn("SECRET",str(summary))
        self.assertEqual(len(internal["cases"]),2)
        self.assertFalse(summary["actual_molecular_matrix_column_samples_verified"])
    def test_crossstudy_PDC_overlap_is_not_RNA_match(self):
        with patch("probe_stage30_gao_biospecimen_overlap.pdc_query",
                   return_value=(self.rows(),"b"*64)):
            a,x=records("PDC000198");b,y=records("PDC000199")
        z=overlap(x,y)
        self.assertEqual(z["common_submitter_case_n"],2)
        self.assertEqual(z["common_submitter_sample_n"],3)
        self.assertFalse(z["RNA_case_link_verified"])
        self.assertNotIn("SECRET",str(z))
        v={"schema":"hcc-stage30-gao-proteome-phospho-biospecimen-v1",
           "studies":[a,b],"proteome_phospho_metadata_intersection":z}
        self.assertEqual(verify(v)["G30_D"],"STOP_NO_MATCHED_INDEPENDENT_RNA_PROTEOME_COHORT")
        z["RNA_case_link_verified"]=True
        with self.assertRaises(ValueError):verify(v)
    def test_sample_type_classification(self):
        self.assertEqual(normalize_kind("Primary Tumor"),"Tumor")
        self.assertEqual(normalize_kind("Solid Tissue Normal"),"Normal")
        self.assertEqual(normalize_kind(None),"OtherOrUnknown")
    def test_api_access_error_not_biological_absence(self):
        with patch("probe_stage30_gao_biospecimen_overlap.pdc_query",
                   side_effect=TimeoutError()):
            s,internal=records("PDC000199")
        self.assertEqual(s["status"],"ACCESS_HOLD")
        self.assertIsNone(internal)
if __name__=="__main__":unittest.main()
