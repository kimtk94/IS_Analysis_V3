"""PDC experimental design source-only and privacy/claim guards."""
import pathlib
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from probe_stage30_gao_experimental_design import source,compare,verify,LABELS

class ExperimentalDesignTests(unittest.TestCase):
    def fake_query(self,query):
        if "biospecimenPerStudy" in query:
            return ({"biospecimenPerStudy":[
                {"aliquot_id":f"private_uuid_{i}",
                 "aliquot_submitter_id":f"private_aliquot_{i}"}
                for i in range(331)]},"a"*64)
        if "studyExperimentalDesign" in query:
            return ({"studyExperimentalDesign":[{
                "plex_dataset_name":"private_plex_one",
                "study_run_metadata_submitter_id":"private_run_one",
                "tmt_126":[{"aliquot_id":"private_uuid_0",
                            "aliquot_submitter_id":"private_aliquot_0"}],
                "tmt_127n":[{"aliquot_id":"private_uuid_1",
                             "aliquot_submitter_id":"private_aliquot_1"}]
            }]},"b"*64)
        raise ValueError("unexpected endpoint")
    def test_actual_source_aliquot_metadata_partial_only(self):
        with patch("probe_stage30_gao_experimental_design.pdc_query",
                   side_effect=self.fake_query):
            a,ids=source("PDC000198")
        self.assertEqual(a["status"],"OFFICIAL_PDC_EXPERIMENTAL_DESIGN_VERIFIED")
        self.assertEqual(a["total_channel_aliquot_assignments_n"],2)
        self.assertEqual(a["aliquot_IDs_also_in_official_biospecimen_n"],2)
        self.assertEqual(a["aliquot_submitter_IDs_also_in_biospecimen_n"],2)
        self.assertNotIn("private_",str(a))
        self.assertFalse(a["quant_file_column_to_channel_mapping_verified"])
        self.assertEqual(len(ids["plexes"]),1)
    def test_crosssource_design_not_RNA(self):
        with patch("probe_stage30_gao_experimental_design.pdc_query",
                   side_effect=self.fake_query):
            a,ids=source("PDC000198")
            b,ids2=source("PDC000199")
        same=compare(ids,ids2)
        self.assertEqual(same["cross_study_common_designed_aliquot_UUID_n"],2)
        self.assertFalse(same["independent_RNA_protein_join"])
        report={"schema":"hcc-stage30-gao-tmt-experimental-design-v1",
                "studies":[a,b],"cross_study_source_design":same}
        self.assertEqual(verify(report)["G30_D"],"STOP_NO_PROTEOME_RNA_PATIENT_JOIN")
        same["independent_RNA_protein_join"]=True
        with self.assertRaises(ValueError):verify(report)
    def test_transport_error_goes_access_hold_not_biological_missing(self):
        with patch("probe_stage30_gao_experimental_design.pdc_query",
                   side_effect=TimeoutError()):
            d,private=source("PDC000198")
        self.assertEqual(d["status"],"ACCESS_HOLD")
        self.assertIsNone(private)
    def test_technical_identity_counts_do_not_imply_matrix(self):
        with patch("probe_stage30_gao_experimental_design.pdc_query",
                   side_effect=self.fake_query):
            a,_=source("PDC000198")
        self.assertFalse(a["source_quant_matrix_channel_alias_verified"])
        self.assertFalse(a["RNA_same_patient_cohort_join_verified"])
        self.assertIn("tmt_126",LABELS)

if __name__=="__main__":
    unittest.main()
