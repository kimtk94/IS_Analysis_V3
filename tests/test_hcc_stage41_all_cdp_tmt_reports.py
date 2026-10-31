import pathlib,sys,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage41_all_cdp_tmt_reports import feature_kind,review,source
from audit_stage40_cdap_crossassay_source_labels import cross

class Stage41Tests(unittest.TestCase):
    def test_report_feature_type_is_candidate_only(self):
        self.assertEqual(feature_kind(["Protein Accession","X Log Ratio"]),"PROTEIN_FEATURE_SCHEMA_CANDIDATE")
        self.assertEqual(feature_kind(["Site Position","X Log Ratio"]),"PHOSPHO_FEATURE_SCHEMA_CANDIDATE")
        self.assertEqual(feature_kind(["Gene","Phosphosite"]),"MIXED_FEATURE_SCHEMA_REVIEW")
    def test_two_independent_reports_never_choose_by_position(self):
        r={"pdc_study_id":"PDC000199","status":"ALL_ORIGINAL_TMT_SOURCE_REPORT_CANDIDATES_AUDITED",
           "source_labels_or_file_names_exported":False,"official_aliquot_ids_certified":False,
           "original_PDC_quant_numeric_parity_certified":False,"RNA_to_protein_same_patient_verified":False,
           "original_TMT_report_candidate_files_n":2,
           "independent_TMT_candidate_report_checks":[
                {"status":"ORIGINAL_TMT_HEADER_VERIFIED","official_biospecimen_join_certified":False,
                 "quantDataMatrix_cells_reconciled":False,"raw_column_names_exported":False},
                {"status":"ORIGINAL_TMT_HEADER_VERIFIED","official_biospecimen_join_certified":False,
                 "quantDataMatrix_cells_reconciled":False,"raw_column_names_exported":False}]}
        other={"pdc_study_id":"PDC000198","status":"ACCESS_HOLD",
               "source_labels_or_file_names_exported":False,"official_aliquot_ids_certified":False,
               "original_PDC_quant_numeric_parity_certified":False,"RNA_to_protein_same_patient_verified":False}
        d={"schema":"hcc-stage41-CDAP-all-TMT-original-report-candidates-v1",
           "studies":[other,r],"cross_original_CDAP_sample_universe":cross(None,None)}
        self.assertEqual(review(d)["G30_D"],"STOP_NO_SOURCE_AUTHORIZED_MATCHED_RNA_PROTEIN")
        r["independent_TMT_candidate_report_checks"].pop()
        with self.assertRaises(ValueError):review(d)
    def test_aggregate_empty_if_one_source_missing(self):
        self.assertFalse(cross(None,{"biological_labels":{"A"},"analytical_samples":set()})["RNA_matching_proven"])
    def test_untrusted_patient_claim_is_not_allowed(self):
        r={"pdc_study_id":"PDC000198","status":"ACCESS_HOLD",
           "source_labels_or_file_names_exported":False,"official_aliquot_ids_certified":True,
           "original_PDC_quant_numeric_parity_certified":False,"RNA_to_protein_same_patient_verified":False}
        d={"schema":"hcc-stage41-CDAP-all-TMT-original-report-candidates-v1",
           "studies":[r,{**r,"pdc_study_id":"PDC000199","official_aliquot_ids_certified":False}],
           "cross_original_CDAP_sample_universe":cross(None,None)}
        with self.assertRaises(ValueError):review(d)

if __name__=="__main__":unittest.main()
