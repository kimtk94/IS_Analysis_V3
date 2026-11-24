"""Synthetic Stage30 PDC/NODE source data gate tests (no live access)."""
import pathlib,sys,unittest
from unittest.mock import patch
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from probe_stage30_gao_catalog import probe_study,probe_catalog,review
from probe_stage30_gao_file_catalog import probe,verify

class Stage30Tests(unittest.TestCase):
    def synthetic_study(self,pid):
        return {"study":{"pdc_study_id":pid,
          "study_id":"01234567-abcd-abcd-abcd-123456789abc",
          "study_name":"HCC Proteome","disease_type":"Hepatocellular Carcinoma",
          "experiment_type":"TMT","cases_count":159,"aliquots_count":350,
          "filesCount":[{"data_category":"Processed Data","file_type":"TSV","files_count":2}]}}
    def test_study_go_only_metadata(self):
        with patch("probe_stage30_gao_catalog.pdc_query",return_value=(self.synthetic_study("PDC000198"),"a"*64)):
            x=probe_study("PDC000198")
        self.assertEqual(x["status"],"VERIFIED_STUDY_METADATA")
        self.assertFalse(x["matched_RNA_protein_patient_count_verified"])
        self.assertFalse(x["full_quantitation_file_downloaded"])
    def test_catalog_missing_study_never_implies_no_research(self):
        with patch("probe_stage30_gao_catalog.pdc_query",return_value=({"studyCatalog":[]},"b"*64)):
            x=probe_catalog()
        self.assertEqual(x["requested_study_ids_found"],{})
    def test_transport_hold(self):
        with patch("probe_stage30_gao_catalog.pdc_query",side_effect=TimeoutError()):
            x=probe_study("PDC000199")
        self.assertEqual(x["status"],"ACCESS_HOLD")
    def test_unverified_matrix_is_not_go(self):
        x={"schema":"hcc-stage30-gao-pdc-node-source-catalog-v1",
           "pdc_studies":[{"pdc_study_id":k,"status":"VERIFIED_STUDY_METADATA",
             "full_quantitation_file_downloaded":False,
             "matched_RNA_protein_patient_count_verified":False}
              for k in ("PDC000198","PDC000199")],
           "node_project_accession":"OEP000321",
           "node":{"status":"LANDING_PAGE_RESPONDED",
             "actual_RNA_quantitation_file_verified":False,
             "source_patient_identifier_bridge_verified":False},
           "case_sample_aliquot_cross_resource_join_completed":False,
           "independent_cohort_gene_effect_replication_completed":False}
        self.assertEqual(review(x)["G30_D"],"STOP_NO_VALIDATED_MATCHED_PATIENT_JOIN")
        x["pdc_studies"][0]["matched_RNA_protein_patient_count_verified"]=True
        with self.assertRaises(ValueError):review(x)
    def test_file_preview_rejects_implied_quant_proof(self):
        s={"schema":"hcc-stage30-pdc-file-catalog-preview-v1",
            "pdc_studies":[{"pdc_study_id":pid,
                "status":"VERIFIED_PDC_FILE_PREVIEW",
                "PDC_study_uuid_not_exported":True,
                "source_file_names_and_identifiers_not_exported":True,
                "actual_gene_by_patient_quantitation_matrix_downloaded":False,
                "patient_aliquot_biospecimen_manifest_join_verified":False,
                "matrix_numerical_cells_verified":False,
                "preview_file_n":2,"total_file_records_reported":4}
                for pid in ("PDC000198","PDC000199")]}
        self.assertEqual(verify(s)["G30_B"],"HOLD_QUANT_MATRIX_NOT_DOWNLOADED")
        s["pdc_studies"][0]["matrix_numerical_cells_verified"]=True
        with self.assertRaises(ValueError):verify(s)
    def test_pdc_file_metadata_preview_no_patient_data(self):
        fake_study=self.synthetic_study("PDC000198")
        files={"getPaginatedFiles":{"total":2,"files":[
                {"file_name":"HCC_protein_abundance.tsv",
                 "data_category":"Processed Data","file_type":"TSV",
                 "file_format":"tsv","md5sum":"abc"}]}}
        with patch("probe_stage30_gao_file_catalog.pdc_query",
                   side_effect=[(fake_study,"c"*64),(files,"d"*64)]):
            v=probe("PDC000198")
        self.assertEqual(v["status"],"VERIFIED_PDC_FILE_PREVIEW")
        self.assertEqual(v["preview_filename_lexical_quant_candidate_n"],1)
        self.assertFalse(v["actual_gene_by_patient_quantitation_matrix_downloaded"])
        self.assertNotIn("HCC_protein_abundance.tsv",str(v))

if __name__=="__main__":unittest.main()
