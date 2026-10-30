"""Public PDC TMT experimental design audit tests (synthetic identifiers only)."""
import pathlib
import sys
import unittest
from unittest.mock import patch

ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"/"hcc"))
from probe_pdc_experimental_design import flatten_design,unique_design_counts,load_design,TMT,normalized_header_fingerprint,reconstruct_from_biospecimen_runs,disambiguate_tmt_run_slots

class TestTMTDesign(unittest.TestCase):
    def test_assignments_and_reference_are_not_independent_patients(self):
        rows=[{"experiment_number":"1","plex_dataset_name":"B1",
               "tmt_126":[{"aliquot_id":"REF"}],
               "tmt_127n":[{"aliquot_id":"A1"}],
               "tmt_127c":[{"aliquot_id":"A2"}]}]
        links=[{"case":"CASE1","type":"Primary Tumor","aliquot_keys":["A1"]},
               {"case":"CASE1","type":"Solid Tissue Normal","aliquot_keys":["A2"]},
               {"case":"COMMON_REFERENCE","type":"Not Reported","aliquot_keys":["REF"]}]
        r=unique_design_counts(rows,links)
        self.assertEqual(r["number_experiment_plex_pairs"],1)
        self.assertEqual(r["number_channel_slots"],3)
        self.assertEqual(r["unique_cases_with_design_assignments"],2)
        self.assertNotIn("CASE1",str(r))
        self.assertEqual(r["ambiguous_channel_slots_multi_aliquot"],0)

    def test_duplicate_and_ambiguous_slots_detected(self):
        rows=[{"experiment_number":"2","plex_dataset_name":"B2",
               "tmt_126":[{"aliquot_id":"A1"},{"aliquot_id":"A2"},{"aliquot_id":"A1"}]}]
        links=[{"case":"C1","type":"Tumor","aliquot_keys":["A1"]},
               {"case":"C2","type":"Tumor","aliquot_keys":["A2"]}]
        r=unique_design_counts(rows,links)
        self.assertEqual(r["ambiguous_channel_slots_multi_aliquot"],1)
        self.assertEqual(r["duplicate_channel_assignments"],1)


    def test_reconstructs_run_metadata_with_known_label(self):
        records=[{"samples":[{"aliquots":[{"aliquot_id":"AL-1",
            "aliquot_run_metadata":[{"experiment_number":"E1","label":"TMT_126","alias":"P1"}]},
            {"aliquot_id":"AL-2","aliquot_run_metadata":[{"experiment_number":"E1",
                "label":"TMT_127N","alias":"P1"}]}]}]}]
        rows,summary=reconstruct_from_biospecimen_runs(records)
        self.assertEqual(summary["distinct_experiment_numbers"],1)
        self.assertEqual(summary["run_metadata_rows_with_valid_tmt_labels"],2)
        self.assertEqual(len(flatten_design(rows)),2)
        self.assertNotIn("AL-1",str(summary))


    def test_fraction_replicate_resolves_ambiguous_aliquot_channel(self):
        cases=[{"samples":[{"aliquots":[
            {"aliquot_id":"A1","aliquot_run_metadata":[{"experiment_number":"E1","label":"tmt_126","fraction":"1","replicate_number":"1"}]},
            {"aliquot_id":"A2","aliquot_run_metadata":[{"experiment_number":"E1","label":"tmt_126","fraction":"2","replicate_number":"2"}]}
        ]}]}]
        result=disambiguate_tmt_run_slots(cases)
        self.assertEqual(result["grouping_sensitivity"]["experiment_label"]["multi_aliquot_slots"],1)
        self.assertEqual(result["grouping_sensitivity"]["experiment_fraction_label"]["multi_aliquot_slots"],0)
        self.assertEqual(result["grouping_sensitivity"]["experiment_replicate_label"]["single_aliquot_slots"],2)
        self.assertNotIn("A1",str(result))

    def test_no_quantitative_tmt_run_is_not_verified(self):
        with self.assertRaises(ValueError):
            reconstruct_from_biospecimen_runs([{"samples":[]}])

    def test_public_api_bad_response_rejected(self):
        with patch("probe_pdc_experimental_design.request_graphql",return_value={"studyExperimentalDesign":[]}):
            with self.assertRaises(ValueError):
                load_design("PDC000198")

    def test_quant_header_does_not_promote_rna_join(self):
        rows=[{"experiment_number":"1","plex_dataset_name":"Plex1"}]
        with patch("probe_pdc_experimental_design.get_file_link",return_value=({},"https://pdc.cancer.gov/a")):
            with patch("probe_pdc_experimental_design.bounded_prefix",
                       return_value=(b"Gene\tLog Ratio (Plex1)\tUnshared Log Ratio (Plex1)\nA\t0.1\t0.2\n",206,True)):
                result=normalized_header_fingerprint("PDC000198",rows)
        self.assertEqual(result["quantitative_header_columns"],2)
        self.assertEqual(result["annotated_ratio_field_counts"]["unshared_log_ratio"],1)
        self.assertFalse(result["RNA_PROTEIN_JOIN_VERIFIED"])
        self.assertNotIn("Plex1",str(result))

if __name__=="__main__":
    unittest.main()
