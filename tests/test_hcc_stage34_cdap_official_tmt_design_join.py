"""Stage34 exact run/plex-to-CDAP identity matching; no order joins."""
import pathlib
import sys
import unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage34_cdap_official_tmt_design_join import analyze,channel_key,verify

def design():
    return [
        {"plex":"Plex_A","run":"Run_A","channel":"tmt_126","aliquot":"POOL_UUID"},
        {"plex":"Plex_A","run":"Run_A","channel":"tmt_127n","aliquot":"A_UUID"},
        {"plex":"Plex_A","run":"Run_A","channel":"tmt_127c","aliquot":"B_UUID"},
    ]

class Stage34Tests(unittest.TestCase):
    def test_channel_alias_is_only_reporter_normalization(self):
        self.assertEqual(channel_key("TMT_127N"),"127N")
        self.assertEqual(channel_key("127c"),"127C")
        self.assertIsNone(channel_key("random"))

    def test_analytical_sample_exact_plex_links_reporter(self):
        x=analyze([{"AnalyticalSample":"Plex_A","126":"POOL","127N":"Biological_1",
           "127C":"Biological_2","Ratios":"127N/126,127C/126"}],
           ["126","127N","127C"],design(),{"POOL_UUID"})
        self.assertEqual(x["plex_exact"]["source_reporter_ratio_positions_with_official_aliquot_n"],2)
        self.assertEqual(x["plex_exact"]["distinct_official_nonreference_aliquots_at_unique_source_positions_n"],2)
        self.assertEqual(x["run_exact"]["source_analytical_samples_uniquely_linked_to_official_run_or_plex_n"],0)

    def test_no_row_order_inference(self):
        x=analyze([{"AnalyticalSample":"Plex_NotInPDC","126":"POOL",
            "127N":"Biological_1","Ratios":"127N/126"}],
            ["126","127N"],design(),{"POOL_UUID"})
        self.assertEqual(x["plex_exact"]["source_reporter_ratio_positions_with_official_aliquot_n"],0)
        self.assertFalse(x["plex_exact"]["provisional_330_of_330_run_channel_link_candidate"])

    def test_casefold_candidate_separate_from_exact(self):
        x=analyze([{"AnalyticalSample":"plex_a","126":"POOL","127N":"Biological_1",
            "Ratios":"127N/126"}],["126","127N"],design(),{"POOL_UUID"})
        self.assertEqual(x["plex_exact"]["source_analytical_samples_uniquely_linked_to_official_run_or_plex_n"],0)
        self.assertEqual(x["plex_casefold"]["source_analytical_samples_uniquely_linked_to_official_run_or_plex_n"],1)

    def test_review_never_certifies_biology(self):
        x={"schema":"hcc-stage34-CDAP-to-PDC-TMT-source-identity-v1",
            "source":{"pdc_study_id":"PDC000198","status":"ACCESS_HOLD",
              "source_names_ids_or_urls_exported":False,
              "quant_column_to_aliquot_certified":False,
              "RNA_protein_patient_join_certified":False}}
        self.assertEqual(verify(x)["G30_D"],"STOP_NO_MATCHED_EXTERNAL_RNA_PROTEIN")
        x["source"]["quant_column_to_aliquot_certified"]=True
        with self.assertRaises(ValueError):verify(x)

if __name__=="__main__":
    unittest.main()
