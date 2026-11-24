"""PDC Stage31 strict TMT reporter/channel source-derived quant header controls."""
import pathlib,sys,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage31_pdc_tmt_header_templates import score,review

class TestMapping(unittest.TestCase):
    def test_exact_plex_channel_matches_without_personal_ID_output(self):
        design=[{"plex":"PRIVATE_PLEX_1","run":"PRIVATE_RUN_1","channel":"tmt_127n","aliquot":"PRIVATE_BIO_1"},
                {"plex":"PRIVATE_PLEX_1","run":"PRIVATE_RUN_1","channel":"tmt_128n","aliquot":"PRIVATE_BIO_2"},
                {"plex":"PRIVATE_PLEX_1","run":"PRIVATE_RUN_1","channel":"tmt_126","aliquot":"PRIVATE_REFERENCE"}]
        hdr=["PRIVATE_PLEX_1_tmt_127n","PRIVATE_PLEX_1_tmt_128n"]
        rows=score(hdr,design,{"PRIVATE_REFERENCE"})
        a=rows[0]
        self.assertEqual(a["quant_headers_with_template_match_n"],2)
        self.assertEqual(a["distinct_unique_nonreference_aliquots_n"],2)
        self.assertNotIn("PRIVATE",str(rows))
    def test_no_overfit_from_position_or_bare_channel(self):
        design=[{"plex":"P1","run":"R1","channel":"tmt_127n","aliquot":"A1"},
                {"plex":"P2","run":"R2","channel":"tmt_127n","aliquot":"A2"}]
        rows=score(["tmt_127n","tmt_127n_EXTRA"],design,set())
        self.assertEqual(max(r["quant_headers_with_template_match_n"] for r in rows),0)
    def test_duplicate_sample_aliquot_not_counted_uniquely(self):
        design=[{"plex":"P1","run":"R1","channel":"tmt_127n","aliquot":"A1"},
                {"plex":"P2","run":"R2","channel":"tmt_127n","aliquot":"A1"}]
        rows=score(["P1_tmt_127n","P2_tmt_127n"],design,set())
        self.assertEqual(rows[0]["quant_headers_with_single_design_aliquot_n"],2)
        self.assertEqual(rows[0]["distinct_unique_nonreference_aliquots_n"],1)
        self.assertFalse(rows[0]["one_to_one_all_330_nonreference_aliquots_candidate"])
    def test_failure_not_biology_negative(self):
        s={"schema":"hcc-stage31-PDC-direct-quant-TMT-template-audit-v1",
           "studies":[{"pdc_study_id":p,"status":"ACCESS_HOLD",
              "quant_column_ids_exported":False,
              "case_or_aliquot_ids_exported":False,
              "source_file_or_gene_values_exported":False,
              "biological_match_certified":False} for p in ("PDC000198","PDC000199")]}
        x=review(s)
        self.assertEqual(x["G30_D"],"STOP_NO_MATCHED_EXTERNAL_RNA_PROTEIN_COHORT")
        s["studies"][0]["biological_match_certified"]=True
        with self.assertRaises(ValueError):review(s)
if __name__=="__main__":unittest.main()
