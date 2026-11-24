"""Stage33 unit tests for CDAP source table semantics and research HOLD gates."""
import pathlib
import sys
import unittest
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts" / "hcc"))

from audit_stage33_original_cdap_sample_bridge import parse_sample_table, compare, review, report_candidates

class Stage33Tests(unittest.TestCase):
    def test_parses_source_reporter_columns_and_ratios(self):
        raw = (b"FileNameRegEx\tAnalyticalSample\t126\t127N\t127C\tRatios\n"
               b"^run01\tplex01\tPOOL\tAL001\tAL002\t127N/126,127C/126\n")
        rows, reporters = parse_sample_table(raw)
        self.assertEqual(reporters, ["126", "127N", "127C"])
        self.assertEqual(rows[0]["AnalyticalSample"], "plex01")
        ids = {"aliquot_id": {}, "aliquot_submitter_id": {
            "AL001": [{"aliquot":"uuid-one"}], "AL002": [{"aliquot":"uuid-two"}]}}
        x = compare(rows, reporters, ["AL001 Log Ratio", "AL002 Log Ratio"], ids)
        self.assertEqual(x["source_ratio_headers_unique_official_aliquot_n"], 2)
        self.assertEqual(x["source_ratio_headers_distinct_unique_aliquots_n"], 2)
        self.assertFalse(x["mapping_proven"])

    def test_repeated_same_header_not_certified(self):
        rows = [{"126":"POOL", "127N":"AL001", "Ratios":"127N/126"},
                {"126":"POOL", "127N":"AL001", "Ratios":"127N/126"}]
        x = compare(rows, ["126", "127N"], ["AL001"], {
            "aliquot_id":{}, "aliquot_submitter_id":{"AL001":[{"aliquot":"uuid-one"}]}})
        self.assertEqual(x["source_ratio_labels_repeated_n"], 1)
        self.assertEqual(x["source_ratio_headers_unique_official_aliquot_n"], 0)

    def test_no_unspecified_ratio_can_generate_match(self):
        rows = [{"126":"POOL","127N":"AL001","Ratios":"unknown"}]
        x = compare(rows, ["126","127N"], ["AL001"], {
            "aliquot_id":{}, "aliquot_submitter_id":{"AL001":[{"aliquot":"uuid-one"}]}})
        self.assertEqual(x["source_ratio_headers_with_direct_CDAP_label_match_n"], 0)
        self.assertEqual(x["source_ratios_unparsed_or_unassigned_n"], 1)

    def test_reject_premature_biological_claim(self):
        rows = [{"pdc_study_id":pid, "status":"ACCESS_HOLD",
                 "sample_ids_or_matrix_header_values_exported":False,
                 "signed_urls_or_file_names_exported":False,
                 "source_to_biospecimen_join_certified":False,
                 "RNA_to_protein_patient_join_certified":False} for pid in ("PDC000198","PDC000199")]
        payload = {"schema":"hcc-stage33-CDAP-original-sample-bridge-v1","studies":rows}
        self.assertEqual(review(payload)["G30_D"], "STOP_NO_MATCHED_EXTERNAL_RNA_PROTEIN_CASES")
        rows[0]["source_to_biospecimen_join_certified"]=True
        with self.assertRaises(ValueError):
            review(payload)

    def test_pool_not_sample_label(self):
        self.assertEqual(report_candidates("POOL"), set())

if __name__ == "__main__":
    unittest.main()
