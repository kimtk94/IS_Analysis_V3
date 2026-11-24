"""Stage32 synthetic direct-catalog test: metadata != source-verified patient join."""
import pathlib
import sys
import unittest
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts" / "hcc"))
from audit_stage32_pdc_direct_file_categories import classify, query, verify

class Stage32Tests(unittest.TestCase):
    def test_pdc_accession_avoids_study_uuid_dependency(self):
        seen = []
        def fake(q):
            seen.append(q)
            return {"filesPerStudy": [
                {"file_name": "sample_map.tsv", "data_category": "Other Metadata", "file_format": "TSV", "md5sum": "a"*32},
                {"file_name": "proteome.gct.gz", "data_category": "Other Metadata", "file_format": "GCT", "md5sum": "b"*32},
            ]}, "c"*64
        x = query("PDC000198", "Other Metadata", fake)
        self.assertEqual(x["status"], "OFFICIAL_CATEGORY_PAGE_VERIFIED")
        self.assertEqual(x["gct_filename_candidates_n"], 1)
        self.assertEqual(x["sample_design_filename_candidates_n"], 1)
        self.assertIn('pdc_study_id:"PDC000198"', seen[0])
        self.assertNotIn("(study_id:", seen[0])
        self.assertEqual(verify({"schema":"hcc-stage32-pdc-direct-file-category-v1","evidence":x})["external_paired_RNA_protein_analysis"], "STOP")

    def test_no_premature_aliquot_certification(self):
        def fake(q): return {"filesPerStudy":[]}, "c"*64
        x = query("PDC000199", "Protein Assembly", fake)
        x["sample_to_aliquot_join_verified"] = True
        with self.assertRaises(ValueError):
            verify({"schema":"hcc-stage32-pdc-direct-file-category-v1","evidence":x})

    def test_filter_ignored_returns_hold(self):
        def fake(q):
            return {"filesPerStudy":[{"file_name":"anything.gct", "data_category":"Raw Mass Spectra"}]}, "c"*64
        x = query("PDC000198", "Other Metadata", fake)
        self.assertEqual(x["status"], "ACCESS_HOLD")
        self.assertFalse(verify({"schema":"hcc-stage32-pdc-direct-file-category-v1","evidence":x})["source_category_page_technical_success"])

    def test_error_does_not_expose_url_or_sample(self):
        def fake(q): raise ValueError("https://private-url/sample_id")
        x = query("PDC000199", "Other Metadata", fake)
        self.assertNotIn("sample_id", str(x))

    def test_filename_classification_is_only_candidate(self):
        self.assertTrue(classify("x.sample.txt")["cdap_sample_suffix"])
        self.assertFalse(classify("ordinary.txt")["gct_suffix"])

if __name__ == "__main__":
    unittest.main()
