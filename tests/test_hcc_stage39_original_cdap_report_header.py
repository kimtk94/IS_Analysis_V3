"""Stage39 exact original CDAP report labels; no inferred aliquot."""
import pathlib,sys,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage39_original_cdap_report_header import expected,compare

class Stage39Tests(unittest.TestCase):
    def test_exact_330_report_columns_candidate_not_biological_proof(self):
        ids=[f"SAMPLE_{i}" for i in range(330)]
        result=compare(["Gene"]+[x+" Log Ratio" for x in ids]+["Description"],ids)
        self.assertTrue(result["source_330_CDAP_report_header_candidate"])
        self.assertEqual(result["exact_source_CDAP_label_in_logratio_header_n"],330)
        self.assertFalse(result["official_biospecimen_join_certified"])
        self.assertFalse(result["quantDataMatrix_cells_reconciled"])
        self.assertNotIn("SAMPLE_0",str(result))
    def test_missing_two_source_labels_not_pass(self):
        ids=[f"SAMPLE_{i}" for i in range(330)]
        r=compare(["Gene"]+[x+" Log Ratio" for x in ids[2:]],ids)
        self.assertFalse(r["source_330_CDAP_report_header_candidate"])
        self.assertEqual(r["missing_expected_labels_n"],2)
    def test_reference_pool_excluded(self):
        rows=[{"AnalyticalSample":"plex_1","126":"POOL","127N":"A","127C":"B",
               "Ratios":"127N/126,127C/126"}]
        self.assertEqual(expected(rows,["126","127N","127C"]),["A","B"])
    def test_unparsed_reporter_fails_closed(self):
        with self.assertRaises(ValueError):
            expected([{"126":"POOL","Ratios":"bad_pattern"}],["126"])
    def test_colliding_biological_ids_not_330_unique(self):
        ids=["A"]*330
        r=compare(["Gene","A Log Ratio"]+["Description"],ids)
        self.assertFalse(r["source_330_CDAP_report_header_candidate"])
        self.assertEqual(r["source_distinct_biological_labels_n"],1)

if __name__=="__main__":
    unittest.main()
