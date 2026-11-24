"""Offline contract tests for public HCC quantitative candidate manifest."""
import pathlib
import sys
import unittest
from unittest.mock import patch

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "hcc"))
import probe_quant_candidates as q


class CandidateProbeTests(unittest.TestCase):
    def test_pdc_file_index_no_patient_identifiers(self):
        with patch.object(q, "request_graphql", return_value={"getPaginatedFiles":{
            "total":2, "files":[{"file_name":"protein_assembly.tsv",
                                 "file_type":"Text","file_format":"TSV",
                                 "data_category":"Protein Assembly","file_id":"f123"}]
        }}):
            obj=q.pdc_protein_assembly()
        self.assertEqual(obj["PDC000198"]["total"],2)
        self.assertEqual(obj["PDC000199"]["files"][0]["file_name"],"protein_assembly.tsv")
        self.assertNotIn("patient_id",str(obj).lower())

    def test_pdc_query_validation_failure_is_unverified(self):
        with patch.object(q, "request_graphql", side_effect=ValueError("validation error")):
            obj=q.pdc_protein_assembly()
        self.assertTrue(all(x["status"] == "UNVERIFIED" for x in obj.values()))

    def test_pride_pagination_bounded(self):
        payload=b'[{"fileName":"P-PROT_results.xlsx","fileSizeBytes":125,"fileCategory":"RESULT"}]'
        with patch.object(q, "read_bounded", return_value=payload) as fn:
            obj=q.pride_pages(max_pages=3)
        self.assertEqual(fn.call_count,2)
        self.assertEqual(obj["PXD025705"]["candidate_count_in_scanned_pages"],1)
        self.assertFalse(obj["PXD025705"]["possibly_truncated"])

    def test_pride_fail_closed(self):
        with patch.object(q, "read_bounded", side_effect=OSError("blocked")):
            obj=q.pride_pages(max_pages=2)
        self.assertTrue(obj["PXD025705"]["errors"])
        self.assertEqual(obj["PXD025705"]["status"],"PARTIAL_METADATA_INDEX")


if __name__=="__main__":
    unittest.main()
