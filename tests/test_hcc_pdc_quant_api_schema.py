"""Offline tests for bounded official PDC gene/aliquot matrix schema probe."""
import pathlib
import sys
import unittest
from unittest.mock import patch
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"/"hcc"))
import probe_pdc_quant_api_schema as probe

class TestPDCQuantAPISchema(unittest.TestCase):
    def test_scans_first_row_only(self):
        data=b'{"data":{"quantDataMatrix":[{"Gene/Aliquot":"ALB","A1":"1.3","A2":"0.2"},'+b'{"Gene/Aliquot":"APOA","BAD":"5"}]}}'
        s=probe.discover_keys(data)
        self.assertEqual(s["status"],"FIRST_RECORD_KEYS_INSPECTED")
        self.assertEqual(s["keys"],["Gene/Aliquot","A1","A2"])
        self.assertNotIn("BAD",str(s))

    def test_prefix_only_does_not_imply_complete_matrix(self):
        s=probe.discover_keys(b'{"data":{"quantDataMatrix":[{"Gene/Aliquot":"ALB","A1":')
        self.assertEqual(s["status"],"FIRST_RECORD_PREFIX_ONLY")
        self.assertIn("A1",s["keys"])

    def test_alias_count_without_identifiers(self):
        r=probe.audit_field_keys(["Gene/Aliquot","ALIQUOT1","ALIQUOT2"],[
            {"samples":[{"sample_type":"Primary Tumor","aliquots":[{"aliquot_id":"ALIQUOT1","aliquot_submitter_id":"PAT1"}]},
                        {"sample_type":"Solid Tissue Normal","aliquots":[{"aliquot_id":"ALIQUOT2"}]}]}
        ])
        self.assertEqual(r["first_record_matched_unique_aliquot_fields"],2)
        self.assertEqual(r["matched_aliquot_sample_type_distribution"]["Primary Tumor"],1)
        self.assertFalse(r["RNA_PROTEIN_SAMPLE_JOIN_VERIFIED"])
        self.assertNotIn("PAT1",str(r))

    def test_does_not_treat_api_response_error_as_matrix(self):
        result=probe.discover_keys(b'{"errors":[{"message":"Field not found"}]}')
        self.assertEqual(result["status"],"MATRIX_FIELD_ABSENT")
        self.assertEqual(result["keys"],[])

if __name__=="__main__":
    unittest.main()
