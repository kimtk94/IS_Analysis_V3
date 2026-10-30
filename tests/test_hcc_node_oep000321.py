"""HCC NODE accession probe checks: no false claims of processed RNA availability."""
import pathlib
import sys
import unittest
from unittest.mock import patch
from io import BytesIO
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"/"hcc"))
import probe_node_oep000321 as node

class FakeResponse:
    status=200
    url=node.URL
    headers={"content-type":"text/html"}
    def __enter__(self):return self
    def __exit__(self,*args):return False
    def read(self, n):return b"<html>OEP000321</html>"

class NODEProbeTests(unittest.TestCase):
    def test_html_shell_does_not_qualify_as_expression(self):
        with patch.object(node.urllib.request,"urlopen",return_value=FakeResponse()):
            r=node.probe()
        self.assertEqual(r["status"],"PUBLIC_LANDING_HTML_ONLY")
        self.assertFalse(r["public_processed_rna_expression_matrix_url_verified"])
        self.assertIsNone(r["rna_protein_matched_patient_count"])
    def test_offline_is_not_access_denied(self):
        import urllib.error
        with patch.object(node.urllib.request,"urlopen",side_effect=urllib.error.URLError("offline")):
            r=node.probe()
        self.assertEqual(r["status"],"NETWORK_UNREACHABLE")
        self.assertFalse(r["downloaded_expression_matrix"])

if __name__=="__main__":
    unittest.main()
