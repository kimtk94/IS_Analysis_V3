"""HCCDB25 alternate public RNA landing accessibility guardrails."""
import pathlib
import sys
import unittest
from unittest.mock import patch
import urllib.error
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"/"hcc"))
import probe_hccdb25_access as h

class Stub:
    url="https://lifeome.net/database/hccdb2"
    status=200
    headers={"content-type":"text/html; charset=utf-8"}
    def __init__(self,contents):self.contents=contents
    def __enter__(self):return self
    def __exit__(self,*args):return False
    def read(self,n):return self.contents[:n]

class Tests(unittest.TestCase):
    def test_html_mention_is_not_expression_matrix(self):
        with patch.object(h.urllib.request,"urlopen",return_value=Stub(b"<html>HCCDB25 OEP000321 download</html>")):
            r=h.inspect_landing(h.SITES["hccdb_v2_canonical"])
        self.assertEqual(r["status"],"LANDING_HTML_ONLY")
        self.assertTrue(r["mentions_HCCDB25"])
        self.assertFalse(r["verified_hccdb25_expression_file"])
        self.assertFalse(r["joined_pdc_aliquot_ids"])
    def test_network_error_is_not_access_refusal(self):
        with patch.object(h.urllib.request,"urlopen",side_effect=urllib.error.URLError("unreachable")):
            r=h.inspect_landing(h.SITES["hccdb_v2_canonical"])
        self.assertEqual(r["status"],"NETWORK_ERROR")
        self.assertFalse(r["downloaded_rna_expression"])
    def test_canonical_sites(self):
        self.assertTrue(all("lifeome.net" in s for s in h.SITES.values()))
        self.assertEqual(h.LIMIT,163840)

if __name__=="__main__":
    unittest.main()
