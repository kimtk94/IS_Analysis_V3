"""Tests for non-invasive HCCDB public JS asset discovery."""
import pathlib
import sys
import unittest
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"/"hcc"))
from probe_hccdb25_assets import Scripts,safe_url,summarize_script_text

class AssetTests(unittest.TestCase):
    def test_extract_js_from_html(self):
        p=Scripts()
        p.feed('<html><script src="/static/js/app.123.js"></script>'
               '<link rel="preload" as="script" href="/static/js/vendor.js"></html>')
        self.assertEqual(len(p.items),2)
    def test_no_foreign_urls(self):
        self.assertIsNone(safe_url("https://notlifeome.net/private"))
        self.assertEqual(safe_url("/static/app.js"),"http://lifeome.net:809/static/app.js")
    def test_mentions_do_not_promote_rna(self):
        r=summarize_script_text(b'const path="/api/download/list"; HCCDB25 download HCCDB25')
        self.assertEqual(r["token_count"]["hccdb25"],2)
        self.assertEqual(r["candidate_api_route_count"],1)
        self.assertTrue(r["does_not_verify_download"])

if __name__=="__main__":
    unittest.main()
