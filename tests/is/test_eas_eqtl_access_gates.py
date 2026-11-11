#!/usr/bin/env python3
"""Offline eQTL source identity tests, never equate source HEAD with coloc."""
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"is"))
from audit_eas_eqtl_access_gates import URLS,query_head

class TestEqtls(unittest.TestCase):
    def test_eas_record_accession(self):
        self.assertEqual(URLS["MAEEA_ZENODO"],"https://zenodo.org/api/records/21296030")
    def test_japanese_healthy_file(self):
        self.assertIn("hum0099.v1.eqtl.v1.zip",URLS["JAPAN_HEALTHY_IMMUNE_EQTL_ARCHIVE"])
    def test_japanese_jctf_readme(self):
        self.assertIn("hum0343",URLS["JAPAN_JCTF_EQTL_README"])
    def test_eur_full_association_not_significant_only(self):
        self.assertTrue(URLS["EUR_BLUEPRINT_MONOCYTE_ALL_ASSOC"].endswith("QTD000021.all.tsv.gz"))
    def test_no_fake_restricted_or_sample_meaning(self):
        self.assertIn("gtex-resources",URLS["GTEX_V8_WHOLE_BLOOD_ALLPAIRS"])
        self.assertNotIn("token",repr(URLS).lower())
if __name__=="__main__":
    unittest.main()
