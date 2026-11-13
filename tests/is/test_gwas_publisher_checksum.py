#!/usr/bin/env python3
"""Offline integrity test; no GWAS downloaded."""
import unittest
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"is"))
from audit_gwas_publisher_checksum import metadata_values,source_paths,INFO

class TestPublisherChecksum(unittest.TestCase):
    def test_eur_eas_ids(self):
        self.assertEqual(set(INFO),{"EUR","EAS"})
        self.assertEqual(INFO["EUR"]["accession"],"GCST90104540")
        self.assertEqual(INFO["EAS"]["accession"],"GCST90104545")
    def test_official_yaml_path(self):
        t,y=source_paths("GCST90104545")
        self.assertTrue(y.endswith("GCST90104545_buildGRCh37.tsv.gz-meta.yaml"))
        self.assertEqual(y,t+"-meta.yaml")
    def test_yaml_parse(self):
        yaml='gwas_id: "GCST90104540"\ngenome_assembly: GRCh37\ndata_file_md5sum: 1bfe4ae8a5042fb24cb0a562b05b0d2f\n'
        d=metadata_values(yaml)
        self.assertEqual(d["gwas_id"],"GCST90104540")
        self.assertEqual(d["data_file_md5sum"],"1bfe4ae8a5042fb24cb0a562b05b0d2f")
    def test_no_fake_checksum_when_metadata_missing(self):
        self.assertNotIn("data_file_md5sum",metadata_values("gwas_id: GCST90104540"))
    def test_balance_effective_n_is_not_total_n(self):
        x=INFO["EAS"]; neff=4*x["cases"]*x["controls"]/(x["cases"]+x["controls"])
        self.assertLess(neff,x["cases"]+x["controls"])
        self.assertGreater(neff,x["cases"])
if __name__=="__main__":
    unittest.main()
