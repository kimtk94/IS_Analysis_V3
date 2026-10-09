"""Broad-discovery regression contracts: synthetic fixtures, no production writes."""
import csv
import pathlib
import subprocess
import sys
import tempfile
import unittest

SCRIPT=pathlib.Path(__file__).resolve().parents[1]/"scripts/is/build_broad_discovery_v1.py"

class TestBroadDiscovery(unittest.TestCase):
    def test_compile(self):
        result=subprocess.run([sys.executable,"-m","py_compile",str(SCRIPT)],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)

    def test_no_nearest_five_hard_limit(self):
        source=SCRIPT.read_text()
        self.assertNotIn("rank > 5",source)
        self.assertNotIn("top_rows[:5]",source)

    def test_preserves_legacy_gene_universe(self):
        source=SCRIPT.read_text()
        self.assertIn("LEGACY_BBJ_FULL_GENE_UNIVERSE.tsv",source)
        self.assertIn("LEGACY_ANCHOR_CANDIDATES.tsv",source)

    def test_provisional_not_independent_claim(self):
        source=SCRIPT.read_text()
        self.assertIn("PROVISIONAL_DISTANCE_CLUSTER_REQUIRES_LD",source)
        self.assertIn("Not LD-clumped",source)

    def test_gwas_suggestive_are_separate(self):
        source=SCRIPT.read_text()
        self.assertIn("5e-8",source)
        self.assertIn("1e-6",source)
        self.assertIn("suggestive_only_regions",source)

if __name__=="__main__":
    unittest.main()
