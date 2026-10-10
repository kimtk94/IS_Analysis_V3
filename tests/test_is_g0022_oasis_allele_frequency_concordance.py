"""G0022 EAS-vs-OASIS AF statistical and fail-closed source lineage contract."""
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_g0022_oasis_allele_frequency_concordance import (
    metrics, read_matches, run
)

class TestOasisAfStats(unittest.TestCase):
    def test_high_same_allele_correlation(self):
        n=30;g=[.08+i*.025 for i in range(n)]
        a=[x+.015 for x in g]
        r=metrics(g,a)
        self.assertEqual(r["closer_to_same_allele_count"],n)
        self.assertEqual(r["closer_to_flipped_allele_count"],0)
        self.assertGreater(r["pearson_alt_eaf_r"],.99)

    def test_flipped_allele_detected(self):
        g=[.04+i*.021 for i in range(30)]
        a=[1-x for x in g]
        r=metrics(g,a)
        self.assertEqual(r["closer_to_flipped_allele_count"],30)
        self.assertLess(r["pearson_alt_eaf_r"],-.99)

    def test_invalid_af_fails_closed(self):
        with self.assertRaises(ValueError):metrics([.2]*11,[1.4]*11)
        with self.assertRaises(ValueError):metrics([.2]*9,[.2]*9)

    def test_bad_source_provenance_blocks_triangulation(self):
        with tempfile.TemporaryDirectory() as root:
            r=Path(root)
            origin=r/"crosswalk.json"
            origin.write_text(json.dumps({"raw_AF_greater_than_half":0,"row_by_variant_match":0}))
            class Args:
                crosswalk_audit=origin
                matched=r/"nonexistent.tsv"
                outdir=r/"out"
            with self.assertRaisesRegex(ValueError,"proof"):
                run(Args)

if __name__=="__main__":
    unittest.main()
