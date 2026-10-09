"""Numerical sanity checks for EAS reference r² estimation."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from build_anchor_ld_proxies import correlation

class TestAnchorLD(unittest.TestCase):
    def test_perfect_and_reverse_alleles(self):
        a=np.tile(np.array([0.,1.,2.]),20)
        self.assertAlmostEqual(correlation(a,a)[0],1.0,places=12)
        self.assertAlmostEqual(correlation(a,2-a)[0],1.0,places=12)
        self.assertEqual(correlation(a,a)[1],60)
    def test_zero_variance_is_not_evidence(self):
        a=np.tile(np.array([0.,1.,2.]),20)
        self.assertIsNone(correlation(a,np.ones(60))[0])
    def test_missing_sample_pairwise(self):
        a=np.tile(np.array([0.,1.,2.]),20)
        b=2-a
        b[:7]=np.nan
        self.assertEqual(correlation(a,b)[1],53)
        self.assertAlmostEqual(correlation(a,b)[0],1,places=12)
        self.assertIsNone(correlation(a,np.full(60,np.nan))[0])
if __name__=="__main__":
    unittest.main()
