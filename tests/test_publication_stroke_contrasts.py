import math
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from publication_stroke_contrasts import contrast,sbp_baseline

class ContrastTest(unittest.TestCase):
    def test_worst_se(self):
        r=contrast(.4,.05,.1,.04)
        self.assertAlmostEqual(r["max_possible_se"],.09)
        self.assertAlmostEqual(r["worst_case_p"],math.erfc((.3/.09)/math.sqrt(2)))
        self.assertTrue(all(x["p"] is None or x["p"]<=r["worst_case_p"]+1e-14 for x in r["rho_sensitivity"]))
    def test_six_family(self):
        r=sbp_baseline()
        self.assertEqual(r["family_n"],6)
        self.assertEqual(len(r["contrasts"]),3)
        self.assertTrue(r["contrasts"][0]["passes_all_covariances"])
    def test_bad_input(self):
        with self.assertRaises(ValueError):contrast(1,0,0,1)
if __name__=="__main__":unittest.main()
