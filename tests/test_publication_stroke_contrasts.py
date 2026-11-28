import math
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from publication_stroke_contrasts import contrast,sbp_baseline,gigastroke_eur_contrasts

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
    def test_eur_sections_not_mixed(self):
        header=["Outcome","Exposure","Method","nsnp","Beta","se","P-value"]
        rows=[["Suppl. Table 24a: MR in Europeans"],header]
        for outcome,beta in [("Large artery stroke",.4),("Cardioembolic stroke",.1),("Small vessel stroke",.2)]:
            for exposure in ("LDL cholesterol","Systolic blood pressure"):
                for method in ("Inverse variance weighted","MR Egger","Weighted median"):
                    rows.append([outcome,exposure,method,10,beta,.1,math.erfc(abs(beta/.1)/math.sqrt(2))])
        rows.extend([["Suppl. Table 24b: MR in East Asians"],header,
                     ["Large artery stroke","LDL cholesterol","Inverse variance weighted",1,99,1,0]])
        class Sheet:
            def iter_rows(self,values_only=True):return iter(rows)
        out=gigastroke_eur_contrasts({"table":Sheet()},"table")
        self.assertEqual(len(out["contrasts"]),6)
        self.assertAlmostEqual(out["contrasts"][0]["primary_ivw"]["difference_log_odds"],.3)
    def test_bad_input(self):
        with self.assertRaises(ValueError):contrast(1,0,0,1)
if __name__=="__main__":unittest.main()
