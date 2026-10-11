import csv,unittest
from pathlib import Path
P=Path("/srv/is-analysis/results/is/audits/is_priority4_unrestricted_cross_alleles_20261011_v1/LEAD_POS_ALLELE_CROSS.tsv")
class TestCrossGWASLeadAlleles(unittest.TestCase):
 def test_all_four_positions_allele_matched(self):
  if not P.exists():self.skipTest("Frozen full-GWAS cross audit unavailable")
  with P.open() as f:r={x["site"]:x for x in csv.DictReader(f,delimiter="\t")}
  self.assertEqual(len(r),4)
  self.assertTrue(all(x["pair_match"]=="True" for x in r.values()))
  self.assertEqual(r["4:81181072"]["palindromic"],"True")
  self.assertEqual(r["4:81181072"]["aligned_beta_sign_agreement"],"PALINDROME_STRAND_UNRESOLVED")
  for k in ("4:81182554","10:105459834","10:105599770"):
   self.assertEqual(r[k]["aligned_beta_sign_agreement"],"True")
  self.assertLess(float(r["4:81181072"]["BBJ_p"]),5e-8)
  self.assertLess(float(r["4:81181072"]["GIGA_p"]),5e-8)
