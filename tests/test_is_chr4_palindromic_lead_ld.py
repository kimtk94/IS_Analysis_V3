import csv,unittest
from pathlib import Path
P=Path("/srv/is-analysis/results/is/audits/is_chr4_palindrome_lead_ld_20261011_v1/IS_CHR4_PALINDROME_REFERENCE_LD.tsv")
class TestChr4PalindromicLD(unittest.TestCase):
 def test_strand_boundary_and_reference_corr(self):
  if not P.exists():self.skipTest("Frozen public genotype reference outputs absent")
  with P.open() as f:r={x["population"]:x for x in csv.DictReader(f,delimiter="\t")}
  self.assertEqual(set(r),{"JPT104","EAS504","otherEAS400","EUR503"})
  self.assertAlmostEqual(float(r["JPT104"]["r2"]),0.9748552298991473,places=8)
  self.assertGreater(float(r["EAS504"]["r2"]),0.94)
  self.assertEqual(int(r["EUR503"]["N"]),503)
  self.assertTrue(all(x["strand_provenance_UNRESOLVED"]=="True" for x in r.values()))
  self.assertTrue(all(x["no_conditional_GWAS_or_causal_inference"]=="True" for x in r.values()))
