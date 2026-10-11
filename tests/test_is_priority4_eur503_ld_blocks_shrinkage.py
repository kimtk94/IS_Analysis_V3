import csv,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_is_priority4_eur503_ld_blocks_shrinkage import SHRINKAGE,EXPECTED
class TestEURBlockShrinkage(unittest.TestCase):
 def test_grid_and_numeric_bounds_when_present(self):
  p=Path("/srv/is-analysis/results/is/audits/is_priority4_eur503_ld_blocks_shrinkage_20261011_v1/IS_PRIORITY4_EUR503_250KB_BLOCK_SHRINKAGE_QC.tsv")
  if not p.exists():self.skipTest("Reference LD grid not mounted")
  with p.open() as f: r=list(csv.DictReader(f,delimiter="\t"))
  self.assertEqual(len(r),100)
  for x in r:
   self.assertIn(x["locus"],EXPECTED)
   a=float(x["shrinkage_alpha"])
   self.assertIn(a,SHRINKAGE)
   self.assertLessEqual(int(x["effective_rank"]),502)
   self.assertEqual(x["rank_deficient"],"True")
   if a>0:
    self.assertGreater(float(x["estimated_lambda_min"]),0)
    self.assertGreaterEqual(float(x["estimated_condition_number"]),1)
 def test_theory_full_locus_not_invertible(self):
  self.assertTrue(all(n>502 for n in EXPECTED.values()))
if __name__=="__main__":unittest.main()
