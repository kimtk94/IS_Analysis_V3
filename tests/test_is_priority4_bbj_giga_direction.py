import csv,unittest
from pathlib import Path
ROOT=Path("/srv/is-analysis/results/is/audits/is_priority4_bbj_giga_direction_20261011_v1")
class TestBBJGigaDirection(unittest.TestCase):
 def test_fixed_source_counts(self):
  p=ROOT/"IS_PRIORITY4_BBJ_GIGA_ALLELE_EXACT_DIRECTION.tsv"
  if not p.exists():self.skipTest("Frozen science output unavailable")
  with p.open() as f:r={x["gene"]:x for x in csv.DictReader(f,delimiter="\t")}
  self.assertEqual(set(r),{"FGF5","C4orf22","CALHM2","NEURL"})
  self.assertEqual([int(r[x]["both_stroke_GWAS_p_lt_5e8"]) for x in ("FGF5","C4orf22","CALHM2","NEURL")],[1,1,0,0])
  self.assertEqual([int(r[x]["same_BBJ_GIGA_ALT_beta_direction"]) for x in ("FGF5","C4orf22","CALHM2","NEURL")],[2216,2215,2292,2281])
  self.assertTrue(all(x["SNPs_not_independent_due_to_LD"]=="True" for x in r.values()))
  self.assertTrue(all(x["cohorts_not_attested_independent"]=="True" for x in r.values()))
