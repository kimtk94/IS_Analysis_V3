import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_is_eur_b37_liftover_coverage import parse_chain,lift
class ChainRegression(unittest.TestCase):
 def test_original_b37_to_b38_known_markers(self):
  p=Path("/srv/is-analysis/data/is/reference/hg19ToHg38.over.chain.gz")
  if not p.exists():self.skipTest("Reference chain unavailable")
  mapper=parse_chain(p)
  for chrom,src,dst in [("12",112241766,111803962),("12",112168009,111730205),("10",104839152,103079395)]:
   with self.subTest(chrom=chrom,variant=src):self.assertEqual(lift(mapper,chrom,src),dst)
if __name__=="__main__":unittest.main()
