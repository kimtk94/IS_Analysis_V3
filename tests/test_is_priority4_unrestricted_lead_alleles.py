import csv,unittest
from pathlib import Path
OUT=Path("/srv/is-analysis/results/is/audits/is_priority4_unrestricted_lead_alleles_20261011_v1/IS_PRIORITY4_UNRESTRICTED_RECIPROCAL_LEAD_SITE_ALLELES.tsv")
class TestReciprocalAlleles(unittest.TestCase):
 def test_full_source_four_leads(self):
  if not OUT.is_file():self.skipTest("Results not mounted")
  with OUT.open() as f:r=list(csv.DictReader(f,delimiter="\t"))
  self.assertEqual(len(r),4)
  self.assertEqual([x["pos"] for x in r],["81182554","81181072","105459834","105599770"])
  self.assertTrue(all(float(x["frequency_abs_difference"])<.004 for x in r))
  self.assertTrue(all(x["same_FastaREF_effect_direction"]=="True" for x in r))
  self.assertEqual(sum(x["palindrome"]=="True" for x in r),1)
  self.assertEqual(sum(x["strand_provenance_verified"]=="False" for x in r),1)
