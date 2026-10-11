import csv, unittest
from pathlib import Path
P=Path("/srv/is-analysis/results/is/audits/is_chr10_two_lead_ld_tags_20261011_v1/IS_CHR10_EAS_JPT_TWO_LEAD_LD_QC.tsv")
class TestChr10LDTags(unittest.TestCase):
 def test_descriptive_not_credible_set(self):
  if not P.exists():self.skipTest("Audit output absent")
  with P.open() as f:r={(x["population"],x["lead"]):x for x in csv.DictReader(f,delimiter="\t")}
  self.assertEqual(len(r),6)
  self.assertEqual(int(r["EAS504","BBJ"]["n_r2_ge_0p8"]),6)
  self.assertEqual(int(r["EAS504","GIGASTROKE"]["n_r2_ge_0p8"]),11)
  self.assertEqual(int(r["JPT104","BBJ"]["n_QTD_CALHM2_nominal_p_lt_1e4_and_r2_ge_0p8"]),3)
  self.assertTrue(all(x["not_independent_signal_or_credible_set"]=="True" for x in r.values()))
