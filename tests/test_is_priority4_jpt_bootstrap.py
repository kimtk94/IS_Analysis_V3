import csv,unittest
from pathlib import Path
ROOT=Path("/srv/is-analysis/results/is/audits/is_priority4_jpt_bootstrap_20261011_v1")
def data(name):
 with (ROOT/name).open() as f:return list(csv.DictReader(f,delimiter="\t"))
class TestJPTBoot(unittest.TestCase):
 def test_reference_counts(self):
  if not ROOT.exists():self.skipTest("No reference output")
  a=data("JPT_REFERENCE_LOCUS_QC.tsv")
  self.assertEqual(len(a),8)
  self.assertEqual({int(x["n_samples"]) for x in a},{104,400,504,503})
 def test_bootstrap_gate(self):
  if not ROOT.exists():self.skipTest("No reference output")
  a={x["population"]:x for x in data("CHR10_JPT_EAS_EUR_BOOTSTRAP.tsv")}
  self.assertEqual(len(a),6)
  self.assertLess(float(a["JPT104"]["bootstrap_r_CI_high"]),0)
  self.assertGreater(float(a["EUR503"]["bootstrap_r_CI_low"]),0)
  self.assertLess(float(a["JPT104_minus_EUR503"]["bootstrap_r_CI_high"]),0)
  self.assertGreater(float(a["JPT104_minus_otherEAS400"]["bootstrap_r_CI_high"]),0)
 def test_jpt_alpha_grid_is_exploratory(self):
  if not ROOT.exists():self.skipTest("No reference output")
  a=data("JPT_EAS_EUR_ALPHA_SCREEN.tsv")
  self.assertEqual(len(a),40)
  self.assertTrue(all(x["inference_gate"]=="REFERENCE_ONLY_NOT_CONDITIONAL_GWAS" for x in a))
  self.assertTrue(all(int(x["heuristic_abs_z_ge_5p45"])==0 for x in a if x["population"]=="JPT104"))
 def test_candidate_alpha_reference(self):
  if not ROOT.exists():self.skipTest("No reference output")
  a={(x["population"],float(x["alpha"])):x for x in data("CHR10_JPT_EAS_EUR_CANDIDATE_ALPHA.tsv")}
  self.assertEqual(len(a),20)
  self.assertAlmostEqual(float(a["JPT104",.05]["signed_r"]),-.24387507877062556)
  self.assertAlmostEqual(float(a["EUR503",.05]["signed_r"]),.1092056318737335)
  self.assertFalse(a["JPT104",.05]["abs_z_ge_5p45"]=="True")
