import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_is_30_abf_joint_sensitivity import run
class TestJointSensitivity(unittest.TestCase):
 def test_bounded_original_joint_counts(self):
  base=Path("/srv/is-analysis/results/is/audits")
  paths=[base/"is_30_abf_raw_input_risk_20261011_v1/IS_30_ABF_RAW_SNP_INPUT_RISKS.tsv",base/"is_646_qtl_scalar_n_sensitivity_20261011_v1/IS_646_COLOC_QTL_N_SCALAR_SENSITIVITY.tsv",base/"legacy_coloc_p12_sensitivity_20261010_v1/IS_LEGACY_646_ABF_P12_GRID.tsv"]
  if not all(p.exists() for p in paths): self.skipTest("source files unavailable")
  with tempfile.TemporaryDirectory() as t:
   x=run(*paths,Path(t))
   self.assertEqual(x["assays"],30)
   self.assertEqual(x["unique_gene_ids"],14)
   self.assertEqual(x["p12_cross_H4_0p5"],18)
   self.assertEqual(x["p12_cross_H4_0p8"],12)
   self.assertEqual(x["scalar_N_cross_H4_0p5"],1)
   self.assertEqual(x["missing_inputs_unchanged"],616)
if __name__=="__main__":unittest.main()
