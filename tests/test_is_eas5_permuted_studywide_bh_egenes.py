"""Validate study-wide (not IS-selected-only) Benjamini-Hochberg eQTL gates."""
import csv,gzip,importlib.util,sys,tempfile,unittest
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"scripts/is/audit_is_eas5_permuted_studywide_bh_egenes.py"
spec=importlib.util.spec_from_file_location("eas5_studywide_BH",p)
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m
spec.loader.exec_module(m)
class TestQTLStudywideBH(unittest.TestCase):
 def test_bh_exact_small_example(self):
  self.assertEqual([round(x,6) for x in m.bh_all([0.01,0.04,0.8])],[0.03,0.06,0.8])
 def test_equal_p_monotonic_not_rank_artifact(self):
  self.assertEqual([round(x,4) for x in m.bh_all([0.01,0.01,0.04,0.2])],
                   [0.02,0.02,0.0533,0.2])
 def test_invalid_p_rejected(self):
  for pp in ([-1.], [1.1], [float("nan")]):
   with self.assertRaises(ValueError):m.bh_all(pp)
 def test_permutation_genomewide_q_not_candidate_filtered(self):
  with tempfile.TemporaryDirectory() as t:
   p=Path(t)/"f.permuted.tsv.gz"
   with gzip.open(p,"wt") as f:
    f.write("molecular_trait_object_id\tn_variants\tp_beta\tvariant\n")
    f.write("ENSG001\t400\t0.01\tchr10_10000_A_G\n")
    f.write("ENSG002\t400\t0.04\tchr10_12000_A_G\n")
    f.write("ENSG003\t400\t0.8\tchr10_13000_A_G\n")
   row,n,valid=m.read_perm(p)
   self.assertEqual((n,valid),(3,3))
   self.assertAlmostEqual(row["ENSG001"]["BH_q"],.03)
   self.assertAlmostEqual(row["ENSG002"]["BH_q"],.06)
if __name__=="__main__":unittest.main()
