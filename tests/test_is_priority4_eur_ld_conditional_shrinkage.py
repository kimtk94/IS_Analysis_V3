import csv,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_is_priority4_eur_ld_conditional_shrinkage import ALPHAS,LEADS
class TestEURLDConditionalDiagnostic(unittest.TestCase):
 def test_reference_gate_and_alpha_counts(self):
  root=Path("/srv/is-analysis/results/is/audits/is_priority4_eur_ref_lead_conditional_20261011_v2")
  fp=root/"IS_PRIORITY4_EUR_REF_LEAD_CONDITIONAL_ALPHA_SUMMARY.tsv"
  if not fp.exists():self.skipTest("Research server reference-only output not present")
  with fp.open() as f:rows=list(csv.DictReader(f,delimiter="\t"))
  self.assertEqual(len(rows),10)
  for x in rows:
   self.assertIn(x["locus"],LEADS)
   self.assertIn(float(x["alpha"]),ALPHAS)
   self.assertEqual(x["no_independence_inference"],"True")
   self.assertGreater(int(x["n_evaluable"]),2800)
  chr4=[r for r in rows if r["locus"]=="BBJ_IS_L001"]
  self.assertTrue(all(int(x["n_absZ_gt_5p45"])==0 for x in chr4))
  chr10=[r for r in rows if r["locus"]=="BBJ_IS_L002"]
  self.assertEqual([int(x["n_absZ_gt_5p45"]) for x in chr10],[1,1,1,1,0])
if __name__=="__main__":unittest.main()
