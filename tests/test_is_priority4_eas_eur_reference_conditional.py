import csv,math,sys,tempfile,unittest,xml.etree.ElementTree as ET
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_is_priority4_eas_eur_reference_conditional import ALPHA,LOCI,CANDIDATE
from render_is_priority4_ancestry_candidate_sensitivity import render
ROOT=Path("/srv/is-analysis/results/is/audits/is_priority4_eas_eur_matched_reference_20261011_v1")
class TestIndependentReferenceAncestryLD(unittest.TestCase):
 def load(self,name):
  f=ROOT/name
  if not f.exists():self.skipTest("Frozen original research audit not mounted")
  with f.open(newline="") as fp:return list(csv.DictReader(fp,delimiter="\t"))
 def test_complete_same_snp_allele_reference(self):
  loc={r["locus"]:r for r in self.load("IS_PRIORITY4_EAS_EUR_LEAD_LD_SUMMARY.tsv")}
  self.assertEqual(set(loc),set(LOCI))
  self.assertEqual({k:int(v["n_exact_shared_SNPs"]) for k,v in loc.items()},
                   {"BBJ_IS_L001":2834,"BBJ_IS_L002":2828})
  self.assertEqual({k:int(v["both_ancestry_lead_r_informative"]) for k,v in loc.items()},
                   {"BBJ_IS_L001":2815,"BBJ_IS_L002":2765})
  self.assertEqual({k:int(v["direction_opposite_abs_r_ge_0p1_both"]) for k,v in loc.items()},
                   {"BBJ_IS_L001":50,"BBJ_IS_L002":91})
  self.assertTrue(all(r["formal_GWAS_conditional_or_SuSiE_performed"]=="False" for r in loc.values()))
 def test_candidate_sign_flip_and_approx_score(self):
  rows=self.load("IS_CHR10_105694487_EAS_EUR_CANDIDATE_SENSITIVITY.tsv")
  self.assertEqual(len(rows),10)
  a={(r["ancestry"],float(r["alpha"])):r for r in rows}
  self.assertEqual(len(a),10)
  self.assertAlmostEqual(float(a["EAS504",0.0]["lead_candidate_signed_r"]),-0.14790281784322776,places=10)
  self.assertAlmostEqual(float(a["EUR503",0.0]["lead_candidate_signed_r"]),0.1092056318737335,places=10)
  self.assertAlmostEqual(float(a["EAS504",0.05]["approx_lead_residual_z"]),4.011050497669197,places=9)
  self.assertAlmostEqual(float(a["EUR503",0.05]["approx_lead_residual_z"]),5.4983144023559225,places=9)
  self.assertTrue(all(a["EAS504",x]["absZ_gt_5p45"]=="False" for x in ALPHA))
  self.assertEqual([a["EUR503",x]["absZ_gt_5p45"] for x in ALPHA],
                   ["True","True","True","True","False"])
 def test_common_reference_screen_counts(self):
  rr=self.load("IS_PRIORITY4_EAS_EUR_APPROX_CONDITIONAL_SCREEN.tsv")
  self.assertEqual(len(rr),10)
  self.assertTrue(all(r["not_formal_conditional_GWAS"]=="True" for r in rr))
  self.assertTrue(all(int(r["EAS504_approx_residual_absZ_gt_5p45"])==0 for r in rr))
  self.assertTrue(all(int(r["EUR503_approx_residual_absZ_gt_5p45"])==0 for r in rr if r["locus"]=="BBJ_IS_L001"))
  self.assertEqual([int(r["EUR503_approx_residual_absZ_gt_5p45"]) for r in rr if r["locus"]=="BBJ_IS_L002"],[1,1,1,1,0])
 def test_visualization_has_labelled_noncausal_boundary(self):
  f=ROOT/"IS_CHR10_105694487_EAS_EUR_CANDIDATE_SENSITIVITY.tsv"
  if not f.exists():self.skipTest("Missing independent audit")
  with tempfile.TemporaryDirectory() as t:
   out=Path(t)/"fig.svg"
   self.assertEqual(render(f,out),10)
   text=" ".join(x.text or "" for x in ET.parse(out).getroot().iter())
   self.assertIn("NOT a formally calibrated conditional test",text)
   self.assertIn("10:105694487:G:T",text)
if __name__=="__main__":unittest.main()
