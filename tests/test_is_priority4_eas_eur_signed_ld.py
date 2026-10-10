import math
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_is_priority4_eas_eur_signed_ld import alt_dosage,signed_r,load_raw

class ReferenceLDSafetyTests(unittest.TestCase):
 def test_ref_count_is_flipped_to_alt_dosage(self):
  counts=np.array([0.,1.,2.,np.nan])
  alt=alt_dosage(counts,"A","A","G")
  np.testing.assert_allclose(alt[:3],[2,1,0])
  self.assertTrue(math.isnan(alt[-1]))
  np.testing.assert_allclose(alt_dosage(counts,"G","A","G")[:3],[0,1,2])
 def test_allele_mismatch_is_fail_closed(self):
  with self.assertRaises(ValueError):
   alt_dosage(np.array([1.]),"T","A","G")
 def test_signed_r_direction_and_monorphic_gate(self):
  x=np.array([0.,0.,1.,1.,2.,2.])
  self.assertAlmostEqual(signed_r(x,2-x,min_pairs=3)[0],-1.)
  self.assertAlmostEqual(signed_r(x,x,min_pairs=3)[0],1.)
  self.assertTrue(math.isnan(signed_r(x,np.ones_like(x),min_pairs=3)[0]))
 def test_real_eur_eas_panel_sample_counts_if_mounted(self):
  root=Path("/srv/is-analysis/results/is/audits/is_priority4_eur_eas_signed_ld_20261011_v1")
  if not (root/"L001.EAS504.raw").exists():self.skipTest("Research server reference genotypes unavailable")
  e=load_raw(root/"L001.EAS504.raw",504)
  u=load_raw(root/"L001.EUR503.raw",503)
  self.assertEqual(len(e),1694)
  self.assertEqual(set(e),set(u))
  v="4:80681087:G:C"
  self.assertAlmostEqual(float(np.mean(e[v])/2),0.8055555556,places=7)
if __name__=="__main__":unittest.main()
