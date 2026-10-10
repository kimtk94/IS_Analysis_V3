"""Conditional signed-LD GWAS diagnostic is not formal conditional association."""
import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_g0022_rs671_conditional_ld import residualized_z,p_from_z
class ConditionalTest(unittest.TestCase):
  def test_positive_signed_ld(self):
    r=residualized_z([4.,2.],[1.,.5],0)
    self.assertTrue(np.isnan(r[0]))
    self.assertAlmostEqual(r[1],0.)
  def test_negative_signed_ld(self):
    r=residualized_z([4.,-2.],[1.,-.5],0)
    self.assertAlmostEqual(r[1],0.)
  def test_no_correlation_preserves_z(self):
    r=residualized_z([5.,3.],[1.,0.],0)
    self.assertAlmostEqual(r[1],3.)
    self.assertAlmostEqual(p_from_z(0.),1.)
  def test_fail_closed_ld_mismatch(self):
    with self.assertRaises(ValueError):
       residualized_z([3.,2.],[.8,.5],0)
    with self.assertRaises(ValueError):
       residualized_z([3.,2.],[1.,1.2],0)
    with self.assertRaises(ValueError):
       residualized_z([3.,2.],[1.],0)
if __name__=="__main__":unittest.main()
