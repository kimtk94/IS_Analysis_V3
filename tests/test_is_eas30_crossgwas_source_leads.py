"""Test same physical variant after source-specific REF/ALT pair reorientation."""
import importlib.util,sys,unittest
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"scripts/is/audit_is_eas30_crossgwas_source_leads.py"
s=importlib.util.spec_from_file_location("lead_eas30",p)
m=importlib.util.module_from_spec(s);sys.modules[s.name]=m;s.loader.exec_module(m)
class TestLeadAlleles(unittest.TestCase):
 def test_swapped_variant_id_reorients_beta_and_eaf(self):
  r={"group_source_lead_present":"True",
    "group_source_lead_variant":"4:81181072:A:T",
    "group_source_lead_beta_ALT":"0.2","group_source_lead_ALT_AF":"0.73",
    "group_source_lead_p":"1e-8"}
  out=m.normalize("4:81181072:T:A",r)
  self.assertEqual(out["status"],"SWAPPED_REF_ALT_PAIR_REVIEW")
  self.assertAlmostEqual(out["beta"],-.2)
  self.assertAlmostEqual(out["ALT_AF"],.27)
 def test_exact_variant_no_reorientation(self):
  r={"group_source_lead_present":"True",
    "group_source_lead_variant":"4:81181072:T:A",
    "group_source_lead_beta_ALT":"-0.2","group_source_lead_ALT_AF":"0.27",
    "group_source_lead_p":"1e-8"}
  out=m.normalize("4:81181072:T:A",r)
  self.assertEqual(out["status"],"EXACT_REF_ALT_PAIR")
  self.assertAlmostEqual(out["beta"],-.2)
 def test_mismatch_variant_fails(self):
  r={"group_source_lead_present":"True","group_source_lead_variant":"4:81181072:T:C"}
  with self.assertRaises(ValueError):m.normalize("4:81181072:T:A",r)
 def test_absent_lead_is_not_negative_association(self):
  r={"group_source_lead_present":"False"}
  self.assertEqual(m.normalize("4:81181072:T:A",r)["status"],"SOURCE_LEAD_NOT_REPORTED")
if __name__=="__main__":unittest.main()
