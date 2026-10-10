"""Regression for EAS PVAR AF orientation and strict reference identity."""
import importlib.util,sys,tempfile,unittest
from pathlib import Path
s=Path(__file__).resolve().parents[1]/"scripts/is/audit_is_8956_maf_palindrome_reference.py"
spec=importlib.util.spec_from_file_location("is_8956_af",s)
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m
spec.loader.exec_module(m)
class TestPVARRead(unittest.TestCase):
 def test_eas_af_is_alt_not_minor_af(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/"x.pvar"
   p.write_text("#CHROM\tPOS\tID\tREF\tALT\tINFO\n"+
     "12\t112090022\t12:112090022:C:A\tC\tA\tEAS_AF=0.999;EUR_AF=0.32\n")
   q=m.pvar(p,{"12:112090022:C:A"})
   self.assertAlmostEqual(q["12:112090022:C:A"],.999)
   self.assertAlmostEqual(min(q["12:112090022:C:A"],1-q["12:112090022:C:A"]),.001)
 def test_pvar_variant_mismatch_rejected(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/"x.pvar"
   p.write_text("#CHROM\tPOS\tID\tREF\tALT\tINFO\n"+
     "12\t112090022\t12:112090022:C:A\tA\tC\tEAS_AF=0.999\n")
   with self.assertRaises(ValueError):m.pvar(p,{"12:112090022:C:A"})
 def test_duplicate_pvar_variant_rejected(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/"x.pvar"
   line="12\t112090022\t12:112090022:C:A\tC\tA\tEAS_AF=0.999\n"
   p.write_text("#CHROM\tPOS\tID\tREF\tALT\tINFO\n"+line+line)
   with self.assertRaises(ValueError):m.pvar(p,{"12:112090022:C:A"})
if __name__=="__main__":unittest.main()
