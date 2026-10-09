"""Expanded regional reference must supersede old BBJ narrow reference safely."""
import tempfile
import unittest
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_gws_gwas_reference_overlap import pvar_path
from audit_gws_lead_eaf import source_prefix

class ExpandedReferenceSelectionTest(unittest.TestCase):
 def test_prefer_complete_expanded_and_block_incomplete(self):
  with tempfile.TemporaryDirectory() as temp:
   root=Path(temp);old=root/"old";new=root/"new";old.mkdir();new.mkdir()
   original=old/"BBJ_IS_L002.1KG_EAS.GRCh37.pvar"
   original.write_text("#CHROM\tPOS\tID\tREF\tALT\n")
   self.assertEqual(pvar_path("IS_XDATA_G0015",new,old),original)
   expanded=new/"IS_XDATA_G0015.stable"
   Path(str(expanded)+".pvar").write_text("#CHROM\tPOS\tID\tREF\tALT\n")
   with self.assertRaisesRegex(RuntimeError,"Incomplete"):
    pvar_path("IS_XDATA_G0015",new,old)
   for ext in (".pgen",".psam"):Path(str(expanded)+ext).write_bytes(b"x")
   self.assertEqual(pvar_path("IS_XDATA_G0015",new,old),Path(str(expanded)+".pvar"))
   self.assertEqual(source_prefix("IS_XDATA_G0015",new,old),expanded)
 def test_other_reference_stays_on_original(self):
  with tempfile.TemporaryDirectory() as td:
   o=Path(td)/"old";o.mkdir();n=Path(td)/"new";n.mkdir()
   legacy=o/"BBJ_IS_L003.1KG_EAS.GRCh37.pvar"
   legacy.write_text("#CHROM\tPOS\tID\tREF\tALT\n")
   self.assertEqual(pvar_path("IS_XDATA_G0022",n,o),legacy)
if __name__=="__main__":unittest.main()
