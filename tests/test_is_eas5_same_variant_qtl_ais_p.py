"""Same physical SNP GWAS/QTL P-threshold QC; reject gene-level min-P misinterpretation."""
import importlib.util,sys,unittest
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"scripts/is/audit_is_eas5_same_variant_qtl_ais_p.py"
sp=importlib.util.spec_from_file_location("same_SNP_gate",p)
m=importlib.util.module_from_spec(sp);sys.modules[sp.name]=m;sp.loader.exec_module(m)
def rec(variant,pqtl,pais,afdiff=.02,bb=True,ais=True,pal=False):
 return {"region_id":"IS_XDATA_G0015","QTL_dataset":"QTD000131",
  "gene_id":"ENSG00001","mapped_variant_GRCh37":variant,
  "BBJ_exact_SNP":str(bb),"EAS_AIS_exact_SNP":str(ais),
  "palindromic_SNP":str(pal),"QTL_p":str(pqtl),
  "EAS_AIS_p":str(pais),"abs_QTL_vs_EAS_AIS_alt_eaf_delta":str(afdiff)}
class TestExactVariantGate(unittest.TestCase):
 def test_two_different_low_p_variants_not_joint(self):
  rows=[rec("1:10:A:G",1e-8,.001),rec("1:20:A:G",.001,1e-9)]
  genes,_,_=m.scan(rows)
  self.assertEqual(genes[("IS_XDATA_G0015","QTD000131","ENSG00001")]["n_same_SNP_QTL_1e4_AIS_5e8"],0)
 def test_joint_same_SNP_with_ancestry_af_discordance_is_flagged(self):
  genes,_,_=m.scan([rec("1:10:A:G",1e-8,1e-9,afdiff=.23)])
  s=genes[("IS_XDATA_G0015","QTD000131","ENSG00001")]
  self.assertEqual(s["n_same_SNP_QTL_1e4_AIS_5e8"],1)
  self.assertEqual(s["n_joint_AF_delta_le_0p1"],0)
 def test_palindrome_never_passes_joint_gate(self):
  genes,_,_=m.scan([rec("1:10:A:T",1e-8,1e-9,pal=True)])
  s=genes[("IS_XDATA_G0015","QTD000131","ENSG00001")]
  self.assertEqual(s["n_SNP_both_GWAS_nonpal"],0)
 def test_exact_variant_duplicates_fail(self):
  r=rec("1:10:A:G",1e-9,1e-10)
  with self.assertRaises(ValueError):m.scan([r,r.copy()])
if __name__=="__main__":unittest.main()
