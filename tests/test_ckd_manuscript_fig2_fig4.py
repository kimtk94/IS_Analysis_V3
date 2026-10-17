import csv,hashlib,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from plot_ckd_manuscript_fig2_fig4 import build,GENES
class Figure24(unittest.TestCase):
 def test_synthetic(self):
  with tempfile.TemporaryDirectory() as d:
   root=Path(d);p=root/"src.tsv"
   fields=["gene_symbol","eas_egfr_wald_beta","eas_egfr_wald_se","eas_egfr_wald_p","top_cell_type","kidney_eqtl_meta686_hits","kidney_eqtl_tubule356_hits","kidney_eqtl_glomerulus303_hits","kidney_pqtl_maintext"]
   with p.open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter="\t");w.writeheader()
    for g in GENES:w.writerow(dict(gene_symbol=g,eas_egfr_wald_beta=".03",eas_egfr_wald_se=".01",eas_egfr_wald_p=".02",top_cell_type="test_cell",kidney_eqtl_meta686_hits="4",kidney_eqtl_tubule356_hits="1",kidney_eqtl_glomerulus303_hits="NA",kidney_pqtl_maintext="1"))
   sha=hashlib.sha256(p.read_bytes()).hexdigest()
   with self.assertRaisesRegex(ValueError,"SHA mismatch"):build(p,root/"bad","0"*64)
   result=build(p,root/"out",sha)
   self.assertEqual(result["mr_estimates"],9)
   self.assertTrue((root/"out"/"CKD_FIG2_EAS_ANCHOR_WALD.svg").is_file())
   self.assertTrue((root/"out"/"CKD_FIG4_KIDNEY_TISSUE_QTL.svg").is_file())
   with self.assertRaisesRegex(ValueError,"overwrite"):build(p,root/"out",sha)
if __name__=="__main__":unittest.main()
