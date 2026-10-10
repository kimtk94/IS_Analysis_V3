"""Korean vs Japanese rs671 blood pressure association direction audit tests."""
import csv,json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_is_rs671_kcps2_bbj_bp_convergence import run
def write(path,rows):
    with path.open("w") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        w.writeheader();w.writerows(rows)
class BPCompare(unittest.TestCase):
    def prepare(self,root):
        write(root/"G0022_RS671_BBJ_BLOOD_PRESSURE_EFFECTS.tsv",
        [dict(phenotype=p,variant_grch37="12:112241766:G:A",effect_allele="ALT_A",
           beta_in_reported_trait_scale="-.063",se=".004",p="1e-30",n_total="145500")
         for p in ["SBP","DBP"]])
        for p in ["SBP","DBP"]:
            (root/f"G0022_RS671_KCPS2_KOREAN_{p}_GWAS_SUMMARY.json").write_text(json.dumps({
              "selected_member_CRC32_verified":True,"selected_member_SHA256_verified":True}))
            write(root/f"G0022_RS671_KCPS2_KOREAN_{p}_GWAS.tsv",[{
              "variant":"12:112241766:G:A","harmonized_effect_allele":"A",
              "KCPS2_beta_A":"-.095","KCPS2_se":".005",
              "KCPS2_p":"1e-10","KCPS2_per_variant_N":"150000"}])
    def test_both_concordant_and_not_claimed_causal(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);self.prepare(root);report=run(root)
            self.assertEqual(report["verified_BP_trait_pairs"],["SBP","DBP"])
            self.assertTrue(report["korean_japanese_BP_beta_direction_concordant"])
            self.assertFalse(report["AIS_mediated_alcohol_effect_identified"])
            self.assertFalse(report["korean_BBJ_BP_beta_units_identical"])
    def test_wrong_allele_fails(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);self.prepare(root)
            write(root/"G0022_RS671_KCPS2_KOREAN_SBP_GWAS.tsv",[{
              "variant":"12:112241766:G:A","harmonized_effect_allele":"G",
              "KCPS2_beta_A":"-.095","KCPS2_se":".005",
              "KCPS2_p":"1e-10","KCPS2_per_variant_N":"150000"}])
            with self.assertRaisesRegex(ValueError,"allele"):
                run(root)
if __name__=="__main__":unittest.main()
