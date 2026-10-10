"""Guarded independent Korean/Japanese rs671 exposure-direction audit."""
import csv,json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_rs671_kcps2_japan_replication import build
def save(path,data):
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter="\t")
        w.writeheader();w.writerows(data)
class KoreanJapanese(unittest.TestCase):
    def fixtures(self,root):
        save(root/"G0022_4CS_KOYANAGI2024_ALCOHOL_INTAKE_EXPOSURE.tsv",[
         {"variant_grch37":"12:112241766:G:A","harmonized_effect_allele_ALT":"A",
          "beta_ALT":-1.32,"se":.0074,"p":0,"ALT_EAF":.25,"metaGWAS_variant_n":154570}])
        save(root/"G0022_RS671_KCPS2_KOREAN_ALCOHOL_AMOUNT.tsv",[
         {"variant":"ALDH2_rs671_GRCh37_12_112241766_G_A",
          "effect_allele_harmonized":"A","KCPS2_alcohol_beta_A":-.59,
          "KCPS2_alcohol_se":.0053,"KCPS2_ALCO_AMOUNT_p":0,"KCPS2_af_A":.155,
          "KCPS2_sample_size_at_variant":131743}])
        (root/"G0022_RS671_KCPS2_KOREAN_SOURCE_AUDIT.json").write_text(json.dumps({
            "selected_ZIP_member_CRC32_verified":True,"selected_member_SHA256_verified":True}))
    def test_ancestry_direction_replication_not_size_meta(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);self.fixtures(root)
            x=build(root)
            self.assertTrue(x["alcohol_association_direction_both_negative"])
            self.assertEqual(x["A_frequency_difference_Japanese_minus_Korean"],.095)
            self.assertFalse(x["Japanese_Korean_beta_magnitude_test_valid"])
            self.assertFalse(x["causal_alcohol_to_BP_to_stroke_mediation_tested"])
            self.assertEqual(x["Korean_BP_traits_verified"],[])
    def test_includes_verified_BP_but_requires_matching_allele(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);self.fixtures(root)
            (root/"G0022_RS671_KCPS2_KOREAN_SBP_GWAS_SUMMARY.json").write_text(json.dumps({
              "selected_member_CRC32_verified":True,"selected_member_SHA256_verified":True}))
            save(root/"G0022_RS671_KCPS2_KOREAN_SBP_GWAS.tsv",[
              {"variant":"12:112241766:G:A","harmonized_effect_allele":"A",
               "KCPS2_beta_A":-.08,"KCPS2_se":.01,"KCPS2_p":1e-9,
               "KCPS2_af_A":.15,"KCPS2_per_variant_N":120000}])
            self.assertEqual(build(root)["Korean_BP_traits_verified"],["SBP"])
            save(root/"G0022_RS671_KCPS2_KOREAN_SBP_GWAS.tsv",[
              {"variant":"12:112241766:G:A","harmonized_effect_allele":"G",
               "KCPS2_beta_A":-.08,"KCPS2_se":.01,"KCPS2_p":1e-9,
               "KCPS2_af_A":.15,"KCPS2_per_variant_N":120000}])
            with self.assertRaisesRegex(ValueError,"allele mismatch"):
                build(root)
if __name__=="__main__":unittest.main()
