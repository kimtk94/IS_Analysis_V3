"""rs671 alcohol/BP/AIS evidence aggregation must never promote mediation."""
import csv,json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_g0022_rs671_allele_triangulation import process
def write(path,lines):
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(lines[0]),delimiter="\t")
        w.writeheader();w.writerows(lines)
class Triangulation(unittest.TestCase):
    def fixture(self,base):
        for trait,coef in (("ALCOHOL_INTAKE","-.4"),("DRINKING_STATUS","-.6")):
            (base/f"G0022_KOYANAGI2024_{trait}_AUDIT.json").write_text(json.dumps({
              "official_MD5_verified":True,"rs671_present":True,
              "drinking_status_event_ever_drinker_verified_from_published_methods":trait=="DRINKING_STATUS"}))
            write(base/f"G0022_4CS_KOYANAGI2024_{trait}_EXPOSURE.tsv",
                [{"variant_grch37":"12:112241766:G:A","harmonized_effect_allele_ALT":"A",
                "beta_ALT":coef,"se":".02","p":"2e-10","metaGWAS_variant_n":145500}])
        write(base/"G0022_RS671_BBJ_BLOOD_PRESSURE_EFFECTS.tsv",
            [{"phenotype":"SBP","variant_grch37":"12:112241766:G:A",
              "beta_in_reported_trait_scale":"-.062","se":".0039",
              "p":"1.2e-56","n_total":145505}])
        write(base/"G0022_RS671_SIX_GWAS_SUBTYPE_EFFECTS.tsv",
             [{"trait":"EAS_AIS","variant":"12:112241766:G:A",
               "ALT_beta":"-.1514","se":".0176","p":"8.426e-18"}])
    def test_all_four_associations_align_without_causal_claim(self):
        with tempfile.TemporaryDirectory() as tmp:
            base=Path(tmp);self.fixture(base)
            x=process(base)
            self.assertTrue(x["alcohol_intake_A_decreases_log2_grams_per_day"])
            self.assertTrue(x["SBP_A_decreases_reported_BP_scale"])
            self.assertTrue(x["EAS_AIS_A_decreases_log_odds"])
            self.assertFalse(x["causal_mediation_percent_calculated"])
            self.assertFalse(x["rs671_A_as_alcohol_only_valid_IV"])
            self.assertTrue(x["drinking_status_case_direction_verified"])
            self.assertAlmostEqual(x["ever_drinker_odds_ratio_per_rs671_A"],.5488116360940264)
    def test_missing_or_nonverified_input_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            base=Path(tmp);self.fixture(base)
            (base/"G0022_KOYANAGI2024_ALCOHOL_INTAKE_AUDIT.json").write_text(json.dumps({
              "official_MD5_verified":False,"rs671_present":True}))
            with self.assertRaisesRegex(ValueError,"MD5"):
                process(base)
            (base/"G0022_4CS_KOYANAGI2024_ALCOHOL_INTAKE_EXPOSURE.tsv").unlink()
            with self.assertRaises(FileNotFoundError):
                process(base)
if __name__=="__main__":unittest.main()
