"""BBJ PheWeb rs671 ALT blood-pressure source and study-overlap barriers."""
import json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_g0022_rs671_bbj_bp import source_variant,extract,PHENOS
def html(pheno):
    v={"chrom":"12","pos":112241766,"ref":"G","alt":"A",
       "rsids":["rs671"],"phenos":pheno}
    return "<html><script>window.variant = "+json.dumps(v)+";</script></html>"
def measure(p):
    return {"phenocode":p,"phenostring":p,"beta":-.06,
     "sebeta":.004,"pval":1.2e-50,"af":.25,"num_samples":145000,
     "citation":"SakaueKanai2021","category":"Blood pressure"}
class BPTests(unittest.TestCase):
    def test_complete_original_bp_effect_and_no_mediation_promotion(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)
            src=root/"html";src.write_text(html([measure(p) for p in PHENOS]))
            r=extract(src,root/"out")
            self.assertEqual(r["BP_traits"],4)
            self.assertTrue(r["all_BP_beta_ALT_A_negative"])
            self.assertFalse(r["mmHg_trait_units_assumed"])
            self.assertFalse(r["mediation_identified"])
            self.assertEqual(r["stroke_cohort_independence"],"UNVERIFIED")
    def test_wrong_variant_allele_and_missing_source_trait_rejected(self):
        records=[measure(p) for p in PHENOS]
        s=html(records)
        self.assertEqual(source_variant(s)["alt"],"A")
        with self.assertRaisesRegex(ValueError,"Wrong variant"):
            source_variant(s.replace('"alt": "A"','"alt": "T"'))
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);src=root/"html"
            src.write_text(html(records[:-1]))
            with self.assertRaisesRegex(ValueError,"missing BBJ BP"):
                extract(src,root/"out")
    def test_wrong_study_source_rejected(self):
        r=[measure(p) for p in PHENOS];r[0]["citation"]="Unknown"
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);src=root/"html"
            src.write_text(html(r))
            with self.assertRaisesRegex(ValueError,"Wrong trait"):
                extract(src,root/"out")
if __name__=="__main__":unittest.main()
