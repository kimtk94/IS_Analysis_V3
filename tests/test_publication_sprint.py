import importlib.util
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location("publication_sprint", Path(__file__).resolve().parents[1]/"scripts/publication_sprint.py")
p = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p)

class SprintTests(unittest.TestCase):
    def test_eas_direction_not_selection_and_missing_stays_family(self):
        base = {"eur_egfr_fdr":".01","support_same_risk_direction_n":"3","support_p05_same_risk_direction_n":"2","eur_egfr_rsid":"rs1","eur_egfr_beta":"1","eas_egfr_available":"1","eas_same_top_snp":"1","eas_egfr_rsid":"rs1","eas_egfr_beta":"1","eas_egfr_p":".02"}
        rows = [dict(base,protein_id="a",gene_symbol="A"),dict(base,protein_id="b",gene_symbol="B",eas_egfr_beta="-1"),dict(base,protein_id="c",gene_symbol="C",eas_egfr_available="0")]
        r = p.ckd_summary(rows)
        self.assertEqual(r["eur_selected_family"],3)
        self.assertEqual(r["eas_same_anchor_tested"],2)
        self.assertEqual(r["eas_bonferroni_direction_concordant"],[])
        self.assertAlmostEqual(r["candidates"][0]["eas_bonferroni"],.06)
    def test_ocular_missing_is_not_nonsignificant(self):
        d = {"a": (2,.01,"A"), "b": (1,.01,"B"),"c": (1,.01,"C")}
        pool = {"a": (2,.01,"A"),"b": (-1,.01,"B"),"c": (1,.01,"C")}
        ext = {"a":{"log2fc":1,"fdr":.01},"b":{"log2fc":-1,"fdr":.01},"c":{"log2fc":1,"fdr":None}}
        r = p.ocular_summary(d,pool,ext)
        self.assertEqual(r["shared_gene_ids"],3)
        self.assertEqual(r["thresholds"][0]["donor_and_pool_fdr05_same_sign"]["external_fdr_missing"],1)
        self.assertEqual(r["three_model_abs1_n"],1)
    def test_bad_probability_rejected(self):
        with self.assertRaises(ValueError): p.probability("nan")
    def test_html_is_not_vcf(self):
        self.assertFalse(p.inspect_vcf_prefix(b"<html>login</html>")["vcf_header"])

if __name__ == "__main__":
    unittest.main()
