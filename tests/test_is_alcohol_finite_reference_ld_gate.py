"""Finite-reference SuSiE RSS experimental science-gate synthetic fixtures."""
import csv,json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_is_alcohol_finite_ld_method_sandbox import build
def fixture(root):
    sandbox=root/"sandbox";study=root/"study"
    (sandbox/"packages").mkdir(parents=True)
    for x in ("susieR_0.16.6","cpp11armadillo_0.5.4"):
        (sandbox/"packages"/f"{x}.tar.gz").write_bytes(x.encode("utf8"))
    for locus,n,s in (("ADH1B",637,.0133),("ALDH2",574,.0908)):
        src=study/locus;src.mkdir(parents=True)
        (src/"input_qc.json").write_text(json.dumps({
           "status":"ALCOHOL_REGIONAL_REAL_GWAS_EAS_LD_INPUT_PASS",
           "reference_EAS_n":504}))
        (src/"ALCOHOL_SUSIE_DIAGNOSTIC_AND_GATE.json").write_text(
            json.dumps({"snps":n,"median_sample_n":154570,"LD_mismatch_s":s}))
        target=sandbox/"models"/locus
        target.mkdir(parents=True)
        for mode in ("FINITE_REF_504","FINITE_REF_504_EB_MISMATCH"):
            diagnostics={"diagnostics_provided":True,"B":504,"r_over_B":.02,
               "R_reliability_flag":True,"R_sensitivity_flag":True,
               "penalty_median":3,"penalty_max":10,
               "lambda_bias":0,"B_corrected":504}
            (target/f"{mode}.qc.json").write_text(json.dumps({
              "finished_without_error":True,
              "converged":True,"susieR_version":"0.16.6",
              "reference_actual_B":504,"R_finite_enabled":True,
              "R_mismatch":"eb" if mode.endswith("EB_MISMATCH") else "none",
              "original_snp_count":n,"source_median_n":154570,
              "source_alleles_signed_and_order_QC":True,
              "finite_LD_reference_correction_sandbox_only":True,
              "credible_sets_exploratory":2,
              "top_variant":"12:1000:G:A","max_PIP":.9,
              "original_no_correction_s":s,
              "warnings":["QC reference reliability warning"],
              "diagnostics":diagnostics}))
    return sandbox,study
class FiniteRSSGate(unittest.TestCase):
    def test_all_risk_flagged_not_causal(self):
        with tempfile.TemporaryDirectory() as temp:
            sandbox,study=fixture(Path(temp))
            result=build(sandbox,study)
            self.assertTrue(result["reliability_warning_all_4_models"])
            self.assertFalse(result["ALDH2_prior_BLOCKED_gate_lifted"])
            self.assertEqual(result["final_causal_genes"],0)
            self.assertEqual(result["total_models"],4)
            self.assertTrue(result["server_default_susieR_package_unmodified"])
    def test_missing_flag_or_wrong_version_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            sandbox,study=fixture(Path(temp))
            p=sandbox/"models"/"ALDH2"/"FINITE_REF_504_EB_MISMATCH.qc.json"
            j=json.loads(p.read_text())
            j["diagnostics"]["R_reliability_flag"]=False
            p.write_text(json.dumps(j))
            with self.assertRaisesRegex(ValueError,"reliability"):
                build(sandbox,study)
            j["diagnostics"]["R_reliability_flag"]=True
            j["susieR_version"]="0.14.2"
            p.write_text(json.dumps(j))
            with self.assertRaisesRegex(ValueError,"version"):
                build(sandbox,study)
    def test_source_variant_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            sandbox,study=fixture(Path(temp))
            p=study/"ADH1B"/"ALCOHOL_SUSIE_DIAGNOSTIC_AND_GATE.json"
            j=json.loads(p.read_text());j["snps"]=636
            p.write_text(json.dumps(j))
            with self.assertRaisesRegex(ValueError,"Variant set"):
                build(sandbox,study)
if __name__=="__main__":unittest.main()
