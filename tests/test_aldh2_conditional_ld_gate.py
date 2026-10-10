import importlib.util
import sys
import unittest
from pathlib import Path

source=Path(__file__).resolve().parents[1]/"scripts/is/audit_aldh2_conditional_ld_gate.py"
spec=importlib.util.spec_from_file_location("aldh2_gate",source)
gate=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=gate
spec.loader.exec_module(gate)

def row(i,jpt_p="0.1",eas_p="0.1",jpt_r2="0.05",eas_r2="0.05"):
    return {"variant_id":f"12:{i}:G:A","rsid":f"rs{i}","jpt_p_cond":jpt_p,
            "eas_p_cond":eas_p,"jpt_r2_to_rs671":jpt_r2,"eas_r2_to_rs671":eas_r2,
            "jpt_z_cond":"2.0","eas_z_cond":"2.0"}

class TestConditionalGate(unittest.TestCase):
    def test_high_ld_residual_not_independent(self):
        rows=[row(1,"0.01","1e-10","0.9","0.98"),
              row(2,"0.15","2e-8","0.89","0.99")]
        result=gate.analyze(rows,"eas")
        self.assertEqual(result["gws_residual"],2)
        self.assertEqual(result["gws_low_ld_residual"],0)

    def test_independent_candidate_flag(self):
        rows=[row(1,"0.01","1e-10","0.9","0.1")]
        self.assertEqual(gate.analyze(rows,"eas")["gws_low_ld_residual"],1)

    def test_fail_on_nonfinite_p(self):
        with self.assertRaises(ValueError):
            gate.analyze([row(1,"nan")],"jpt")

    def test_gate_preserves_blocked(self):
        rows=[row(i+1) for i in range(3033)]
        rows[0]["eas_p_cond"]="1e-10"
        rows[0]["eas_r2_to_rs671"]="0.98"
        alcohol={"ALDH2":{"source_input_variants":574,"LD_mismatch_s":0.0908385}}
        finite={"total_models":4,"reliability_warning_all_4_models":True,
                "sensitivity_warning_all_4_models":True}
        stable={"results":{"ALDH2":{"status":"BLOCKED",
                "R_reliability_flag":[True,True],
                "R_sensitivity_flag":[True,True],"credible_set_counts":[2,2]}}}
        result=gate.science_gate(rows,alcohol,finite,stable)
        self.assertFalse(result["ALDH2_BLOCKED_lifted"])
        self.assertEqual(result["sensitivity_EAS"]["gws_low_ld_residual"],0)

if __name__=="__main__":unittest.main()
