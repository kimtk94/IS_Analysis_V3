import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from ocular_frozen_241_external import main, summarize

class SummaryTest(unittest.TestCase):
    def test_local_replay_checks_digest_without_downloading(self):
        wrong_source = Path(__file__).parent / "fixtures" / "gigastroke.tsv"
        with tempfile.TemporaryDirectory() as output:
            argv = ["replay", "--panel", "unused.zip", "--cpm", str(wrong_source), "--out", output]
            with patch.object(sys, "argv", argv), patch("ocular_frozen_241_external.read_panel", return_value={}), patch("ocular_frozen_241_external.get_source") as download:
                with self.assertRaisesRegex(ValueError, "7-day CPM hash mismatch"):
                    main()
                download.assert_not_called()
            self.assertEqual(list(Path(output).iterdir()), [])

    def test_missing_and_strata(self):
        p={"a":{"ecm_member":True},"b":{"ecm_member":False},"c":{"ecm_member":True}}
        x={"a":{"min_arm_mean_cpm":2,"mean_7d_log2_cpm_ratio":1,"acute_log2fc":1,"n_strains_same_sign":3,"n_strains_same_sign_abs_ge_025":2,"all_leave_one_strain_out_same_sign":True},"b":{"min_arm_mean_cpm":.5,"mean_7d_log2_cpm_ratio":-1,"acute_log2fc":1,"n_strains_same_sign":1,"n_strains_same_sign_abs_ge_025":1,"all_leave_one_strain_out_same_sign":False}}
        r=summarize(p,x)
        self.assertEqual(r["all_241"]["unavailable"],1)
        self.assertEqual(r["ECM_16"]["frozen_total"],2)
        self.assertEqual(r["all_241"]["strata"][0]["all_three_strains_same_sign"],1)
        self.assertEqual(r["all_241"]["strata"][0]["all_three_strains_same_sign_abs_ge_025"],0)
        self.assertEqual(r["all_241"]["strata"][1]["evaluable"],1)
        self.assertEqual(r["all_241"]["strata"][2]["evaluable"],0)
if __name__=="__main__":unittest.main()
