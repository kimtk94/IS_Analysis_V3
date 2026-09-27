import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/run_ckd_stage4_eas_anchor_mr.py"

spec = importlib.util.spec_from_file_location(
    "ckd_eas_anchor_mr",
    SCRIPT,
)

mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class TestCKDEASAnchorMR(unittest.TestCase):

    def test_direct_harmonization(self):
        sign, status = mod.harmonize(
            "T", "C", "T", "C"
        )
        self.assertEqual(sign, 1.0)
        self.assertEqual(status, "direct")

    def test_reverse_harmonization(self):
        sign, status = mod.harmonize(
            "G", "A", "A", "G"
        )
        self.assertEqual(sign, -1.0)
        self.assertEqual(status, "reverse")

    def test_strand_harmonization(self):
        sign, status = mod.harmonize(
            "G", "A", "C", "T"
        )
        self.assertEqual(sign, 1.0)
        self.assertEqual(status, "strand_direct")

    def test_wald_ratio(self):
        x = mod.wald_ratio(
            beta_x=0.5,
            se_x=0.01,
            beta_y=0.1,
            se_y=0.02,
        )

        self.assertAlmostEqual(
            x["wald_beta"],
            0.2,
            places=10,
        )

        self.assertAlmostEqual(
            x["wald_se_first_order"],
            0.04,
            places=10,
        )

    def test_kidney_direction(self):
        self.assertEqual(
            mod.kidney_direction("eGFRcrea", 0.1),
            "favorable",
        )

        self.assertEqual(
            mod.kidney_direction("BUN", -0.1),
            "favorable",
        )


if __name__ == "__main__":
    unittest.main()
