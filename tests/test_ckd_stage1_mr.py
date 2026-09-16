import importlib.util
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_ckd_stage1_mr.py"
SPEC = importlib.util.spec_from_file_location("ckd_stage1", SCRIPT)
M = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(M)


class TestCKDStage1MR(unittest.TestCase):
    def exp(self, ref="A", alt="G", af="0.2"):
        return {"ref_allele": ref, "alt_allele": alt, "alt_freq": af}

    def out(self, ea="G", oa="A", eaf="0.2"):
        return {"effect_allele": ea, "other_allele": oa, "eaf": eaf}

    def test_direct(self):
        self.assertEqual(M.harmonize(self.exp(), self.out(), .42, .1), ("direct", 1))

    def test_swapped(self):
        self.assertEqual(
            M.harmonize(self.exp(), self.out("A", "G", "0.8"), .42, .1),
            ("swapped", -1),
        )

    def test_strand(self):
        self.assertEqual(
            M.harmonize(self.exp(), self.out("C", "T", "0.2"), .42, .1),
            ("strand", 1),
        )

    def test_palindrome_frequency_keep(self):
        self.assertEqual(
            M.harmonize(self.exp("A", "T", "0.1"), self.out("T", "A", "0.11"), .42, .1),
            ("palindrome_freq_keep", 1),
        )

    def test_palindrome_frequency_flip(self):
        self.assertEqual(
            M.harmonize(self.exp("A", "T", "0.1"), self.out("T", "A", "0.89"), .42, .1),
            ("palindrome_freq_flip", -1),
        )

    def test_palindrome_without_frequency_is_dropped(self):
        self.assertEqual(
            M.harmonize(self.exp("A", "T", "0.1"), self.out("T", "A", ""), .42, .1),
            ("palindrome_no_frequency", None),
        )

    def test_bh(self):
        q = M.bh_adjust([0.01, 0.04, 0.03, None])
        self.assertAlmostEqual(q[0], 0.03)
        self.assertAlmostEqual(q[1], 0.04)
        self.assertAlmostEqual(q[2], 0.04)
        self.assertIsNone(q[3])

    def test_ivw(self):
        x = [
            {"beta_exposure": 0.2, "beta_outcome": 0.1, "se_outcome": 0.02},
            {"beta_exposure": 0.1, "beta_outcome": 0.05, "se_outcome": 0.01},
        ]
        beta, se, p, q = M.ivw(x)
        self.assertAlmostEqual(beta, 0.5)
        self.assertGreater(se, 0)
        self.assertLess(p, 1e-6)
        self.assertAlmostEqual(q, 0.0)


if __name__ == "__main__":
    unittest.main()
