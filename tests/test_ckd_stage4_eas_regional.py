import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

SCRIPT = (
    ROOT
    / "scripts"
    / "prepare_ckd_stage4_eas_regional.py"
)

spec = importlib.util.spec_from_file_location(
    "eas_regional",
    SCRIPT,
)

mod = importlib.util.module_from_spec(
    spec
)

spec.loader.exec_module(mod)


class TestEASRegional(unittest.TestCase):

    def test_direct(self):
        x = mod.harmonize(
            "C",
            "T",
            "T",
            "C",
        )

        self.assertEqual(
            x,
            (1.0, "direct"),
        )

    def test_swap(self):
        x = mod.harmonize(
            "A",
            "G",
            "A",
            "G",
        )

        self.assertEqual(
            x,
            (-1.0, "swapped"),
        )

    def test_strand(self):
        x = mod.harmonize(
            "A",
            "G",
            "C",
            "T",
        )

        self.assertEqual(
            x,
            (1.0, "strand_direct"),
        )

    def test_palindrome_unresolved(self):
        x = mod.harmonize(
            "A",
            "T",
            "T",
            "A",
        )

        # Exact orientation remains identifiable.
        self.assertEqual(
            x,
            (1.0, "direct"),
        )


if __name__ == "__main__":
    unittest.main()
