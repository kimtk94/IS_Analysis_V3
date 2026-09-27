import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

SCRIPT = (
    ROOT
    / "scripts"
    / "prepare_ckd_stage4_dual_ld.py"
)

spec = importlib.util.spec_from_file_location(
    "dual_ld",
    SCRIPT,
)

mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class TestDualLD(unittest.TestCase):

    def test_subset_window(self):

        rows = [
            {
                "chr37": "1",
                "pos37": "900",
                "snp": "a",
            },
            {
                "chr37": "1",
                "pos37": "1000",
                "snp": "b",
            },
            {
                "chr37": "1",
                "pos37": "1100",
                "snp": "c",
            },
            {
                "chr37": "2",
                "pos37": "1000",
                "snp": "d",
            },
        ]

        out, lo, hi = mod.subset_window(
            rows,
            "1",
            1000,
            100,
        )

        self.assertEqual(lo, 900)
        self.assertEqual(hi, 1100)

        self.assertEqual(
            [x["snp"] for x in out],
            ["a", "b", "c"],
        )

    def test_unique_by_snp(self):

        rows = [
            {
                "snp": "rs1",
                "p_pqtl": "0.1",
            },
            {
                "snp": "rs1",
                "p_pqtl": "0.01",
            },
            {
                "snp": "rs2",
                "p_pqtl": "0.2",
            },
        ]

        out = mod.unique_by_snp(rows)

        self.assertEqual(
            len(out),
            2,
        )

        self.assertEqual(
            out["rs1"]["p_pqtl"],
            "0.01",
        )


if __name__ == "__main__":
    unittest.main()
