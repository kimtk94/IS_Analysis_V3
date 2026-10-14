import csv
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "is"))
from inventory_snp_inputs import discover

class TestInputInventory(unittest.TestCase):
    def test_finds_candidates_without_claiming_validity(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d) / "source"
            root.mkdir()
            (root / "bbj_gwas.tsv").write_text("not a verified table\n")
            (root / "unrelated.tsv").write_text("irrelevant\n")
            out = Path(d) / "out"
            result = discover([root], out)
            self.assertEqual(result["possible_inputs"], 1)
            with (out / "IS_POSSIBLE_SNP_INPUT_INVENTORY.tsv").open() as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
            self.assertEqual(rows[0]["eligible_for_coloc"], "NO")
            self.assertEqual(rows[0]["source_verified"], "NO")
            self.assertEqual((root / "bbj_gwas.tsv").read_text(), "not a verified table\n")

    def test_scan_cap(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "root"
            p.mkdir()
            (p / "a").write_text("a")
            (p / "b").write_text("b")
            with self.assertRaisesRegex(ValueError, "cap"):
                discover([p], Path(d) / "out", max_files=1)

if __name__ == "__main__":
    unittest.main()
