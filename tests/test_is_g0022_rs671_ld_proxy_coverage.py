"""Fixture tests: signed EAS LD proxy coverage vs allele-harmonized EUR source."""
import csv
import struct
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/is"))
from audit_g0022_rs671_ld_proxy_ancestry_coverage import audit, RS671, read_ld_row


def tsv(path, rows):
    with path.open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class TestR671LDProxy(unittest.TestCase):
    def make_inputs(self, root):
        variants = root / "variants.tsv"
        coverage = root / "coverage.tsv"
        ld = root / "signed_ld.bin"
        ids = [RS671, "12:112168009:G:A", "12:112468206:C:T", "12:113111111:A:G"]
        tsv(variants, [{"variant_id": v, "status": "ALT_SIGNED_CONFIRMED",
                        "reference_alt_af": "0.2"} for v in ids])
        tsv(coverage, [{"variant_id": v,
                        "eur_source_coverage": "EXACT_ALLELE_PAIR_MATCHED" if i == 3 else "POSITION_ABSENT_EUR_SOURCE",
                        "eur_gwas_alt_eaf": "0.2" if i == 3 else "",
                        "eur_p": "0.5" if i == 3 else ""} for i, v in enumerate(ids)])
        row = (1, .99, .82, .1)
        matrix = [row, (.99, 1, .81, .11), (.82, .81, 1, .10), (.1, .11, .1, 1)]
        with ld.open("wb") as f:
            f.write(struct.pack("<16d", *(x for r in matrix for x in r)))
        return variants, coverage, ld

    def test_proxy_threshold_and_no_eur_ld_inference(self):
        with tempfile.TemporaryDirectory() as t:
            v, c, ld = self.make_inputs(Path(t))
            rows, j = audit(v, c, ld)
            self.assertEqual(j["thresholds"]["0.8"]["eas_ld_variants"], 2)
            self.assertEqual(j["thresholds"]["0.8"]["eur_source_exact_matched"], 0)
            self.assertEqual(j["thresholds"]["0.5"]["eas_ld_variants"], 3)
            self.assertEqual(j["thresholds"]["0.5"]["eur_source_exact_matched"], 0)
            self.assertEqual(j["thresholds"]["0.1"]["eur_source_exact_matched"], 0)
            self.assertFalse(j["multi_ancestry_causal_gain_demonstrated"])
            self.assertEqual(len(rows), 4)

    def test_dimension_mismatch_must_fail(self):
        with tempfile.TemporaryDirectory() as t:
            v, c, ld = self.make_inputs(Path(t))
            with ld.open("ab") as f:
                f.write(b"xyz")
            with self.assertRaises(ValueError):
                audit(v, c, ld)

    def test_bad_ld_diagonal_must_fail(self):
        with tempfile.TemporaryDirectory() as t:
            v, c, ld = self.make_inputs(Path(t))
            with ld.open("r+b") as f:
                f.write(struct.pack("<d", 0.7))
            with self.assertRaises(ValueError):
                audit(v, c, ld)


if __name__ == "__main__":
    unittest.main()
