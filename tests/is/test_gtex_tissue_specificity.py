#!/usr/bin/env python3
"""Synthetic tests only; true gene expression is evaluated in separate CI task."""
import csv
import gzip
import sys
import tempfile
import unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "is"))
from audit_gtex_tissue_specificity import (  # noqa: E402
    tau, parse_gct, CANDIDATES, TISSUES, SOURCE
)


class GtexSpecificityTests(unittest.TestCase):
    def test_source_is_public_gtex_v8_median(self):
        self.assertIn("GTEx_Analysis_2017-06-05_v8", SOURCE)
        self.assertTrue(SOURCE.endswith("_gene_median_tpm.gct.gz"))

    def test_prespecified_genes_no_duplicates(self):
        self.assertEqual(len(CANDIDATES), len(set(CANDIDATES)))
        self.assertEqual(len(CANDIDATES), 15)

    def test_tau_is_descriptive(self):
        self.assertAlmostEqual(tau([0, 0, 10]), 1)
        self.assertAlmostEqual(tau([10, 10, 10]), 0)
        self.assertIsNone(tau([0, 0, 0]))

    def test_parse_gct_tissue_matrix(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "synthetic.gct.gz"
            columns = list(TISSUES.values()) + [f"other_{i}" for i in range(35)]
            with gzip.open(path, "wt", encoding="utf-8") as f:
                f.write("#1.2\n2 44\n")
                writer = csv.writer(f, delimiter="\t")
                writer.writerow(["Name", "Description", *columns])
                writer.writerow(["ENSG00000000001.1", "ABO", *([1] * 44)])
                writer.writerow(["ENSG00000000002.1", "IL6RA", *([2] * 44)])
            cols, matches, n = parse_gct(path)
            self.assertEqual(len(cols), 44)
            self.assertEqual(n, 2)
            self.assertEqual(len(matches["ABO"]), 1)
            self.assertEqual(len(matches["IL6RA"]), 1)


if __name__ == "__main__":
    unittest.main()
