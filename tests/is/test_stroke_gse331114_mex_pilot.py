#!/usr/bin/env python3
"""Real MEX audit uses open public sources; tests synthetic matrices only."""
import gzip
import sys
import tempfile
import unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "is"))
from audit_stroke_gse331114_mex_pilot import (
    classify_file, group_files, count_rows, inspect_mtx
)


class StrokeMEXPilotTests(unittest.TestCase):
    def test_standard_triplet_filenames(self):
        self.assertEqual(classify_file("https://ftp.ncbi.nlm.nih.gov/a/matrix.mtx.gz"), "matrix")
        self.assertEqual(classify_file("https://ftp.ncbi.nlm.nih.gov/a/barcodes.tsv.gz"), "barcodes")
        self.assertEqual(classify_file("https://ftp.ncbi.nlm.nih.gov/a/features.tsv.gz"), "features")

    def test_gsm_prefix_10x_filenames(self):
        self.assertEqual(classify_file("ftp://ftp.ncbi.nlm.nih.gov/a/GSM1_matrix.mtx.gz"), "matrix")
        self.assertEqual(classify_file("ftp://ftp.ncbi.nlm.nih.gov/a/GSM1_barcodes.tsv.gz"), "barcodes")
        self.assertEqual(classify_file("ftp://ftp.ncbi.nlm.nih.gov/a/GSM1_features.tsv.gz"), "features")

    def test_group_file_missing_is_explicit(self):
        result = group_files({"files": ["ftp://ftp.ncbi.nlm.nih.gov/matrix.mtx.gz"]})
        self.assertEqual(set(result), {"matrix"})

    def test_integer_matrix_content_and_dimensions(self):
        with tempfile.TemporaryDirectory() as d:
            matrix = Path(d) / "matrix.mtx.gz"
            features = Path(d) / "features.tsv.gz"
            barcodes = Path(d) / "barcodes.tsv.gz"
            with gzip.open(matrix, "wt") as f:
                f.write("%%MatrixMarket matrix coordinate integer general\n%\n2 3 2\n1 1 4\n2 3 2\n")
            with gzip.open(features, "wt") as f:
                f.write("ENSG1\tGENE1\nENSG2\tGENE2\n")
            with gzip.open(barcodes, "wt") as f:
                f.write("AA\nBB\nCC\n")
            summary = inspect_mtx(matrix)
            self.assertEqual(summary["genes_rows"], count_rows(features, min_fields=2))
            self.assertEqual(summary["barcodes_columns"], count_rows(barcodes))
            self.assertEqual(summary["nonzero_entries"], 2)
            self.assertEqual(summary["sum_matrix_entries"], 6)

    def test_invalid_matrix_fails(self):
        with tempfile.TemporaryDirectory() as d:
            matrix = Path(d) / "bad.mtx.gz"
            with gzip.open(matrix, "wt") as f:
                f.write("%%MatrixMarket matrix coordinate integer general\n2 2 1\n3 1 4\n")
            with self.assertRaises(ValueError):
                inspect_mtx(matrix)


if __name__ == "__main__":
    unittest.main()
