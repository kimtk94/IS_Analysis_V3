#!/usr/bin/env python3
"""Offline synthetic MEX tests; no synthetic data are biological evidence."""
import csv
import gzip
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "is"))
from qc_gse331114_mex import (
    SAMPLES, CANDIDATES, SOFT, safe_file_url, parse_soft_samples,
    select_files, load_features, count_barcodes, qc_matrix,
)


class GSE331114MatrixTests(unittest.TestCase):
    def test_eight_predefined_independent_gsm_records(self):
        self.assertEqual(len(SAMPLES), 8)
        self.assertEqual(sum(x[0] == "AIS" for x in SAMPLES.values()), 4)
        self.assertEqual(sum(x[0] == "Control" for x in SAMPLES.values()), 4)
        self.assertEqual(len(CANDIDATES), 15)
        self.assertEqual(len(set(CANDIDATES)), 15)

    def test_only_ncbi_sample_urls(self):
        url = "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9739nnn/GSM9739251/suppl/GSM9739251_AIS_brain1_matrix.mtx.gz"
        self.assertEqual(safe_file_url(url, "GSM9739251").split(":")[0], "https")
        with self.assertRaises(ValueError):
            safe_file_url("https://evil.example/GSM9739251/suppl/x", "GSM9739251")

    def test_sample_file_manifest_check(self):
        sample = {
            "organism": "Homo sapiens",
            "title": "AIS brain",
            "files": [
                "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9739nnn/GSM9739251/suppl/"
                f"GSM9739251_AIS_brain1_{suffix}"
                for suffix in ("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")
            ],
        }
        result = select_files(sample, "GSM9739251", "AIS", 1)
        self.assertEqual(set(result), {"matrix", "features", "barcodes"})

    def test_missing_matrix_is_fatal(self):
        sample = {"organism": "Homo sapiens", "files": []}
        with self.assertRaises(ValueError):
            select_files(sample, "GSM9739251", "AIS", 1)

    def test_gzip_sof_sample_metadata(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "family.soft.gz"
            with gzip.open(path, "wt") as handle:
                handle.write(
                    "^SERIES = GSE331114\n"
                    "^SAMPLE = GSM9739251\n"
                    "!Sample_title = AIS brain1\n"
                    "!Sample_organism_ch1 = Homo sapiens\n"
                    "!Sample_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/suppl/s1\n"
                    "^SAMPLE = GSM9739255\n"
                    "!Sample_title = Control brain1\n"
                )
            found = parse_soft_samples(path)
            self.assertEqual(len(found), 2)
            self.assertEqual(len(found["GSM9739251"]["files"]), 1)

    def test_synthetic_complete_matrix(self):
        with tempfile.TemporaryDirectory() as folder:
            directory = Path(folder)
            genes_path = directory / "features.tsv.gz"
            cells_path = directory / "barcodes.tsv.gz"
            matrix_path = directory / "matrix.mtx.gz"
            with gzip.open(genes_path, "wt") as out:
                for i in range(1100):
                    symbol = "ABO" if i == 0 else ("MT-CO1" if i == 1 else f"GENE_{i}")
                    out.write(f"ENSG{i:011d}\t{symbol}\tGene Expression\n")
            with gzip.open(cells_path, "wt") as out:
                for i in range(120):
                    out.write(f"AAAG{i:04d}-1\n")
            # Cell 1: 250 detected genes, UMI 750, mito=0; all other cells fail min gene cutoff.
            with gzip.open(matrix_path, "wt") as out:
                out.write("%%MatrixMarket matrix coordinate integer general\n%\n1100 120 250\n")
                for i in range(250):
                    out.write(f"{i+1} 1 3\n")
            genes = load_features(genes_path)
            ncells = count_barcodes(cells_path)
            summary, panel = qc_matrix(matrix_path, genes, ncells)
            self.assertEqual(summary["source_matrix"], {"features": 1100, "barcodes": 120, "nnz": 250})
            self.assertEqual(summary["barcodes_qc_pass"], 1)
            self.assertEqual(summary["analysis_readiness"], "HOLD_TOO_FEW_BARCODES")
            abo = [x for x in panel if x["gene_symbol"] == "ABO"][0]
            self.assertEqual(abo["total_umi_qc_barcodes"], 3)
            self.assertAlmostEqual(abo["sample_pseudobulk_cpm"], 4000)

    def test_mismatched_matrix_dimensions_fail(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "broken.mtx.gz"
            with gzip.open(path, "wt") as handle:
                handle.write("%%MatrixMarket matrix coordinate integer general\n1000 120 200\n")
            features = [(f"ENSG{i:011d}", f"GENE{i}") for i in range(1100)]
            with self.assertRaises(ValueError):
                qc_matrix(path, features, 120)


if __name__ == "__main__":
    unittest.main()
