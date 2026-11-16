#!/usr/bin/env python3
"""Offline regression checks of full original GIGASTROKE data QC."""
import csv
import gzip
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "is"))
from audit_gwas_full_integrity import (  # noqa: E402
    full_url, read_numeric, validate_row, scan_source_gzip
)


def row(**extra):
    values = {
        "chromosome": "1", "base_pair_location": "1010",
        "effect_allele_frequency": "0.3", "beta": "0.01",
        "standard_error": "0.01", "p_value": "0.05",
        "effect_allele": "A", "other_allele": "C"
    }
    values.update(extra)
    return values


class SourceFullQcTest(unittest.TestCase):
    def test_expected_url(self):
        self.assertEqual(
            full_url("GCST90104540"),
            "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
            "GCST90104001-GCST90105000/GCST90104540/"
            "GCST90104540_buildGRCh37.tsv.gz"
        )

    def test_valid_row(self):
        self.assertIsNone(validate_row(row()))

    def test_zero_p_is_allowed_as_possible_underflow(self):
        self.assertIsNone(validate_row(row(p_value="0")))

    def test_bad_p_fails(self):
        self.assertEqual(validate_row(row(p_value="-0.03")), "INVALID_P_VALUE")

    def test_bad_se_fails(self):
        self.assertEqual(validate_row(row(standard_error="0")), "INVALID_SE")

    def test_non_snp_fails(self):
        self.assertEqual(validate_row(row(effect_allele="AT")), "INVALID_ALLELES")

    def test_eaf_outside_range_fails(self):
        self.assertEqual(validate_row(row(effect_allele_frequency="1.2")), "INVALID_FREQUENCY")

    def test_x_chromosome_is_excluded_not_malformed(self):
        self.assertEqual(validate_row(row(chromosome="X")), "CHROMOSOME_OUTSIDE_AUTOSOMES")

    def test_scan_complete_gzip_with_no_sample_count_or_rsid(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.tsv.gz"
            fields = list(row())
            with gzip.open(path, "wt", newline="", encoding="utf-8") as out:
                writer = csv.DictWriter(out, delimiter="\t", fieldnames=fields)
                writer.writeheader()
                writer.writerow(row())
                writer.writerow(row(chromosome="X"))
                writer.writerow(row(standard_error="0"))
            result = scan_source_gzip(path)
            self.assertEqual(result["full_gzip_crc_integrity"], "PASS")
            self.assertEqual(result["total_variant_rows"], 3)
            self.assertEqual(result["valid_autosomal_snp_rows"], 1)
            self.assertEqual(result["excluded_nonautosomal_rows_not_invalid"], 1)
            self.assertEqual(result["invalid_row_counts"]["INVALID_SE"], 1)
            self.assertFalse(result["rsid_column_present"])
            self.assertFalse(result["per_variant_sample_size_column_present"])
            self.assertEqual(result["row_qc_status"], "STRUCTURE_QC_REVIEW")

    def test_bad_gzip_is_not_a_pass(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.tsv.gz"
            path.write_bytes(b"not-gzip")
            with self.assertRaises((OSError, EOFError)):
                scan_source_gzip(path)


if __name__ == "__main__":
    unittest.main()
