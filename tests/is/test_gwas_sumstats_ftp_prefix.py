#!/usr/bin/env python3
"""Network-free unit tests for bounded GIGASTROKE source discovery."""
import gzip
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "is"))
from audit_gwas_sumstats_ftp_prefix import (  # noqa: E402
    accession_directory, extract_directory_names, choose_data_file,
    classify_filename, sniff_text, inspect_table_header
)


class SourceAuditTest(unittest.TestCase):
    def test_ftp_bin_boundary(self):
        self.assertEqual(accession_directory("GCST90104540").split("/")[-3],
                         "GCST90104001-GCST90105000")
        self.assertEqual(accession_directory("GCST90104545").split("/")[-2], "GCST90104545")

    def test_invalid_accession(self):
        with self.assertRaises(ValueError):
            accession_directory("../secrets")

    def test_listing_scope_does_not_follow_offsite(self):
        url = accession_directory("GCST90104540")
        html = '<a href="GCST90104540.tsv.gz">safe</a>' + (
            '<a href="https://evil.example/GCST90104540.tsv.gz">bad</a>'
        ) + '<a href="../other.txt.gz">outside</a>'
        self.assertEqual(extract_directory_names(html, url), ["GCST90104540.tsv.gz"])

    def test_select_original_study_data(self):
        names = ["harmonised/", "GCST90104540.tsv.gz", "GCST90104540.tsv.gz-meta.yaml"]
        self.assertEqual(choose_data_file(names, "GCST90104540"), "GCST90104540.tsv.gz")

    def test_yaml_is_not_data(self):
        self.assertEqual(classify_filename("GCST90104540.tsv.gz-meta.yaml"), "META")

    def test_short_gzip_prefix_is_decodable(self):
        text = "chromosome\tbase_pair_location\tp_value\tbeta\tstandard_error\teffect_allele\tother_allele\n1\t2345\t0.02\t0.1\t0.1\tA\tC\n"
        gz = gzip.compress(text.encode("utf-8"))
        decoded, err = sniff_text(gz, "file.tsv.gz")
        self.assertIsNone(err)
        self.assertIn("base_pair_location", decoded)
        self.assertEqual(inspect_table_header(decoded)["header_status"], "HEADER_CORE_FIELDS_PRESENT")

    def test_missing_effect_fields_are_review_not_valid(self):
        header = "chr\tpos\tp_value\n1\t123\t0.03\n"
        self.assertEqual(inspect_table_header(header)["header_status"], "HEADER_REVIEW")

    def test_no_input_marks_unavailable(self):
        self.assertEqual(inspect_table_header("")["header_status"], "UNAVAILABLE")

    def test_bad_gzip_does_not_pass(self):
        decoded, err = sniff_text(b"not actually gzip", "file.tsv.gz")
        self.assertIsNone(decoded)
        self.assertIn("INVALID_GZIP_PREFIX", err)


if __name__ == "__main__":
    unittest.main()
