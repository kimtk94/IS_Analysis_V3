import csv
import gzip
import importlib.util
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts/is/audit_eas5_original_archive_gate.py"
spec = importlib.util.spec_from_file_location("is_eas5_gate", MODULE_PATH)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

HEADERS = ["mapped_genomic_variant_GRCh37", "BBJ_beta_ALT",
           "BBJ_p", "BBJ_allele_exact_Fasta37"]

def write_table(file, rows, *, headers=HEADERS):
    with gzip.open(file, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

def entry(key="1:100:A:G", beta="0.2", p="1e-8", exact="True"):
    return dict(zip(HEADERS, [key,beta,p,exact]))

class EAS5OriginalArchiveTests(unittest.TestCase):
    def test_missing_archive_is_scientific_hold(self):
        result = mod.audit()
        self.assertEqual(result["file_gate"], "HOLD_ARCHIVE_NOT_AVAILABLE_ON_CI_RUNNER")
        self.assertFalse(result["archive_bytes_available"])
        self.assertIsNone(result["archive_sha256"])
        self.assertFalse(result["coloc_executed"])

    def test_real_file_calculates_sha_without_trusting_it(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "source.tsv.gz"
            write_table(path, [entry()])
            result = mod.audit(path)
            self.assertEqual(result["archive_sha256"], mod.sha256_file(path))
            self.assertEqual(result["bounded_qc"]["distinct_source_eligible_variant_keys"], 1)
            self.assertEqual(result["file_gate"], "HOLD_NO_TRUSTED_EXPECTED_SHA256")
            self.assertFalse(result["provided_expected_sha256_matched"])
            self.assertFalse(result["expected_sha256_independently_authenticated"])

    def test_wrong_sha_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "source.tsv.gz"
            write_table(path, [entry()])
            report = mod.audit(path, expected_sha256="0"*64)
            self.assertEqual(report["file_gate"], "HOLD_SOURCE_SHA256_MISMATCH")
            self.assertFalse(report["study_sample_size_adopted"])

    def test_matching_sha_does_not_override_key_count(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "source.tsv.gz"
            write_table(path, [entry()])
            report = mod.audit(path, expected_sha256=mod.sha256_file(path))
            self.assertEqual(report["file_gate"], "HOLD_SOURCE_ELIGIBLE_VARIANT_COUNT_MISMATCH")
            self.assertFalse(report["BBJ_effect_allele_independently_verified"])

    def test_conflicting_duplicate_is_a_blocker(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "source.tsv.gz"
            write_table(path, [entry(), entry(beta="0.4")])
            report = mod.audit(path, expected_sha256=mod.sha256_file(path))
            self.assertEqual(report["file_gate"], "HOLD_CONFLICTING_ARCHIVE_DUPLICATES")
            self.assertEqual(report["bounded_qc"]["conflicting_duplicate_keys"], 1)

    def test_missing_header_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "source.tsv.gz"
            write_table(path, [{"mapped_genomic_variant_GRCh37":"1:100:A:G"}],
                        headers=["mapped_genomic_variant_GRCh37"])
            self.assertEqual(mod.audit(path)["file_gate"], "HOLD_REQUIRED_COLUMNS_MISSING")

    def test_corrupt_gzip_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "source.tsv.gz"
            path.write_bytes(b"not gzip")
            report = mod.audit(path)
            self.assertEqual(report["file_gate"], "HOLD_ARCHIVE_PARSE_OR_GZIP_ERROR")
            self.assertEqual(report["scientific_status"], "HOLD")

    def test_unverified_allele_does_not_enter_candidate_keys(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "source.tsv.gz"
            write_table(path, [entry(exact="False"),entry(key="1:101:C:T")])
            report = mod.audit(path)
            self.assertEqual(report["bounded_qc"]["distinct_source_eligible_variant_keys"],1)
            self.assertFalse(report["original_build_and_ref_alt_independently_verified"])

if __name__ == "__main__":
    unittest.main()
