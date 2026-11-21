"""Offline tests for user-authorized GWAS archive intake (no web, no real GWAS)."""
import gzip
import io
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch
import zipfile

from scripts.ckm_stroke.audit_authorized_sbp_t2d_file import (
    PUBLISHED_T2D_MD5, check, digest_file, examine_gzip_header, examine_zip,
)


def write_ledger(path: Path):
    path.write_text(
        "trait\tstudy\tapproved_as_ldsc_input\tverified_landing_url\n"
        "SBP\tICBP_Evangelou2018_original_GRASP\tNO\thttps://example.org/source\n"
        "SBP\tIEU_OpenGWAS_ieu-b-38\tNO\thttps://example.org/secondary\n"
        "T2D\tT2DGGI_Suzuki2024_EUR\tNO\thttps://example.org/terms\n",
        encoding="utf-8",
    )


class TestAuthorizedInput(unittest.TestCase):
    def test_sbp_gzip_reports_sha_without_publisher_checksum(self):
        with TemporaryDirectory() as d:
            base = Path(d)
            f = base / "SBP.txt.gz"
            f.write_bytes(gzip.compress(b"CHR\tPOS\tBETA\tSE\n1\t10\t0.2\t0.01\n"))
            ledger = base / "ledger.tsv"
            write_ledger(ledger)
            result = check("SBP", "ICBP_Evangelou2018_original_GRASP", f, ledger, repo_root=base / "repo")
            self.assertEqual(result["source_status"], "HOLD_NO_PUBLISHER_CHECKSUM")
            self.assertFalse(result["ready_for_ldsc"])
            self.assertFalse(result["raw_local_path_recorded"])
            self.assertEqual(result["container_schema"]["type"], "TSV")

    def test_sbp_vcf_header_read(self):
        with TemporaryDirectory() as d:
            f = Path(d) / "ieu-b-38.vcf.gz"
            f.write_bytes(gzip.compress(b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"))
            self.assertEqual(examine_gzip_header(f, is_vcf=True)["type"], "VCF")

    def test_zip_source_and_checksum_mismatch_hold(self):
        with TemporaryDirectory() as d:
            base = Path(d)
            f = base / "t2d_2024.zip"
            with zipfile.ZipFile(f, "w", zipfile.ZIP_DEFLATED) as z:
                z.writestr("EUR_summary.txt", "SNP\tP\nrs1\t0.1\n")
            ledger = base / "ledger.tsv"
            write_ledger(ledger)
            result = check("T2D", "T2DGGI_Suzuki2024_EUR", f, ledger, repo_root=base / "repo")
            self.assertEqual(result["source_status"], "HOLD_PUBLISHED_MD5_MISMATCH")
            self.assertFalse(result["container_schema"]["zip_members_extracted"])

    def test_publisher_checksum_match_can_only_pass_archive(self):
        with TemporaryDirectory() as d:
            base = Path(d)
            f = base / "t2d_2024.zip"
            with zipfile.ZipFile(f, "w") as z:
                z.writestr("data.txt", "SNP\tP\nrs1\t0.1\n")
            ledger = base / "ledger.tsv"
            write_ledger(ledger)
            fake = {"bytes": f.stat().st_size, "md5": PUBLISHED_T2D_MD5, "sha256": "a" * 64}
            with patch("scripts.ckm_stroke.audit_authorized_sbp_t2d_file.digest_file", return_value=fake):
                result = check("T2D", "T2DGGI_Suzuki2024_EUR", f, ledger, repo_root=base / "repo")
            self.assertEqual(result["source_status"], "PASS_ARCHIVE_MD5_ONLY")
            self.assertFalse(result["ready_for_ldsc"])

    def test_reject_archive_inside_repository(self):
        with TemporaryDirectory() as d:
            base = Path(d)
            repo = base / "repo"
            repo.mkdir()
            f = repo / "source.gz"
            f.write_bytes(gzip.compress(b"CHR\tPOS\tBETA\tSE\n"))
            ledger = base / "ledger.tsv"
            write_ledger(ledger)
            with self.assertRaisesRegex(ValueError, "outside Git"):
                check("SBP", "ICBP_Evangelou2018_original_GRASP", f, ledger, repo_root=repo)

    def test_unsafe_zip_member_rejected(self):
        with TemporaryDirectory() as d:
            f = Path(d) / "file.zip"
            with zipfile.ZipFile(f, "w") as z:
                z.writestr("../outside", b"hello")
            with self.assertRaisesRegex(ValueError, "unsafe"):
                examine_zip(f)


if __name__ == "__main__":
    unittest.main()
