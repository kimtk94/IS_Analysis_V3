"""Offline fail-closed tests for pinned Zenodo EUR LD Score reference mirror."""
import gzip
import hashlib
import io
import json
import tarfile
import unittest
from unittest.mock import patch

from scripts.ckm_stroke.audit_eur_ldscore_zenodo import (
    FILE_URL, RECORD_API, EXPECTED_MD5, read_record, safe_member, scan_tar,
)


class Response(io.BytesIO):
    def __init__(self, raw, url=FILE_URL, content_length=None):
        super().__init__(raw)
        self._url = url
        self.headers = {"Content-Length": str(len(raw) if content_length is None else content_length)}
        self.status = 200

    def geturl(self):
        return self._url


def fixture(omit_chrom=None, include_M=True):
    output = io.BytesIO()
    with tarfile.open(fileobj=output, mode="w:gz") as tar:
        for chrom in range(1, 23):
            if chrom == omit_chrom:
                continue
            prefix = f"eur_w_ld_chr/{chrom}.l2."
            text = "CHR SNP BP L2\n" + "".join(
                f"{chrom} rs{i + 1} {i + 2} 1.2\n" for i in range(1001)
            )
            parts = [
                (prefix + "ldscore.gz", gzip.compress(text.encode())),
                (prefix + "M", b"100\n"),
                (prefix + "M_5_50", b"70\n"),
            ]
            if not include_M:
                parts = [item for item in parts if not item[0].endswith(".l2.M")]
            for name, blob in parts:
                info = tarfile.TarInfo(name)
                info.size = len(blob)
                tar.addfile(info, io.BytesIO(blob))
    return output.getvalue()


class ReferenceTests(unittest.TestCase):
    def test_full_22_chromosomes_passes_archive_only(self):
        raw = fixture()
        with patch("scripts.ckm_stroke.audit_eur_ldscore_zenodo.SIZE", len(raw)), patch(
            "scripts.ckm_stroke.audit_eur_ldscore_zenodo.EXPECTED_MD5",
            hashlib.md5(raw).hexdigest(),
        ):
            got = scan_tar(Response(raw))
        self.assertEqual(got["technical_status"], "PASS_ZENODO_ARCHIVE_AND_22_CHR_LD_REFERENCE")
        self.assertEqual(len(got["chromosome_sets"]["ldscore"]), 22)
        self.assertEqual(got["ldscore_headers"]["22"]["row_count_lower_bound"], 1001)
        self.assertEqual(got["archive_sha256"], hashlib.sha256(raw).hexdigest())
        self.assertFalse(got["ready_for_ldsc"])

    def test_optional_M_absent_still_passes_standard_ldsc_mode(self):
        raw = fixture(include_M=False)
        with patch("scripts.ckm_stroke.audit_eur_ldscore_zenodo.SIZE", len(raw)), patch(
            "scripts.ckm_stroke.audit_eur_ldscore_zenodo.EXPECTED_MD5",
            hashlib.md5(raw).hexdigest(),
        ):
            result = scan_tar(Response(raw))
        self.assertEqual(result["technical_status"], "PASS_ZENODO_ARCHIVE_AND_22_CHR_LD_REFERENCE")
        self.assertTrue(result["ldsc_default_M_5_50_available"])
        self.assertFalse(result["ldsc_not_M_5_50_option_supported"])

    def test_missing_chromosome_is_hold(self):
        raw = fixture(omit_chrom=17)
        with patch("scripts.ckm_stroke.audit_eur_ldscore_zenodo.SIZE", len(raw)), patch(
            "scripts.ckm_stroke.audit_eur_ldscore_zenodo.EXPECTED_MD5",
            hashlib.md5(raw).hexdigest(),
        ):
            got = scan_tar(Response(raw))
        self.assertEqual(got["technical_status"], "HOLD_MISSING_REQUIRED_LD_COMPONENT")
        self.assertEqual(got["missing_chromosomes"]["M"], [17])

    def test_archive_checksum_mismatch_is_error(self):
        raw = fixture()
        with patch("scripts.ckm_stroke.audit_eur_ldscore_zenodo.SIZE", len(raw)):
            with self.assertRaisesRegex(ValueError, "MD5"):
                scan_tar(Response(raw))

    def test_rejects_unsafe_tar_member(self):
        with self.assertRaisesRegex(ValueError, "unsafe"):
            safe_member("eur_w_ld_chr/../passwd")

    def test_rejects_stale_zenodo_record(self):
        item = {
            "id": 8182036,
            "doi": "10.5281/zenodo.8182036",
            "files": [
                {"key": "eur_w_ld_chr.tar.gz", "size": 33357890,
                 "checksum": "md5:" + "0" * 32,
                 "links": {"self": FILE_URL}}
            ],
        }
        def opener(request, timeout):
            self.assertEqual(request.full_url, RECORD_API)
            return Response(json.dumps(item).encode(), RECORD_API)
        with self.assertRaisesRegex(ValueError, "differs"):
            read_record(opener)


if __name__ == "__main__":
    unittest.main()
