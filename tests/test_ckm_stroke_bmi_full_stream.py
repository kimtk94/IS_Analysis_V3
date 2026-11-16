"""Offline bounded source integrity checks. Do not touch the network."""
import gzip
import hashlib
import io
import pathlib
import tempfile
import unittest

from scripts.ckm_stroke.audit_bmi_full_stream import (
    EXPECTED, SOURCE, BoundedHashReader, audit_one, scan,
)


class FakeResponse(io.BytesIO):
    def __init__(self, data, url=SOURCE, length=None, status=200):
        super().__init__(data)
        self.url = url
        self.status = status
        self.headers = {"Content-Length": str(len(data) if length is None else length)}

    def geturl(self):
        return self.url


def fixture(*rows, header=EXPECTED):
    return gzip.compress(("\t".join(header) + "\n" +
                          "\n".join(rows) + "\n").encode("utf-8"))


VALID = (
    "7\t92383888\trs10\tA\tC\t0.06431\t0.0013\t0.0042\t7.5e-01\t598895",
    "12\t126890980\trs1000000\tA\tG\t0.2219\t0.0001\t0.0021\t9.6e-01\t689928",
)


class FullBmiOfflineTests(unittest.TestCase):
    def test_clean_gzip_stream_hashes_every_byte(self):
        raw = fixture(*VALID)
        result = scan(FakeResponse(raw))
        self.assertEqual(result["source_status"], "PASS_STREAM_STRUCTURE_ONLY")
        self.assertEqual(result["rows_total"], 2)
        self.assertEqual(result["valid_autosome_rows"], 2)
        self.assertEqual(result["invalid_rows"], 0)
        self.assertEqual(result["source_sha256"], hashlib.sha256(raw).hexdigest())
        self.assertEqual(result["full_compressed_bytes_read"], len(raw))
        self.assertFalse(result["ready_for_ldsc"])
        self.assertTrue(result["no_original_publisher_checksum"])

    def test_reject_broken_header(self):
        raw = fixture(*VALID, header=["SNP", "P"])
        with self.assertRaisesRegex(ValueError, "field names"):
            scan(FakeResponse(raw))

    def test_fails_if_content_length_incomplete(self):
        raw = fixture(*VALID)
        with self.assertRaisesRegex(ValueError, "Content-Length"):
            scan(FakeResponse(raw, length=len(raw) + 1))

    def test_enforce_compressed_size_cap(self):
        raw = fixture(*VALID)
        with self.assertRaisesRegex(ValueError, "cap"):
            scan(FakeResponse(raw), max_bytes=20)

    def test_detect_bad_se_and_alleles(self):
        se_bad = VALID[0].replace("\t0.0042\t", "\t0\t")
        allele_bad = VALID[1].replace("\tA\tG\t", "\tA\tA\t")
        result = scan(FakeResponse(fixture(se_bad, allele_bad, VALID[0])))
        self.assertEqual(result["source_status"], "HOLD_INVALID_RECORDS")
        self.assertEqual(result["invalid_rows"], 2)
        self.assertIn("SE_nonpositive", result["errors_by_type"])
        self.assertIn("alleles_invalid", result["errors_by_type"])

    def test_underflow_p_zero_allowed(self):
        row = VALID[0].replace("\t7.5e-01\t", "\t0\t")
        result = scan(FakeResponse(fixture(row)))
        self.assertEqual(result["source_status"], "PASS_STREAM_STRUCTURE_ONLY")
        self.assertEqual(result["p_underflow_zero_rows"], 1)

    def test_non_autosome_excluded_from_ldsc_autosome_count(self):
        row = VALID[0].replace("7\t92383888", "23\t92383888")
        result = scan(FakeResponse(fixture(row, VALID[1])))
        self.assertEqual(result["non_autosomal_rows"], 1)
        self.assertEqual(result["valid_autosome_rows"], 1)

    def test_truncated_gzip_errors(self):
        raw = fixture(*VALID)
        with self.assertRaises((EOFError, OSError, ValueError)):
            scan(FakeResponse(raw[:-8]))

    def test_pinned_manifest_unapproved(self):
        from scripts.ckm_stroke.audit_bmi_full_stream import pinned_source
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "source.tsv"
            text = ("trait\tstudy\texact_original_url\tgenome_build\tapproved_as_ldsc_input\n" +
                    "BMI\tGIANT_Yengo2018_corrected\t" + SOURCE + "\thg19_GRCh37\tNO\n")
            path.write_text(text, encoding="utf-8")
            result = pinned_source(path)
            self.assertEqual(result["approved_as_ldsc_input"], "NO")
            path.write_text(text.replace("\tNO\n", "\tYES\n"), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "contract"):
                pinned_source(path)

    def test_authored_audit_fails_closed_on_network(self):
        from scripts.ckm_stroke.audit_bmi_full_stream import SOURCE
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "source.tsv"
            path.write_text(
                "trait\tstudy\texact_original_url\tgenome_build\tapproved_as_ldsc_input\n" +
                "BMI\tGIANT_Yengo2018_corrected\t" + SOURCE + "\thg19_GRCh37\tNO\n",
                encoding="utf-8",
            )
            def opener(req, timeout):
                raise ConnectionError("network denied")
            result = audit_one(path, opener=opener)
            self.assertEqual(result["result"]["source_status"], "HOLD_SOURCE_ERROR")
            self.assertFalse(result["ready_for_ldsc"])


if __name__ == "__main__":
    unittest.main()
