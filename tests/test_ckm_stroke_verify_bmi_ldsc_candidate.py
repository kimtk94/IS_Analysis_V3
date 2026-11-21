"""Offline tests for independent candidate digest/count/allele validator."""
import gzip
import hashlib
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from scripts.ckm_stroke.verify_bmi_ldsc_candidate import inspect


def provenance(raw):
    return {"technical_status":"PASS_CONSERVATIVE_CANDIDATE_CREATED_NOT_LDSC_APPROVED",
            "scientific_status":"HOLD","ready_for_ldsc":False,
            "counts":{"candidate_rows":len(raw.splitlines())-1},
            "candidate_plaintext_sha256":hashlib.sha256(raw).hexdigest()}


class CandidateVerifierTests(unittest.TestCase):
    def test_valid_source_never_promotes(self):
        raw=b"SNP\tN\tZ\tA1\tA2\nrs1\t600000\t-4\tA\tC\nrs2\t500000\t2.5\tG\tA\n"
        with TemporaryDirectory() as d:
            path=Path(d)/"bmi.sumstats.gz"
            path.write_bytes(gzip.compress(raw))
            result=inspect(path,provenance(raw))
        self.assertEqual(result["candidate_rows_verified"],2)
        self.assertEqual(result["candidate_plaintext_sha256"],hashlib.sha256(raw).hexdigest())
        self.assertFalse(result["ready_for_ldsc"])

    def test_duplicate_rsid_is_rejected(self):
        raw=b"SNP\tN\tZ\tA1\tA2\nrs1\t10\t1\tA\tC\nrs1\t20\t2\tA\tG\n"
        with TemporaryDirectory() as d:
            path=Path(d)/"dup.gz"
            path.write_bytes(gzip.compress(raw))
            with self.assertRaisesRegex(ValueError,"duplicated"):
                inspect(path,provenance(raw))

    def test_palindromic_alleles_are_rejected(self):
        raw=b"SNP\tN\tZ\tA1\tA2\nrs1\t10\t1\tA\tT\n"
        with TemporaryDirectory() as d:
            path=Path(d)/"pal.gz"
            path.write_bytes(gzip.compress(raw))
            with self.assertRaisesRegex(ValueError,"ambiguous"):
                inspect(path,provenance(raw))

    def test_hash_mismatch_rejected(self):
        raw=b"SNP\tN\tZ\tA1\tA2\nrs1\t10\t1\tA\tC\n"
        with TemporaryDirectory() as d:
            path=Path(d)/"hash.gz"
            path.write_bytes(gzip.compress(raw))
            expected=provenance(raw)
            expected["candidate_plaintext_sha256"]="0"*64
            with self.assertRaisesRegex(ValueError,"hash"):
                inspect(path,expected)


if __name__=="__main__":
    unittest.main()
