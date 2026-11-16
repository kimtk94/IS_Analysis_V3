"""Offline regression for fail-closed GIGASTROKE original-byte audit."""
import gzip
import hashlib
import io
import unittest

from scripts.ckm_stroke.verify_gigastroke_raw import (
    audit_one, check_sidecars, digest_stream, source_urls,
)

ROW = {
    "phenotype": "LAA", "ancestry": "EUR",
    "accession": "GCST90104542", "reported_trait": "Large artery stroke",
    "build": "GRCh37", "sample_total": "1241207",
    "published_raw_md5": "", "verification": "METADATA_VERIFIED_RAW_HASH_UNCHECKED",
    "metadata_url": "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
                    "GCST90104001-GCST90105000/GCST90104542/"
                    "GCST90104542_buildGRCh37.tsv.gz-meta.yaml",
}


class Response(io.BytesIO):
    def __init__(self, body, length=None):
        super().__init__(body)
        self.headers = {"Content-Length": str(len(body) if length is None else length)}
        self.status = 200


class TestIntegrity(unittest.TestCase):
    def setUp(self):
        self.data = gzip.compress(b"chromosome\tbase_pair_location\tbeta\n1\t12\t0.4\n")
        self.row = {**ROW, "published_raw_md5": hashlib.md5(self.data).hexdigest()}

    def metadata(self):
        r = self.row
        return ("gwas_id: GCST90104542\ntrait_description:\n  - Large artery stroke\n"
                "genome_assembly: GRCh37\nsample_size: 1241207\n"
                "data_file_name: GCST90104542_buildGRCh37.tsv.gz\n"
                "data_file_md5sum: " + r["published_raw_md5"] + "\n")

    def sidecar(self):
        return self.row["published_raw_md5"] + " GCST90104542_buildGRCh37.tsv.gz\n"

    def fetcher(self, request, timeout=30):
        u = request.full_url
        if u.endswith("-meta.yaml"):
            return Response(self.metadata().encode())
        if u.endswith("md5sum.txt"):
            return Response(self.sidecar().encode())
        return Response(self.data)

    def test_source_urls_are_strictly_pinned(self):
        urls = source_urls(self.row)
        self.assertEqual(len(urls), 3)
        bad = {**self.row, "metadata_url": "https://example.com/other"}
        with self.assertRaises(ValueError):
            source_urls(bad)

    def test_full_valid_stream_passes(self):
        result = audit_one(self.row, fetcher=self.fetcher)
        self.assertEqual(result["status"], "PASS_RAW_SOURCE_INTEGRITY")
        self.assertEqual(result["scientific_status"], "HOLD")
        self.assertTrue(result["gzip_magic_verified"])

    def test_incomplete_http_stream_fails(self):
        with self.assertRaisesRegex(ValueError, "incomplete"):
            digest_stream(Response(self.data, len(self.data) + 2), max_bytes=10000)

    def test_byte_limit_fails(self):
        with self.assertRaisesRegex(ValueError, "cap"):
            digest_stream(Response(self.data), max_bytes=10)

    def test_sidecar_tampering_fails(self):
        with self.assertRaisesRegex(ValueError, "md5sum"):
            check_sidecars(self.row, self.metadata(), "0"*32 + " GCST90104542_buildGRCh37.tsv.gz")

    def test_non_gzip_fails(self):
        with self.assertRaisesRegex(ValueError, "gzip"):
            digest_stream(Response(b"badbytes"), max_bytes=10000)


if __name__ == "__main__":
    unittest.main()
