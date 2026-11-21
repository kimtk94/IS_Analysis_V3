import unittest
from pathlib import Path
from urllib.error import URLError
from scripts.ckm_stroke import audit_public_sources as audit

ROOT = Path(__file__).resolve().parents[1]


class Response:
    def __init__(self, url, data=b"GCST90104542", status=200):
        self.url = url
        self.data = data
        self.status = status
        self.headers = {"Content-Length": "12345", "Accept-Ranges": "bytes"}
    def __enter__(self):
        return self
    def __exit__(self, *args):
        return False
    def read(self, limit):
        return self.data[:limit]
    def geturl(self):
        return self.url


class TestSourceAudit(unittest.TestCase):
    def test_inventory_all_14_outcomes_and_exposures(self):
        rows = audit.read_inventory(ROOT / "config/ckm_stroke/public_sources.tsv")
        self.assertEqual(len(rows), 14)
        self.assertEqual(sum(x["kind"] == "pending" for x in rows), 6)

    def test_landing_marker_checked(self):
        x = {"kind": "page", "url": "https://www.ebi.ac.uk/gwas/studies/GCST90104542",
             "needle": "GCST90104542"}
        r = audit.check_one(x, opener=lambda req, timeout: Response(req.full_url))
        self.assertEqual(r["status"], "PAGE_VERIFIED")

    def test_missing_marker_holds(self):
        x = {"kind": "page", "url": "https://www.ebi.ac.uk/gwas/studies/GCST90104542",
             "needle": "not-on-page"}
        r = audit.check_one(x, opener=lambda req, timeout: Response(req.full_url))
        self.assertEqual(r["status"], "HOLD")

    def test_binary_head_only(self):
        x = {"kind": "file", "url": "https://ckdgen.imbi.uni-freiburg.de/test.gz", "needle": ""}
        def mock(req, timeout):
            self.assertEqual(req.get_method(), "HEAD")
            return Response(req.full_url)
        r = audit.check_one(x, opener=mock)
        self.assertEqual(r["status"], "BINARY_REACHABLE")
        self.assertNotIn("checksum_verified", r)

    def test_pending_holds(self):
        self.assertEqual(audit.check_one({"kind": "pending"})["status"], "HOLD")

    def test_external_failure_is_not_success(self):
        x = {"kind": "page", "url": "https://www.ebi.ac.uk/gwas/studies/GCST90104542",
             "needle": "GCST90104542"}
        def deny(req, timeout):
            raise URLError("blocked")
        self.assertEqual(audit.check_one(x, opener=deny)["status"], "SOURCE_UNAVAILABLE")

    def test_forbid_rogue_host(self):
        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as d:
            p = Path(d) / "sources.tsv"
            p.write_text("trait\tancestry\tsource\tkind\turl\tneedle\texpected_build\n"
                         "BMI\tEUR\tfake\tfile\thttps://evil.invalid/a.gz\t\tGRCh37\n")
            with self.assertRaises(ValueError):
                audit.read_inventory(p)

if __name__ == "__main__":
    unittest.main()
