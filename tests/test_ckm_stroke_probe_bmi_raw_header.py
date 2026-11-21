"""Tests that a bounded BMI header probe cannot approve a source."""
import gzip
import io
import unittest

from scripts.ckm_stroke.probe_bmi_raw_header import audit, MAX_PREFIX, SOURCE, README


class Fake(io.BytesIO):
    def __init__(self, body, url, status=200):
        super().__init__(body)
        self._url = url
        self.status = status
        self.headers = {"Content-Length": str(len(body))}

    def geturl(self):
        return self._url


def good_opener(req, timeout):
    from scripts.ckm_stroke.probe_bmi_raw_header import EXPECTED
    if req.full_url == README:
        text = "Genome build hg19\nHealth and Retirement Study\n"
        text += "\n".join(col + ": field" for col in EXPECTED)
        return Fake(text.encode(), README)
    assert req.full_url == SOURCE
    data = gzip.compress(("\t".join(EXPECTED) + "\n" + "1\t2\t3\n").encode())
    return Fake(data, SOURCE, 206)


class TestBmiHeader(unittest.TestCase):
    def test_reads_valid_prefix_without_promoting(self):
        result = audit(good_opener)
        self.assertEqual(result["technical_header_status"], "PASS_HEADER_ONLY")
        self.assertFalse(result["approved_for_ldsc"])
        self.assertFalse(result["full_gwas_downloaded"])
        self.assertFalse(result["raw_md5_verified"])
        self.assertEqual(result["scientific_status"], "HOLD")

    def test_redirect_produces_hold(self):
        def opener(req, timeout):
            return Fake(b"bad", "https://example.net/data", 200)
        self.assertEqual(audit(opener)["technical_header_status"], "SOURCE_UNAVAILABLE_OR_INVALID")

    def test_wrong_header_produces_hold(self):
        def opener(req, timeout):
            if req.full_url == README:
                return good_opener(req, timeout)
            return Fake(gzip.compress(b"A\tB\n1\t2\n"), SOURCE, 200)
        self.assertEqual(audit(opener)["technical_header_status"], "HOLD_SCHEMA_DISAGREEMENT")


if __name__ == "__main__":
    unittest.main()
