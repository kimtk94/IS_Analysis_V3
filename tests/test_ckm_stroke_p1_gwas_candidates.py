"""P1 candidate manifest cannot be confused with approved full GWAS inputs."""
import io
import pathlib
import tempfile
import unittest

from scripts.ckm_stroke.audit_p1_gwas_candidates import audit, manifest, probe


class Response(io.BytesIO):
    def __init__(self, text, url):
        super().__init__(text.encode())
        self._url = url
        self.status = 200

    def geturl(self):
        return self._url


class Tests(unittest.TestCase):
    def sample(self, **changes):
        row = {
            "trait": "BMI", "priority": "primary", "study": "GIANT",
            "ancestry": "EUR", "landing_url": "https://giant-consortium.web.broadinstitute.org/GIANT_consortium_data_files",
            "metadata_source_url": "https://giant-consortium.web.broadinstitute.org/GIANT_consortium_data_files",
            "reported_build": "REVIEW_REQUIRED", "reported_sample_size": "",
            "published_md5": "", "page_marker": "BMI and Height",
            "access_and_qc_notes": "pending", "status": "CANDIDATE_ONLY",
        }
        return {**row, **changes}

    def test_landing_match_never_approves(self):
        row = self.sample()
        def opener(req, timeout):
            return Response("BMI and Height", req.full_url)
        result = audit([row], opener=opener)
        self.assertEqual(result["observations"][0]["status"], "LANDING_MARKER_FOUND")
        self.assertFalse(result["ready_for_ldsc"])
        self.assertFalse(result["source_approved"])

    def test_redirect_outside_allowlist_hold(self):
        def opener(req, timeout):
            return Response("BMI and Height", "https://example.org/data")
        result = probe(self.sample(), opener=opener)
        self.assertEqual(result["status"], "HOLD_REDIRECT")

    def test_missing_marker_hold(self):
        def opener(req, timeout):
            return Response("Loading...", req.full_url)
        result = probe(self.sample(), opener=opener)
        self.assertEqual(result["status"], "HOLD_MARKER_NOT_FOUND")

    def test_reject_promoted_candidate_and_duplicate(self):
        rows = [self.sample(status="APPROVED")]
        with tempfile.TemporaryDirectory() as d:
            p = pathlib.Path(d) / "c.tsv"
            import csv
            with p.open("w", encoding="utf-8", newline="") as f:
                wr = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
                wr.writeheader(); wr.writerows(rows)
            with self.assertRaises(ValueError):
                manifest(p)


if __name__ == "__main__":
    unittest.main()
