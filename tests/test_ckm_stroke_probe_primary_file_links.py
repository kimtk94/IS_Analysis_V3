"""Source distribution link probes must remain source-candidate/HOLD only."""
import io
import unittest

from scripts.ckm_stroke.probe_primary_file_links import LinkParser, inspect_one, pull_html


class Fake(io.BytesIO):
    def __init__(self, txt, url):
        super().__init__(txt.encode("utf-8"))
        self.url = url

    def geturl(self):
        return self.url


class Tests(unittest.TestCase):
    def test_html_anchor_extraction(self):
        p = LinkParser()
        p.feed('<a href="/file.gz">Download <strong>Updated</strong></a>')
        self.assertEqual(p.links, [("Download Updated", "/file.gz")])

    def test_giant_links_confirm_only_distribution(self):
        page = """<a href="/corrected.gz">Download Updated Meta-analysis Locke et al + UKBiobank 2018 GZIP</a>
<a href="/README.txt">Download README Summary Statistics for Yengo et al 2018</a>"""
        def opener(req, timeout):
            return Fake(page, req.full_url)
        got = inspect_one("BMI_GIANT", opener)
        self.assertTrue(all(x["status"] == "ANCHOR_FOUND" for x in got["links_checked"]))
        self.assertFalse(got["binary_downloaded"])
        self.assertFalse(got["ready_for_ldsc"])

    def test_missing_t2d_anchor_is_not_missing_dataset(self):
        def opener(req, timeout):
            return Fake("<div>T2DGGI GWAS EUR ancestry meta-analysis summary statistics</div>", req.full_url)
        got = inspect_one("T2D_DIAGRAM", opener)
        self.assertTrue(got["requires_user_terms"])
        self.assertEqual(got["links_checked"][0]["status"], "HOLD_NO_UNIQUE_ANCHOR")

    def test_redirect_is_rejected(self):
        def opener(req, timeout):
            return Fake("<body>Test</body>", "https://not-the-official-host.com")
        with self.assertRaisesRegex(ValueError, "unapproved source redirect"):
            pull_html("https://www.diagram-consortium.org/downloads.html", opener)


if __name__ == "__main__":
    unittest.main()
