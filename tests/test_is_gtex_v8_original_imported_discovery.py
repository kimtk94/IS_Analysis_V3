import gzip
import importlib.util
from pathlib import Path
import unittest

PATH=Path(__file__).resolve().parents[1]/"scripts/is/audit_gtex_v8_original_imported_discovery.py"
spec=importlib.util.spec_from_file_location("gtex_v8_discover",PATH)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

QTD_HTML=b'<html><a href="QTD000161.all.tsv.gz">all</a><a href="QTD000161.all.tsv.gz.tbi">index</a></html>'
IMP_HTML=b'<html><a href="GTEx_V8/">GTEx v8</a></html>'
CHILD_HTML=b'<html><a href="Brain_Cerebellar_Hemisphere.all.tsv.gz">original</a></html>'
HEADER=b"molecular_trait_id\tvariant\tpvalue\tbeta\nENSG0\tchr1_42_A_C\t1e-2\t0.3\n"

class GTExImportedDiscoveryTests(unittest.TestCase):
    def mock(self,url):
        value={
            m.QTL_DIR:QTD_HTML,m.IMPORTED:IMP_HTML,
            m.IMPORTED+"GTEx_V8/":CHILD_HTML,
            m.QTL_URL:gzip.compress(HEADER,mtime=0),
        }
        return value[url]

    def test_success_discovery_does_not_approve_source_identity(self):
        out=m.audit(self.mock)
        self.assertEqual(out["source_discovery_gate"],
            "PUBLIC_SOURCE_INDEXES_DISCOVERED_ORIGINAL_QTL_COMPARISON_PENDING")
        self.assertEqual(out["GTEx_v8_imported_child_index"]["gate"],
                         "OFFICIAL_IMPORTED_CHILD_INDEX_FOUND")
        self.assertFalse(out["original_GTEx_v8_FastQTL_effects_compared"])
        self.assertEqual(out["scientific_status"],"HOLD")

    def test_qtd_header_bounded_gzip(self):
        out=m.audit(self.mock)
        self.assertIn("molecular_trait_id",out["reprocessed_QTD_gzip_source_prefix"]["header"])
        self.assertEqual(out["reprocessed_QTD_gzip_source_prefix"]["gate"],
                         "SOURCE_QTD_HEADER_FOUND_NOT_RELEASE_PIN")

    def test_missing_qtd_index_remains_hold(self):
        def bad(url):
            if url==m.QTL_DIR:raise OSError("FTP unreachable")
            return self.mock(url)
        out=m.audit(bad)
        self.assertEqual(out["source_discovery_gate"],"HOLD_SOURCE_INDEX_INCOMPLETE")

    def test_empty_imported_index_fails_closed(self):
        def empty(url):
            return b"<html>no listing</html>" if url==m.IMPORTED else self.mock(url)
        out=m.audit(empty)
        self.assertEqual(out["source_discovery_gate"],"HOLD_SOURCE_INDEX_INCOMPLETE")
        self.assertFalse(out["GTEx_v8_original_imported_retrieved"])

    def test_non_gzip_qtd_header_does_not_pass(self):
        def plain(url):
            return b"not a gzip source" if url==m.QTL_URL else self.mock(url)
        out=m.audit(plain)
        self.assertEqual(out["reprocessed_QTD_gzip_source_prefix"]["gate"],
                         "HOLD_UNREADABLE_GZIP_HEADER")

    def test_public_source_url_allowlist(self):
        for url in ["https://evil.example/QTD000161",
                    "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/../../private"]:
            with self.assertRaises(ValueError):
                m.verified_url(url)

    def test_bounded_documentation_index(self):
        with self.assertRaises(ValueError):
            m.fetch(m.QTL_DIR,4,fetcher=lambda _:b"123456")

    def test_imported_multiple_folders_does_not_guess_source(self):
        def two(url):
            if url==m.IMPORTED:
                return b'<a href="GTEx_V8/">GTEx v8</a><a href="Other/">Other</a>'
            return self.mock(url)
        out=m.audit(two)
        self.assertNotIn("GTEx_v8_imported_child_index",out)
        self.assertFalse(out["original_GTEx_v8_FastQTL_effects_compared"])

if __name__=="__main__":
    unittest.main()
