import unittest
from scripts.ckm_stroke import glgc_pysam_probe as pilot
from unittest.mock import patch

HEADER="rsID\tCHROM\tPOS_b37\tREF\tALT\tN\tEFFECT_SIZE\tSE\tpvalue"
ROW="rs1\t22\t20000010\tA\tG\t100000\t0.1\t0.03\t0.0008581"

class FakeTabix:
    def __init__(self,url):
        self.url=url
        self.header=[HEADER]
    def __enter__(self): return self
    def __exit__(self,*args): return False
    def fetch(self,chrom,start,end):
        assert (chrom,start,end)==("22",19999999,20050000)
        return iter([ROW])

class NoHeader(FakeTabix):
    header=[]
    def __init__(self,url):
        super().__init__(url)
        self.header=[]

class TestPysamRegion(unittest.TestCase):
    def test_complete_locus_window(self):
        result=pilot.probe("LDL_C",factory=FakeTabix)
        self.assertEqual(result["status"],"PASS")
        self.assertEqual(result["rows_region"],1)
        self.assertEqual(result["scientific_status"],"NOT_RUN")

    def test_source_header_absent_is_hold(self):
        result=pilot.probe("HDL_C",factory=NoHeader,
                           header_fetcher=lambda url: {"status":"HOLD","reason":"no header"})
        self.assertEqual(result["status"],"HOLD")

    def test_source_prefix_header_fallback_validates_real_rows(self):
        result=pilot.probe(
            "HDL_C",factory=NoHeader,
            header_fetcher=lambda url: {"status":"PASS","header":HEADER,
                                        "origin":"official_gzip_prefix","header_bytes_read":4096})
        self.assertEqual(result["status"],"PASS")
        self.assertEqual(result["rows_region"],1)
        self.assertEqual(result["header_evidence"]["origin"],"official_gzip_prefix")
        self.assertEqual(result["scientific_status"],"NOT_RUN")

    def test_official_header_prefix_rejects_wrong_schema(self):
        with patch.object(pilot,"inspect",return_value={
                "status":"HEADER_OBSERVED","first_lines":["WRONG\tPOS\tALLELE"],
                "bytes_read":1024}):
            r=pilot.fetch_official_header("https://csg.sph.umich.edu/x.gz")
        self.assertEqual(r["status"],"HOLD")

    def test_official_header_prefix_accepts_actual_columns(self):
        with patch.object(pilot,"inspect",return_value={
                "status":"HEADER_OBSERVED","first_lines":["#"+HEADER],
                "bytes_read":4096,"http_status":206}):
            r=pilot.fetch_official_header("https://csg.sph.umich.edu/x.gz")
        self.assertEqual(r["status"],"PASS")
        self.assertEqual(r["header"],HEADER)

    def test_factory_failure_is_not_success(self):
        def bad(url):
            raise OSError("index not reachable")
        result=pilot.probe("TG",factory=bad)
        self.assertEqual(result["status"],"SOURCE_UNAVAILABLE")

    def test_invalid_trait(self):
        with self.assertRaises(ValueError):
            pilot.probe("fake",factory=FakeTabix)

if __name__=="__main__":
    unittest.main()
