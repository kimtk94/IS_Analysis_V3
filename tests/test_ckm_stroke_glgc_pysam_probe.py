import unittest
from scripts.ckm_stroke import glgc_pysam_probe as pilot

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
        result=pilot.probe("HDL_C",factory=NoHeader)
        self.assertEqual(result["status"],"HOLD")

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
