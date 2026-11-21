import gzip
import hashlib
import io
import tempfile
import unittest
from pathlib import Path
from scripts.ckm_stroke import ingest_uacr as pilot

class Response:
    status=200
    def __init__(self,payload,length=None):
        self.stream=io.BytesIO(payload)
        self.headers={"Content-Length":str(length if length is not None else len(payload))}
    def __enter__(self): return self
    def __exit__(self,*args): return False
    def read(self,size): return self.stream.read(size)

class Tests(unittest.TestCase):
    def source(self):
        return gzip.compress(b"Chr Pos_b37 RSID Allele1 Allele2 Freq1 Effect StdErr P-value n_total_sum\n"
                             b"1 736289 rs79010578 a t 0.14 0.000315 0.0007081 0.6563 347599\n")
    def test_real_algorithm_on_synthetic_gz_with_expected_digest(self):
        data=self.source()
        actual=hashlib.md5(data).hexdigest()
        prior=pilot.UACR["md5"]
        try:
            pilot.UACR["md5"]=actual
            with tempfile.TemporaryDirectory() as d:
                x=pilot.ingest(Path(d),max_mb=1,opener=lambda req,timeout: Response(data))
                self.assertEqual(x["status"],"PASS")
                self.assertTrue(x["source_md5_verified"])
                self.assertEqual(x["normalized_rows"],1)
        finally:
            pilot.UACR["md5"]=prior
    def test_mismatched_checksum_fails(self):
        data=self.source()
        with tempfile.TemporaryDirectory() as d:
            x=pilot.ingest(Path(d),max_mb=1,opener=lambda req,timeout: Response(data))
            self.assertEqual(x["status"],"FAIL")
            self.assertEqual(x["reason"],"original MD5 mismatch")
    def test_hard_cap_is_enforced(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):
                pilot.stream_to_file(Response(b"123456"),Path(d)/"file.gz",3)

if __name__=="__main__": unittest.main()
