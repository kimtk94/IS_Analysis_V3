import gzip
import io
import unittest
from scripts.ckm_stroke import probe_headers as p

class Response:
    status=206
    headers={"Content-Range":"bytes 0-100/5000"}
    def __init__(self,data,url): self.data,self.url=data,url
    def __enter__(self): return self
    def __exit__(self,*args): pass
    def geturl(self): return self.url
    def read(self,n): return self.data[:n]

class Tests(unittest.TestCase):
    def test_plain_header(self):
        x=p.extract_first_lines(b"CHR\tPOS\n1\t20\n")
        self.assertEqual(x["encoding"],"plain")
        self.assertTrue(x["complete_first_line"])
    def test_gzip_partial(self):
        data=gzip.compress(b"Chr\tPos_b37\tAllele1\tAllele2\n1\t33\tA\tT\n")
        x=p.extract_first_lines(data[:len(data)-7])
        self.assertEqual(x["first_lines"][0],"Chr\tPos_b37\tAllele1\tAllele2")
    def test_range_request(self):
        url="https://ckdgen.imbi.uni-freiburg.de/demo.gz"
        def opener(req,timeout):
            self.assertIn("bytes=0-",req.headers["Range"])
            return Response(gzip.compress(b"Chr\tPos_b37\n1\t30\n"),url)
        x=p.inspect(url,opener=opener)
        self.assertEqual(x["status"],"HEADER_OBSERVED")
    def test_host_rejected(self):
        with self.assertRaises(ValueError):
            p.inspect("https://example.org/bad.gz",opener=lambda req,timeout: None)

if __name__=="__main__": unittest.main()
