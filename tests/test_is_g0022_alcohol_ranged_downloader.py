"""Synthetic exact-206 ranged transport, integrity, and fail-closed promotion."""
import hashlib,json,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from download_g0022_japanese_alcohol_ranges import fetch_range
from download_g0022_japanese_alcohol_ranges_parallel import acquire_parallel

class FakeResponse:
    status_code=206
    def __init__(self,lo,hi,total,payload):
        self.headers={"Content-Range":f"bytes {lo}-{hi}/{total}",
                      "Content-Length":str(len(payload))}
        self.payload=payload
    def __enter__(self):return self
    def __exit__(self,*_):return False
    def iter_content(self,chunk_size):
        for pos in range(0,len(self.payload),chunk_size):
            yield self.payload[pos:pos+chunk_size]
class FakeSession:
    def __init__(self,response):self.response=response
    def get(self,*args,**kwargs):return self.response

class RangeIntegrityTests(unittest.TestCase):
    def test_exact_range_success_and_no_200_fallback(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/"x.bin"
            s=FakeSession(FakeResponse(4,7,12,b"abcd"))
            self.assertEqual(fetch_range(s,"https://zenodo.org/",4,7,12,out,1),"NEW")
            self.assertEqual(out.read_bytes(),b"abcd")
            self.assertEqual(fetch_range(s,"https://zenodo.org/",4,7,12,out,1),"ALREADY_STAGED")
            invalid=FakeResponse(0,3,12,b"abcd")
            with self.assertRaisesRegex(RuntimeError,"Range response mismatch"):
                fetch_range(FakeSession(invalid),"u",4,7,12,Path(tmp)/"no",1)
    def test_truncated_and_http_200_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"x"
            fake=FakeResponse(0,3,4,b"ab")
            fake.headers["Content-Length"]="4"
            with self.assertRaisesRegex(RuntimeError,"Range truncated"):
                fetch_range(FakeSession(fake),"u",0,3,4,path,1)
            self.assertFalse(path.is_file())
            fake.status_code=200
            with self.assertRaisesRegex(RuntimeError,"requires exact 206"):
                fetch_range(FakeSession(fake),"u",0,3,4,path,1)
    def test_full_file_reassembled_only_after_original_md5(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            payload=b"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            name="1_Alcohol_intake_Unstratified.tsv.gz"
            record={"id":10038152,"files":[{"key":name,"size":len(payload),
                 "checksum":"md5:"+hashlib.md5(payload).hexdigest()}]}
            (root/"ZENODO_10038152_METADATA.json").write_text(json.dumps(record))
            def fake(job):
                _name,_url,start,end,size,path=job
                path.write_bytes(payload[start:end+1])
                return start,end,"SYNTHETIC"
            with patch("download_g0022_japanese_alcohol_ranges_parallel.fetch_one",side_effect=fake):
                x=acquire_parallel("alcohol",root,workers=2,block_size=8)
            self.assertEqual(x["status"],"OFFICIAL_COMPLETE_MD5_VERIFIED")
            self.assertEqual((root/name).read_bytes(),payload)
    def test_bad_md5_never_promotes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);name="1_Alcohol_intake_Unstratified.tsv.gz"
            (root/"ZENODO_10038152_METADATA.json").write_text(json.dumps({
                "id":10038152,"files":[{"key":name,"size":11,"checksum":"md5:"+"0"*32}]}))
            def fake(job):
                _,_,start,end,total,path=job
                path.write_bytes(b"x"*(end-start+1))
                return start,end,"SYNTHETIC"
            with patch("download_g0022_japanese_alcohol_ranges_parallel.fetch_one",side_effect=fake):
                with self.assertRaisesRegex(ValueError,"MD5 mismatch"):
                    acquire_parallel("alcohol",root,workers=1,block_size=6)
            self.assertFalse((root/name).exists())
if __name__=="__main__":unittest.main()
