"""Synthetic KCPS2 public large-ZIP central directory and selected-member CRC gate."""
import csv,gzip,hashlib,io,json,os,struct,sys,tempfile,unittest,zipfile
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from acquire_kcps2_selective_zip_member import FILES,central,info,locate_member,fetch_trait,ROOT,TAIL,META

class KCPS2Member(unittest.TestCase):
    def make_fake(self,root):
        raw=io.BytesIO()
        payload=gzip.compress(b"CHR\tPOS\tBETA\n12\t112241766\t-.59\n",mtime=0)
        with zipfile.ZipFile(raw,"w",compression=zipfile.ZIP_DEFLATED) as z:
            for trait in ("alcohol","sbp","dbp"):
                z.writestr(FILES[trait],payload)
            for i in range(33):
                z.writestr(f"SAIGE_FAKE_{i:02d}_INFO.txt.gz",gzip.compress(b"header\n",mtime=0))
        binary=raw.getvalue()
        self.assertLess(len(binary),131072)
        (root/TAIL.name).write_bytes(bytes(131072-len(binary))+binary)
        (root/META.name).write_text(json.dumps({
            "id":15132424,"files":[{"key":"KCPS2_SAIGE_sumstats.zip",
               "size":len(binary),"checksum":"md5:0c6e3f69ca73af97417ce8c4ea8d307c"}]}))
        return binary,payload
    def test_central_directory_lists_36_and_schemas(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);blob,_=self.make_fake(root)
            x=central(root)
            self.assertEqual(len(x),36)
            self.assertEqual(x[FILES["alcohol"]]["compression_method"],8)
            self.assertEqual(info(root),len(blob))
    def test_selected_member_verified_crc_not_zip_md5(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);blob,payload=self.make_fake(root)
            ref=central(root)[FILES["alcohol"]]
            offset=ref["local_header_offset"]
            header=blob[offset:offset+512]
            with patch("acquire_kcps2_selective_zip_member.get_exact",return_value=header):
                start=locate_member(FILES["alcohol"],{FILES["alcohol"]:ref},len(blob))
                compressed=blob[start:start+ref["compressed_bytes"]]
                stage=root/(FILES["alcohol"]+".zipdeflate.ranges")
                stage.mkdir()
                (stage/f"{start:012d}-{start+len(compressed)-1:012d}.bin").write_bytes(compressed)
                result=fetch_trait(root,"alcohol",workers=1,chunk_mib=1)
            self.assertTrue(result["complete_member_retrieved"])
            self.assertEqual(result["selected_member_only"],True)
            self.assertFalse(result["original_archive_md5_verified"])
            self.assertEqual((root/FILES["alcohol"]).read_bytes(),payload)
            reused=fetch_trait(root,"alcohol",workers=1,chunk_mib=1)
            self.assertEqual(reused["inner_gzip_sha256"],hashlib.sha256(payload).hexdigest())
    def test_reject_changed_zip_central_directory(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);self.make_fake(root)
            p=root/TAIL.name;p.write_bytes(b"0"*131072)
            with self.assertRaisesRegex(ValueError,"Missing ZIP"):
                central(root)
if __name__=="__main__":unittest.main()
