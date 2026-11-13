#!/usr/bin/env python3
"""Zip range source contract tests; no individual Japanese donor data."""
import io
from pathlib import Path
import struct
import sys
import unittest
import zipfile

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"is"))
from audit_japanese_eqtl_zip_central_dir import (
    SOURCE_URL, parse_eocd, list_zip_members, source_url_ok,
)


class JapaneseEqtlsZipManifestTests(unittest.TestCase):
    def test_archive_source_identity(self):
        self.assertIn("hum0099.v1.eqtl.v1.zip",SOURCE_URL)
        source_url_ok(SOURCE_URL)
        with self.assertRaises(ValueError):
            source_url_ok("https://attacker.com/mydata.zip")

    def make_sample_zip(self):
        b=io.BytesIO()
        with zipfile.ZipFile(b,"w",compression=zipfile.ZIP_DEFLATED) as archive:
            archive.writestr("hum0099/README.txt","Gene cis QTL associations\n")
            archive.writestr("hum0099/test_assoc.tsv","chr\tpos\tgene\tbeta\tse\n")
        return b.getvalue()

    def test_small_real_zip_syntax_can_be_audited_from_tail(self):
        raw=self.make_sample_zip()
        length=min(2048,len(raw))
        tail=raw[-length:]
        start=len(raw)-length
        d=parse_eocd(tail,start,len(raw))
        self.assertEqual(d["central_directory_entries"],2)
        self.assertTrue(d["no_zip_members_downloaded"])
        p=d["central_directory_offset"]
        cd=raw[p:p+d["central_directory_bytes"]]
        files=list_zip_members(cd,d["central_directory_entries"])
        self.assertEqual(len(files["files"]),2)
        self.assertIn("README",files["files"][0]["member"])
        self.assertTrue(all(x["member_data_not_fetched"] for x in files["files"]))

    def test_corrupt_zip_fails(self):
        with self.assertRaisesRegex(ValueError,"EOCD"):
            parse_eocd(b"x"*100,0,100)

    def test_central_directory_truncation_fails(self):
        raw=self.make_sample_zip()
        end=raw.rfind(b"PK\x05\x06")
        d=parse_eocd(raw[:],0,len(raw))
        start=d["central_directory_offset"]
        segment=raw[start:start+d["central_directory_bytes"]]
        with self.assertRaisesRegex(ValueError,"Truncated|Unexpected|Malformed"):
            list_zip_members(segment[:-5],d["central_directory_entries"])

    def test_restricted_liftover_not_claimed_by_members(self):
        raw=self.make_sample_zip()
        d=parse_eocd(raw,0,len(raw))
        segment=raw[d["central_directory_offset"]:d["central_directory_end"]+1]
        ret=list_zip_members(segment,2)
        self.assertTrue(ret["all_zip_member_names_covered"])
        self.assertEqual(ret["files"][1]["uncompressed_bytes"],len("chr\tpos\tgene\tbeta\tse\n"))
        # File headers alone never imply full cis-eQTL completeness.

if __name__=="__main__":
    unittest.main()
