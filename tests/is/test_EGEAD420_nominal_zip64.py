#!/usr/bin/env python3
"""Synthetic ZIP64 and manifest tests; 43GB public source fetched only in CI."""
import gzip
import io
import struct
import sys
import tarfile
import unittest
import zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"is"))
from audit_EGEAD420_nominal_zip64 import (
    EXPECTED_ZIP_SIZE,EXPECTED_PUBLISHED_ZIP_MD5,
    source_filelist,zip64_eocd,zip_directory_members,nested_gzip_tar_head,
    check_url,ZIP,FILELIST,
)


def mock_manifest():
    return ("\n".join([
        "#Archive/File\tName\tTime\tSize\tMD5",
        "File\tE-GEAD-420.idf.txt\t2021-06-09\t2858\t0672e92775fd103fd4a848d653de1f1a",
        f"Archive\tE-GEAD-420.processed.zip\t2021-06-09\t{EXPECTED_ZIP_SIZE}\t{EXPECTED_PUBLISHED_ZIP_MD5}",
        "File\tnominal_eQTL_2.tar.gz\t2021-06-09\t43495814577\t41544fb671177aeee7fd2d0634b90ae2",
    ])+"\n").encode()


class NominalZip64Tests(unittest.TestCase):
    def test_exact_NBDC_source(self):
        check_url(ZIP);check_url(FILELIST)
        with self.assertRaisesRegex(ValueError,"Unapproved|Source"):
            check_url("https://evil.example/public/E-GEAD-420.processed.zip")

    def test_publisher_manifest(self):
        d=source_filelist(mock_manifest())
        self.assertEqual(d["E-GEAD-420.processed.zip"]["size_bytes"],43507433191)
        self.assertEqual(d["nominal_eQTL_2.tar.gz"]["size_bytes"],43495814577)

    def test_changed_publisher_filelist_hash_fails(self):
        bad=mock_manifest().replace(EXPECTED_PUBLISHED_ZIP_MD5.encode(),b"0"*32)
        with self.assertRaisesRegex(ValueError,"manifest drift"):
            source_filelist(bad)

    def test_small_zip_eocd(self):
        buf=io.BytesIO()
        with zipfile.ZipFile(buf,"w") as writer:
            writer.writestr("nominal_eQTL_2.tar.gz",b"test")
        b=buf.getvalue()
        meta={"start":0,"end":len(b)-1,"total":len(b)}
        parsed=zip64_eocd(b,meta)
        self.assertEqual(parsed["members_expected"],1)
        self.assertFalse(parsed["ZIP64_footer_required"])
        start=parsed["central_directory_start"]
        end=parsed["central_directory_end_inclusive"]
        members=zip_directory_members(b[start:end+1],1)
        self.assertEqual(members[0]["name"],"nominal_eQTL_2.tar.gz")
        self.assertEqual(members[0]["uncompressed_size"],4)

    def test_zip64_eocd_locator(self):
        rec=struct.pack("<IQHHIIQQQQ",0x06064B50,44,45,45,0,0,1,1,100,700)
        loc=struct.pack("<IIQI",0x07064B50,0,0,1)
        eocd=struct.pack("<IHHHHIIH",0x06054B50,0,0,0xffff,0xffff,0xffffffff,0xffffffff,0)
        tail=rec+loc+eocd
        data=zip64_eocd(tail,{"start":0,"end":len(tail)-1,"total":1000})
        self.assertTrue(data["ZIP64_footer_required"])
        self.assertEqual(data["central_directory_start"],700)
        self.assertEqual(data["central_directory_byte_count"],100)

    def test_zip64_member_extended_sizes(self):
        name=b"nominal_eQTL_2.tar.gz"
        extra=struct.pack("<HHQQQ",1,24,43495814577,43495700000,1234567890)
        packed=struct.pack("<IHHHHHHIIIHHHHHII",
            0x02014B50,45,45,0,0,0,0,0,0xffffffff,0xffffffff,
            len(name),len(extra),0,0,0,0,0xffffffff)+name+extra
        record=zip_directory_members(packed,1)[0]
        self.assertEqual(record["uncompressed_size"],43495814577)
        self.assertEqual(record["compressed_size"],43495700000)
        self.assertEqual(record["local_header_offset"],1234567890)

    def test_gzip_nested_tar_first_member(self):
        fileobj=io.BytesIO()
        with tarfile.open(fileobj=fileobj,mode="w:gz") as archive:
            payload=b"variant\tgene\tbeta\tp\n"
            info=tarfile.TarInfo("nominal_eQTL_CD4.tsv")
            info.size=len(payload)
            archive.addfile(info,io.BytesIO(payload))
        row=nested_gzip_tar_head(fileobj.getvalue())[0]
        self.assertEqual(row["member"],"nominal_eQTL_CD4.tsv")
        self.assertTrue(row["contents_not_validated"])

    def test_bad_nested_gzip_not_validated(self):
        with self.assertRaisesRegex(ValueError,"not gzip"):
            nested_gzip_tar_head(b"invalid")

if __name__=="__main__":
    unittest.main()
