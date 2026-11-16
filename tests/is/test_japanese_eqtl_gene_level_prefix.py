#!/usr/bin/env python3
"""Synthetic ZIP/nested TAR/method DOCX tests, no human research subject data."""
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
from audit_japanese_eqtl_gene_level_prefix import (
    zip_directory_entries,first_tar_records_from_nested_gzip,methods_docx_text,
    MAX_PREFIX,
)
from audit_japanese_eqtl_zip_central_dir import parse_eocd


def nested_tar_fixture():
    stream=io.BytesIO()
    with tarfile.open(fileobj=stream,mode="w:gz") as tf:
        content=b"gene\tvariant\tbeta\tse\nABO\t1_123_A_G\t0.2\t0.1\n"
        info=tarfile.TarInfo("gene_level/cell_type/nominal.tsv")
        info.size=len(content)
        tf.addfile(info,io.BytesIO(content))
    return stream.getvalue()


def minimal_docx():
    out=io.BytesIO()
    with zipfile.ZipFile(out,mode="w",compression=zipfile.ZIP_DEFLATED) as doc:
        doc.writestr("word/document.xml",
            '<?xml version="1.0"?><w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
            '<w:body><w:p><w:r><w:t>QTL methods from original study.</w:t></w:r></w:p></w:body></w:document>')
    return out.getvalue()


def zipped_source():
    out=io.BytesIO()
    with zipfile.ZipFile(out,mode="w",compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr("hum0099.v1.eqtl.v1/eQTL_gene_level.tar.gz",nested_tar_fixture())
        z.writestr("hum0099.v1.eqtl.v1/eQTL_method_riken_ishigaki.docx",minimal_docx())
    return out.getvalue()


class EASJapaneseGeneLevelTests(unittest.TestCase):
    def test_nested_tar_header_in_gzip(self):
        parsed=first_tar_records_from_nested_gzip(nested_tar_fixture())
        self.assertIn("nominal.tsv",parsed["nested_tar_first_member_headers"][0]["name"])
        self.assertFalse(parsed["full_variant_cis_pairs_examined"])

    def test_methods_wordprocessingml(self):
        result=methods_docx_text(minimal_docx())
        self.assertEqual(result["paragraphs"],1)
        self.assertTrue(result["source_method_docx_crc_verified"])
        self.assertIn("QTL methods",result["method_text_first_12000_chars"])

    def test_source_directory_member_offsets(self):
        source=zipped_source()
        eocd=parse_eocd(source,0,len(source))
        d=source[eocd["central_directory_offset"]:eocd["central_directory_end"]+1]
        records=zip_directory_entries(d,eocd["central_directory_entries"])
        self.assertEqual(len(records),2)
        self.assertEqual(records[0]["compression_method"],8)
        self.assertIn("gene_level",records[0]["name"])
        for x in records:
            self.assertGreater(x["compressed_size"],0)
            self.assertGreaterEqual(x["local_header_offset"],0)

    def test_empty_zip_methods_not_valid(self):
        buf=io.BytesIO()
        with zipfile.ZipFile(buf,"w") as f:
            f.writestr("README.txt","hello")
        with self.assertRaises(ValueError):
            methods_docx_text(buf.getvalue())

    def test_invalid_nested_tar_not_association(self):
        with self.assertRaisesRegex(ValueError,"not gzip"):
            first_tar_records_from_nested_gzip(b"not gzip "+b"x"*3000)

    def test_corrupt_undecodable_zip_directory(self):
        with self.assertRaises(ValueError):
            zip_directory_entries(b"not-a-zip"*100,2)


if __name__=="__main__":
    unittest.main()
