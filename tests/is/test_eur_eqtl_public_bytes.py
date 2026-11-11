#!/usr/bin/env python3
"""Offline and synthetic-only safety tests for public eQTL source byte prefix."""
import gzip
import io
import sys
import tarfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"is"))
from audit_eur_eqtl_public_bytes import (
    GTEX_EUR_URL,EBI_BLUEPRINT_URL,REQUIRED_ALL_STATS,SOURCE_BUILD,
    first_tar_member,gzip_prefix_header,safe_url,
)


class EqtlsActualSourcePrepTests(unittest.TestCase):
    def test_correct_source_object_names(self):
        self.assertTrue(GTEX_EUR_URL.endswith("GTEx_Analysis_v8_eQTL_EUR.tar"))
        self.assertTrue(EBI_BLUEPRINT_URL.endswith("QTD000021.all.tsv.gz"))
        self.assertIn("GRCh37",SOURCE_BUILD["GIGASTROKE_EUR_EAS"])
        self.assertIn("GRCh38",SOURCE_BUILD["eQTL_Catalogue_BLUEPRINT"])

    def test_safe_approved_https_sources(self):
        safe_url(GTEX_EUR_URL)
        safe_url(EBI_BLUEPRINT_URL)
        with self.assertRaises(ValueError):
            safe_url("https://not-ebi.example/file")
        with self.assertRaises(ValueError):
            safe_url("http://ftp.ebi.ac.uk/file")

    def test_realistic_tar_first_header(self):
        t=io.BytesIO()
        with tarfile.open(fileobj=t,mode="w") as tar:
            d=b"Header\nFirst line\n"
            member=tarfile.TarInfo("GTEx_Analysis_v8_eQTL_EUR/README.txt")
            member.size=len(d)
            tar.addfile(member,io.BytesIO(d))
        data=t.getvalue()
        result=first_tar_member(data[:262144])
        self.assertTrue(result["ustar_magic_present"])
        self.assertEqual(result["first_member_size_bytes"],len(d))
        self.assertIn("README",result["first_member_name"])
        self.assertTrue(result["does_not_enumerate_entire_tar"])

    def test_not_tar_is_rejected(self):
        with self.assertRaises(ValueError):
            first_tar_member(b"not tar"*80)

    def test_gzip_nominal_association_schema(self):
        cols=[
          "variant","r2","pvalue","molecular_trait_id","maf","gene_id",
          "beta","se","chromosome","position","ref","alt","an","ac"
        ]
        txt="\t".join(cols)+"\nchr1_101_A_G\t0.99\t0.1\tENSG1\t0.22\tENSG1\t0.1\t0.02\t1\t101\tA\tG\t200\t88\n"
        result=gzip_prefix_header(gzip.compress(txt.encode("utf-8")))
        self.assertTrue(result["full_nominal_required_fields_present"])
        self.assertEqual(result["missing_required_fields"],[])
        self.assertTrue(result["still_need_association_rows_and_variant_coverage"])

    def test_significant_only_schema_does_not_qualify(self):
        txt="variant\tgene_id\tpvalue\nchr1_101_A_G\tENSG1\t0.01\n"
        result=gzip_prefix_header(gzip.compress(txt.encode()))
        self.assertFalse(result["full_nominal_required_fields_present"])
        self.assertIn("se",result["missing_required_fields"])

    def test_duplicate_headers_are_not_acceptable(self):
        with self.assertRaisesRegex(ValueError,"Duplicate"):
            gzip_prefix_header(gzip.compress(b"variant\tbeta\tbeta\tse\tpvalue\n"))

    def test_truncated_gzip_is_rejected(self):
        with self.assertRaises((ValueError,EOFError,OSError)):
            gzip_prefix_header(b"not gzip content")

if __name__=="__main__":
    unittest.main()
