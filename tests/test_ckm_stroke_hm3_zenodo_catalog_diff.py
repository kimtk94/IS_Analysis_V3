"""Offline fixtures for pinned HapMap3 MRCIEU-vs-Zenodo+LD score comparison."""
import gzip
import hashlib
import io
import json
import tarfile
import unittest
from unittest.mock import patch

from scripts.ckm_stroke.audit_bmi_hm3_overlap import allele_pair
from scripts.ckm_stroke.audit_hm3_zenodo_catalog_diff import (
    read_bundle_hm3, parse_archive_response, parse_ldscore_member, summarize,
    audit, BUNDLED_NAME, FILE_URL
)


class Response(io.BytesIO):
    def __init__(self, content, length=None):
        super().__init__(content)
        self.status=200
        self.headers={"Content-Length": str(len(content) if length is None else length)}
    def geturl(self):return FILE_URL


def fake_bundle():
    return (b"SNP\tA1\tA2\n"
            b"rs1\tA\tC\n"
            b"rs2\tT\tC\n"
            b"rs3\tA\tT\n"
            b"rs4\tC\tG\n")


def bundle_fixture(missing_m=None):
    output=io.BytesIO()
    member_data=fake_bundle()
    with tarfile.open(fileobj=output,mode="w:gz") as tar:
        def add(name, body):
            info=tarfile.TarInfo(name)
            info.size=len(body)
            tar.addfile(info,io.BytesIO(body))
        add(BUNDLED_NAME,member_data)
        for chr in range(1,23):
            body=("CHR\tSNP\tBP\tL2\n"+f"{chr}\trs1\t100\t1.1\n"+
                  f"{chr}\trs2\t101\t1.2\n").encode()
            add(f"eur_w_ld_chr/{chr}.l2.ldscore.gz",gzip.compress(body))
            if chr != missing_m:
                add(f"eur_w_ld_chr/{chr}.l2.M_5_50",b"100\n")
    return output.getvalue(),member_data


class CatalogTests(unittest.TestCase):
    def test_bundled_header_and_alleles(self):
        data=fake_bundle()
        catalog,meta=read_bundle_hm3(io.BytesIO(data),len(data))
        self.assertEqual(len(catalog),4)
        self.assertEqual(meta["sha256"],hashlib.sha256(data).hexdigest())
        self.assertEqual(catalog["rs1"],allele_pair("A","C"))

    def test_truncated_bundle_fails_closed(self):
        with self.assertRaisesRegex(ValueError,"truncated"):
            read_bundle_hm3(io.BytesIO(fake_bundle()),len(fake_bundle())+1)

    def test_duplicate_label_allele_conflict_excluded(self):
        data=fake_bundle()+b"rs1\tG\tT\n"
        catalog,meta=read_bundle_hm3(io.BytesIO(data),len(data))
        self.assertNotIn("rs1",catalog)
        self.assertEqual(meta["conflicting_distinct_rsids_excluded"],1)

    def test_ldscore_chr_and_member_read(self):
        lines=("CHR SNP BP L2\n"+"1 rs1 100 2\n"+"1 rs3 101 2.1\n").encode()
        ids,meta=parse_ldscore_member(io.BytesIO(gzip.compress(lines)),1,min_rows=1)
        self.assertEqual(ids,{"rs1","rs3"})
        self.assertEqual(meta["rows"],2)

    def test_archive_22_chrs_and_hash(self):
        raw,bundle=bundle_fixture()
        digest=hashlib.md5(raw).hexdigest()
        with patch("scripts.ckm_stroke.audit_hm3_zenodo_catalog_diff.BUNDLED_SIZE",len(bundle)):
            catalog,hm,ld,meta=parse_archive_response(
                Response(raw),expected_size=len(raw),expected_md5=digest,min_ld_rows=1)
        self.assertEqual(len(catalog),4)
        self.assertEqual(len(meta["canonical_ld_chromosomes"]),22)
        self.assertEqual(len(meta["canonical_M_5_50_chromosomes"]),22)
        self.assertEqual(ld,{"rs1","rs2"})
        self.assertEqual(meta["archive_md5"],digest)

    def test_archive_missing_regression_weight_fails(self):
        raw,bundle=bundle_fixture(missing_m=10)
        with patch("scripts.ckm_stroke.audit_hm3_zenodo_catalog_diff.BUNDLED_SIZE",len(bundle)):
            with self.assertRaisesRegex(ValueError,"incomplete"):
                parse_archive_response(Response(raw),expected_size=len(raw),
                                       expected_md5=hashlib.md5(raw).hexdigest(),min_ld_rows=1)

    def test_mirror_hash_mismatch_is_error(self):
        raw,bundle=bundle_fixture()
        with patch("scripts.ckm_stroke.audit_hm3_zenodo_catalog_diff.BUNDLED_SIZE",len(bundle)):
            with self.assertRaisesRegex(ValueError,"MD5"):
                parse_archive_response(Response(raw),expected_size=len(raw),
                                       expected_md5="0"*32,min_ld_rows=1)

    def test_comparison_differences_not_interchangeable(self):
        mrc={"rs1":allele_pair("A","C"),"rs2":allele_pair("A","G"),
             "rs3":allele_pair("A","T"),"mrc_only":allele_pair("A","C")}
        zen={"rs1":allele_pair("A","C"),"rs2":allele_pair("T","C"),
             "rs3":allele_pair("A","T"),"zen_only":allele_pair("A","G")}
        report=summarize(mrc,zen,{"rs1","rs2","zen_only"}, {}, {}, {})
        c=report["comparison"]
        self.assertEqual(c["shared_rsids"],3)
        self.assertEqual(c["mrcieu_only"],1)
        self.assertEqual(c["zenodo_only"],1)
        self.assertEqual(c["common_rsid_allele_categories"]["pair_equal_unordered"],2)
        self.assertEqual(c["common_rsid_allele_categories"]["strand_complement_possible"],1)
        self.assertEqual(c["mrcieu_zenodo_ldscore_triple_intersection"],2)
        self.assertEqual(c["triple_direct_unordered_pair_agreement"],1)
        self.assertFalse(c["strict_catalog_interchangeable"])
        self.assertFalse(report["ready_for_ldsc"])

    def test_pinned_hm3_sha_mismatch_keeps_hold(self):
        def fake_hm3(opener):
            return {"rs1":allele_pair("A","C")},{"ref_sha256":"0"*64}
        with patch("scripts.ckm_stroke.audit_hm3_zenodo_catalog_diff.hm3_catalog",fake_hm3):
            r=audit()
        self.assertEqual(r["technical_status"],"HOLD_REFERENCE_OR_SCHEMA")
        self.assertFalse(r["ready_for_ldsc"])


if __name__=="__main__":
    unittest.main()
