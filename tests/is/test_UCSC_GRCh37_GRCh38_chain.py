#!/usr/bin/env python3
"""Synthetic UCSC source manifest and chain mapping contract tests."""
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"is"))
from audit_UCSC_GRCh37_GRCh38_chain import (
    checksum_from_publisher_manifest,unique_same_strand_mapping,
    map_source_interval,CHAIN_IDS,
)


class MockLifter:
    def __init__(self,results=None):
        self.results=results or {}
    def convert_coordinate(self,c,p0):
        if (c,p0) in self.results:return self.results[(c,p0)]
        return [(c,p0+120000,"+",5000)]


class UcscChainQC(unittest.TestCase):
    def test_manifest_verified_filename_checksum(self):
        fixture=(
            "9e19f5942bc0352d470b2242912fe8b1  hg19ToHg38.over.chain.gz\n"
            "00000000000000000000000000000000  other.chain.gz\n"
        ).encode()
        self.assertEqual(checksum_from_publisher_manifest(
            fixture,"hg19ToHg38.over.chain.gz"),
            "9e19f5942bc0352d470b2242912fe8b1")

    def test_unknown_chain_publisher_md5_fails_closed(self):
        with self.assertRaisesRegex(ValueError,"missing"):
            checksum_from_publisher_manifest(b"fake contents\n","hg38ToHg19.over.chain.gz")

    def test_1based_mapping_does_not_lose_first_base(self):
        x=unique_same_strand_mapping(MockLifter(),"chr1",1000000)
        self.assertEqual(x["mapped_pos1"],1120000)
        self.assertEqual(x["orientation"],"+")
        self.assertEqual(x["status"],"UNIQUE_CHAIN_COORDINATE_MAPPING_NOT_ALLELE_VALIDATED")

    def test_multimapping_is_never_accepted(self):
        m=MockLifter({("chr2",999):[("chr2",2000,"+",12),("chr2",3000,"+",5)]})
        self.assertEqual(unique_same_strand_mapping(m,"chr2",1000)["status"],"UNMAPPED_OR_MULTIMAPPED")

    def test_absent_mapping_is_hold(self):
        m=MockLifter({("chr3",999):[]})
        self.assertEqual(unique_same_strand_mapping(m,"chr3",1000)["status"],"UNMAPPED_OR_MULTIMAPPED")

    def test_interval_geometry_only_not_indel_allele(self):
        x=map_source_interval(MockLifter(),"chr1",54712,54716)
        self.assertTrue(x["coordinate_geometry_consistent"])
        self.assertEqual(x["REF_ALT_base_compatibility"],"NOT_TESTED")
        self.assertFalse(x["GWAS_effect_allele_harmonized"])
        self.assertFalse(x["can_directly_run_stroke_coloc"])

    def test_indel_breakpoint_geometry_fail(self):
        m=MockLifter({("chr1",54715):[("chr1",54716+500000,"+",1000)]})
        x=map_source_interval(m,"chr1",54712,54716)
        self.assertFalse(x["coordinate_geometry_consistent"])

    def test_reverse_orientation_indel_geometry(self):
        m=MockLifter({("chr1",54711):[("chr1",120000,"-",9000)],
                      ("chr1",54715):[("chr1",119996,"-",9000)]})
        x=map_source_interval(m,"chr1",54712,54716)
        self.assertTrue(x["coordinate_geometry_consistent"])
        self.assertFalse(x["GWAS_effect_allele_harmonized"])

    def test_only_autosomal_valid_coordinates(self):
        with self.assertRaises(ValueError):
            unique_same_strand_mapping(MockLifter(),"chrX",1000)
        with self.assertRaises(ValueError):
            unique_same_strand_mapping(MockLifter(),"chr1",0)

    def test_bidirectional_chain_identities(self):
        self.assertEqual(len(CHAIN_IDS),2)
        self.assertIn("GRCh38_to_GRCh37",CHAIN_IDS)

if __name__=="__main__":
    unittest.main()
