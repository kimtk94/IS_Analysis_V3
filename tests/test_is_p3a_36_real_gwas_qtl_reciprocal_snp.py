"""Fail-closed 36 assay SNP build conversion and signed GIGASTROKE allele tests."""
import importlib
import json
import math
import unittest
from pathlib import Path
from unittest.mock import patch

M=importlib.import_module("scripts.is.audit_p3a_36_real_gwas_qtl_reciprocal_snp")
MANIFEST=Path("docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json")


class FakeSeq:
    def __init__(self,expected_pos,allele="A"):
        self.expected_pos=expected_pos
        self.allele=allele
    def __len__(self):return 150000000
    def __getitem__(self,index):
        if index==self.expected_pos-1:return ord(self.allele)
        return ord("A")


class FakeLift:
    def __init__(self,shift):
        self.shift=shift
    def convert_coordinate(self,chrom,zero):
        return [("chr10",zero+self.shift,"+",1000)]


def raw_fields():
    d={k:"." for k in M.SOURCE_HEADER}
    d.update({
      "chromosome":"10","position":"103000000",
      "ref":"A","alt":"C","variant":"chr10_103000000_A_C",
      "molecular_trait_id":"ENSG00000065613","gene_id":"ENSG00000065613",
      "beta":"0.23","se":"0.05","maf":"0.20","pvalue":"0.00001",
      "type":"SNP","an":"300","ac":"65",
      "median_tpm":"8","ma_samples":"25","rsid":"rs123",
    })
    return d


class SNPAlignmentTests(unittest.TestCase):
    def test_grch38_to_grch37_to_grch38_plus_strand(self):
        pos,why=M.reciprocal_allele_check(103000000,"A",
               FakeSeq(104759757),FakeSeq(103000000),
               FakeLift(1759757),FakeLift(-1759757))
        self.assertEqual(pos,104759757)
        self.assertIsNone(why)

    def test_wrong_new_reference_rejected(self):
        pos,why=M.reciprocal_allele_check(103000000,"C",
               FakeSeq(104759757),FakeSeq(103000000),
               FakeLift(1759757),FakeLift(-1759757))
        self.assertIsNone(pos)
        self.assertEqual(why,"HOLD_GRCh38_REFERENCE_MISMATCH")

    def test_wrong_old_reference_rejected(self):
        pos,why=M.reciprocal_allele_check(103000000,"A",
               FakeSeq(104759757,"C"),FakeSeq(103000000),
               FakeLift(1759757),FakeLift(-1759757))
        self.assertEqual(why,"HOLD_GRCh37_REFERENCE_MISMATCH")

    def test_nonreciprocal_rejected(self):
        pos,why=M.reciprocal_allele_check(103000000,"A",
               FakeSeq(104759757),FakeSeq(103000000),
               FakeLift(1759757),FakeLift(-1759756))
        self.assertEqual(why,"HOLD_RECIPROCAL_LIFTOVER_DISCORDANCE")

    def test_duplicate_liftover_mapping_rejected(self):
        class Duplicated(FakeLift):
            def convert_coordinate(self,c,p):
                return [("chr10",p+self.shift,"+",100),
                        ("chr10",p+self.shift+1,"+",20)]
        pos,why=M.reciprocal_allele_check(103000000,"A",
              FakeSeq(104759757),FakeSeq(103000000),
              Duplicated(1759757),FakeLift(-1759757))
        self.assertEqual(why,"HOLD_NO_UNIQUE_CHAIN_MAPPING")

    def test_negative_strand_is_held(self):
        class Reverse(FakeLift):
            def convert_coordinate(self,c,p):
                return [("chr10",p+self.shift,"-",120)]
        pos,why=M.reciprocal_allele_check(103000000,"A",
              FakeSeq(104759757),FakeSeq(103000000),
              Reverse(1759757),FakeLift(-1759757))
        self.assertEqual(why,"HOLD_NONCHR10_OR_NEGATIVE_STRAND")

    def test_original_qtl_valid_se_and_ac_over_an(self):
        x=M.qtl_record(raw_fields())
        self.assertEqual(x["pos38"],103000000)
        self.assertAlmostEqual(x["alt_af"],65/300)
        self.assertEqual(x["se"],0.05)

    def test_original_qtl_maf_not_effect_af(self):
        x=M.qtl_record(raw_fields())
        self.assertNotEqual(x["maf"],x["alt_af"])

    def test_observed_qtl_nonpositive_se_rejected(self):
        d=raw_fields();d["se"]="0"
        self.assertIsNone(M.qtl_record(d))

    def test_invalid_variant_key_rejected(self):
        d=raw_fields();d["variant"]="chr9_103000000_A_C"
        self.assertIsNone(M.qtl_record(d))

    def test_QTL_indels_not_treated_as_valid_SNP(self):
        d=raw_fields()
        d.update({"ref":"AA","variant":"chr10_103000000_AA_C","type":"INDEL"})
        self.assertIsNone(M.qtl_record(d))

    def test_falsified_effect_frequency_is_not_gene_N(self):
        d=raw_fields();d["ac"]="."
        x=M.qtl_record(d)
        self.assertIsNone(x["alt_af"])

    def test_oriented_original_gwas_alt_effect_is_preserved(self):
        s,b,af=M.harmonize("A","C","C","A",-0.23,0.13)
        self.assertEqual(s,"MATCH_ALT_EFFECT")
        self.assertAlmostEqual(b,-0.23)
        self.assertAlmostEqual(af,0.13)

    def test_oriented_original_gwas_ref_effect_sign_is_flipped(self):
        s,b,af=M.harmonize("A","C","A","C",0.23,0.87)
        self.assertEqual(s,"MATCH_REF_EFFECT_SIGN_FLIPPED")
        self.assertAlmostEqual(b,-0.23)
        self.assertAlmostEqual(af,0.13)

    def test_palindromic_G_C_always_hold(self):
        s,b,af=M.harmonize("G","C","C","G",-0.23,0.13)
        self.assertEqual(s,"HOLD_PALINDROMIC_STRAND")
        self.assertIsNone(b)
        self.assertIsNone(af)

    def test_palindromic_A_T_always_hold(self):
        self.assertEqual(M.harmonize("A","T","T","A",1,0.5)[0],
                         "HOLD_PALINDROMIC_STRAND")

    def test_complement_only_not_sign_approved(self):
        s,b,af=M.harmonize("A","C","T","G",0.23,0.13)
        self.assertEqual(s,"HOLD_COMPLEMENT_ONLY_NO_ORIENTATION_APPROVAL")
        self.assertIsNone(b)

    def test_non_matching_allele_pair_hold(self):
        self.assertEqual(M.harmonize("A","C","A","G",0.23,0.13)[0],
                         "HOLD_ALLELE_PAIR_MISMATCH")

    def test_record_group_count_is_thirty_six(self):
        d=json.loads(MANIFEST.read_text(encoding="utf8"))
        self.assertEqual(len(d["groups"]),36)
        self.assertEqual(len({x["gene_id"] for x in d["groups"]}),20)
        self.assertEqual(M.GWAS_ANCESTRIES,("EAS","EUR"))

    def test_new_schema_is_not_formal_coloc(self):
        self.assertEqual(M.PASS_STATUS,
                         "PASS_QTL_GWAS_SOURCE_SNP_HARMONIZATION_AGGREGATES_ONLY")
        self.assertEqual(M.REVERSE_CHAIN_MD5,
                         "ff3031d93792f4cbb86af44055efd903")


if __name__=="__main__":
    unittest.main()
