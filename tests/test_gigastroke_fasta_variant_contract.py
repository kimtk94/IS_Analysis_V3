"""Regression coverage: historical GIGASTROKE IDs are not necessarily FASTA REF:ALT."""
import importlib.util
import unittest
from pathlib import Path

FILE=(Path(__file__).resolve().parents[1]/
      "scripts/is/gigastroke_fasta_variant_contract.py")
spec=importlib.util.spec_from_file_location("giga_contract",FILE)
tool=importlib.util.module_from_spec(spec)
import sys
sys.modules[spec.name]=tool
spec.loader.exec_module(tool)


class TestGigaContract(unittest.TestCase):
    def test_rs1229984_swapped_id_and_effect_ref(self):
        r=tool.normalize_record({
            "chr":"4","pos":"100239319","build":"GRCh37",
            "effect_allele":"T","other_allele":"C",
            "beta":"-0.0124","se":"0.0137","p":"0.3659","eaf":"0.7546",
            "variant_id":"4:100239319:C:T"
        },"T")
        self.assertEqual(r.variant_id,"4:100239319:T:C")
        self.assertEqual(r.historical_id_status,"OTHER_EFFECT_ORDER")
        self.assertEqual(r.effect_orientation,"EFFECT_REF_FLIPPED")
        self.assertAlmostEqual(r.beta_alt,0.0124)
        self.assertAlmostEqual(r.eaf_alt,0.2454)

    def test_rs671_direct_id(self):
        r=tool.normalize_record({
            "chr":"12","pos":"112241766","build":"GRCh37",
            "effect_allele":"A","other_allele":"G",
            "beta":"-0.1514","eaf":"0.2311","variant_id":"12:112241766:G:A"
        },"G")
        self.assertEqual(r.variant_id,"12:112241766:G:A")
        self.assertEqual(r.historical_id_status,"FASTA_REF_ALT")
        self.assertEqual(r.effect_orientation,"EFFECT_ALT")
        self.assertAlmostEqual(r.beta_alt,-0.1514)

    def test_bad_alleles_raise(self):
        with self.assertRaises(tool.VariantContractError):
            tool.normalize_record({"chr":"4","pos":"100","effect_allele":"A",
                                   "other_allele":"C","beta":"0.1"},"T")

    def test_historical_id_unknown_raises(self):
        with self.assertRaises(tool.VariantContractError):
            tool.normalize_record({"chr":"4","pos":"100","effect_allele":"A",
                                   "other_allele":"C","beta":"0.1",
                                   "variant_id":"4:100:G:T"},"C")

    def test_reject_wrong_build(self):
        with self.assertRaises(tool.VariantContractError):
            tool.normalize_record({"chr":"4","pos":"100","build":"GRCh38",
                                   "effect_allele":"A","other_allele":"C","beta":"0.1"},"C")

    def test_reverse_frequency_when_ref_is_effect(self):
        r=tool.normalize_record({"chr":"4","pos":"1","effect_allele":"A",
                                 "other_allele":"G","beta":"0.25",
                                 "effect_allele_frequency":"0.13"},"A")
        self.assertEqual(r.variant_id,"4:1:A:G")
        self.assertAlmostEqual(r.beta_alt,-0.25)
        self.assertAlmostEqual(r.eaf_alt,0.87)

if __name__=="__main__": unittest.main()
