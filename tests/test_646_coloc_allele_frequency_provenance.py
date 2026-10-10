"""Fail-closed allele-frequency provenance tests, with explicit AF vs MAF."""
import importlib.util
import sys
import unittest
from pathlib import Path

SCRIPT=Path(__file__).resolve().parents[1]/"scripts/is/audit_646_coloc_allele_frequency_provenance.py"
spec=importlib.util.spec_from_file_location("is_646_af",SCRIPT)
mod=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=mod
spec.loader.exec_module(mod)

def row(**overrides):
    r={
      "match_key":"12:111803962:G:A",
      "variant_id_hg38":"12:111803962:G:A",
      "ref":"G","alt":"A","locus":"BBJ_IS_L003",
      "dataset_key":"GTEx_V8__Artery_Aorta",
      "gene_base":"ENSG00000111275",
      "gwas_eaf":"0.6",
      "gwas_maf":"0.4",
      "eqtl_maf":"0.02",
      "ac":"8","an":"400",
      "gwas_beta":"0.08","eqtl_beta":"0.05",
      "harmonization":"EXACT_REF_ALT_GRCH38",
    }
    r.update(overrides)
    return r

class TestAlleleFrequencyProvenance(unittest.TestCase):
    def test_distinguish_minor_maf_and_effect_allele_frequency(self):
        a=mod.audit_row(row())
        self.assertAlmostEqual(a["qtl_ac_an_ratio"],0.02)
        self.assertAlmostEqual(a["maf_abs_delta"],0.38)
        self.assertGreater(a["eaf_vs_AC_AN_delta"],0.5)
        self.assertFalse(a["palindromic"])

    def test_palindromic_and_near_half(self):
        a=mod.audit_row(row(match_key="12:111803962:A:T",
            variant_id_hg38="12:111803962:A:T",
            ref="A",alt="T",gwas_eaf="0.49",
            gwas_maf="0.49"))
        self.assertTrue(a["palindromic"])
        self.assertTrue(a["palindromic_near_half"])

    def test_reference_pair_mismatch_rejected(self):
        with self.assertRaises(ValueError):
            mod.audit_row(row(ref="C"))
        with self.assertRaises(ValueError):
            mod.audit_row(row(variant_id_hg38="12:111803963:G:A"))

    def test_frequencies_must_be_mathematically_consistent(self):
        with self.assertRaises(ValueError):
            mod.audit_row(row(gwas_maf="0.6"))
        with self.assertRaises(ValueError):
            mod.audit_row(row(eqtl_maf="0.17"))
        with self.assertRaises(ValueError):
            mod.audit_row(row(ac="401"))

    def test_label_is_not_effect_allele_orientation_confirmation(self):
        a=mod.audit_row(row(gwas_eaf="0.9",gwas_maf="0.1"))
        self.assertGreater(a["eaf_vs_AC_AN_delta"],
                           a["eaf_complement_vs_AC_AN_delta"])
        self.assertEqual(set(a).intersection(
          {"inferred_effect_allele_flip","strand_flip_applied"}),set())

if __name__=="__main__":
    unittest.main()
