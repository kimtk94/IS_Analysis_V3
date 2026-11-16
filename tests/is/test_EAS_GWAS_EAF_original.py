#!/usr/bin/env python3
"""Synthetic EAS GWAS outcome EAF and 22-chromosome denominator guardrail tests."""
import gzip
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"is"))
import audit_EAS_GWAS_EAF_original as q

HEADER=("chromosome\tbase_pair_location\teffect_allele_frequency\t"
        "beta\tstandard_error\tp_value\teffect_allele\tother_allele\n")


class OutcomeFrequencyTests(unittest.TestCase):
    def test_outcome_eaf_numeric_and_missing(self):
        self.assertIsNone(q.classify_freq("NA"))
        self.assertIsNone(q.classify_freq("."))
        self.assertEqual(q.classify_freq("0.8"),0.8)
        with self.assertRaises(ValueError):
            q.classify_freq("1.5")
        with self.assertRaises(ValueError):
            q.classify_freq("-0.4")

    def test_minor_frequency_symmetry_and_bin(self):
        self.assertEqual(q.maf_bin(0.001),"0_001_to_0_01")
        self.assertEqual(q.maf_bin(0.49),"0_4_to_0_5")
        self.assertEqual(min(0.9,1-0.9),min(0.1,1-0.1))

    def sample(self,missing_eaf=False):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/"outcome.tsv.gz"
            with gzip.open(path,"wt",encoding="utf-8") as f:
                f.write(HEADER)
                for c in range(1,23):
                    f.write(f"{c}\t123456\t0.15\t0.1\t0.02\t0.01\tA\tG\n")
                    f.write(f"{c}\t123457\t{'NA' if missing_eaf else '0.49'}\t0.1\t0.02\t5e-9\tA\tT\n")
            with patch.object(q,"EXPECTED_ROWS",44):
                result=q.audit_original_eaf(path)
            return result

    def test_all_22_original_autosomes_conserved(self):
        d=self.sample()
        self.assertEqual(d["rows"]["original_autosomal_rows"],44)
        self.assertEqual(d["rows"]["eaf_valid"],44)
        self.assertEqual(d["rows"]["palindromic_maf_ge_0_4_frequency_ambiguous"],22)
        self.assertEqual(d["rows"]["gwas_significant_eaf_valid_rows"],22)
        self.assertFalse(d["maf_is_from_immuNexUT_donors"])
        self.assertFalse(d["palindromic_frequency_resolved"])

    def test_missing_eaf_not_imputed(self):
        d=self.sample(True)
        self.assertEqual(d["rows"]["eaf_missing"],22)
        self.assertEqual(d["rows"]["eaf_valid"],22)
        self.assertEqual(d["rows"]["palindromic_eaf_valid"],0)

    def test_wrong_original_row_count_must_fail(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/"small.gz"
            with gzip.open(path,"wt") as f:
                f.write(HEADER+"1\t123456\t0.5\t0.1\t0.02\t0.1\tA\tG\n")
            with self.assertRaisesRegex(ValueError,"22-autosome"):
                q.audit_original_eaf(path)

    def test_per_cell_qtl_n_never_inferred(self):
        d=self.sample()
        self.assertEqual(d["n_ancestry_specific_sample_frequency_source"],"OUTCOME_GWAS_EAS_NOT_EQTL")
        self.assertEqual(d["maf_reference_for_coloc"],"NOT_YET_VALIDATED")
        self.assertFalse(d["GWAS_to_eQTL_effect_direction_harmonized"])

if __name__=="__main__":
    unittest.main()
