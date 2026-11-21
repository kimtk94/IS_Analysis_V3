"""Fail-closed tests for TPMI original SAIGE phenotype 433.21 sample fields."""
import unittest
from collections import Counter
from scripts.ckm_stroke.qc_tpmi_43321_saige_full_casecontrol import (
    EXPECTED_HEADER,parse_n,record_fields,
)

def synthetic_row():
    out=["1","102030","1:102030_A/G","A","G","1000","0.003",
         "0.9","0.12","0.04","3.0","0.005","0.0027","0.003",
         "TRUE","0.003","0.002","9249","304660","PASS"]
    return out

class SourceContract(unittest.TestCase):
    def test_header_has_exact_casecontrol_indexes(self):
        self.assertEqual(EXPECTED_HEADER[17:19],("N_case","N_ctrl"))
        self.assertEqual(EXPECTED_HEADER[8:10],("BETA","SE"))
    def test_reject_noninteger_casecounts(self):
        for v in ("","NA","-1","0","9249.5","1e3"):
            self.assertIsNone(parse_n(v))
    def test_known_source_valid_casecontrol(self):
        counts,pairs,chrom=Counter(),Counter(),Counter()
        record_fields(synthetic_row(),counts,pairs,chrom)
        self.assertEqual(counts["records"],1)
        self.assertEqual(pairs[(9249,304660)],1)
        self.assertEqual(counts["finite_beta_positive_se_valid_p"],1)
        self.assertEqual(counts["matches_phenotype_page_N"],1)
    def test_source_samplecount_change_is_flag_not_silent_replacement(self):
        r=synthetic_row()
        r[17]="8500"
        counts,pairs,chrom=Counter(),Counter(),Counter()
        record_fields(r,counts,pairs,chrom)
        self.assertEqual(counts["records"],1)
        self.assertEqual(counts["matches_phenotype_page_N"],0)
        self.assertEqual(pairs[(8500,304660)],1)
    def test_source_qc_fail_must_not_be_approved(self):
        r=synthetic_row()
        r[19]="FAIL"
        counts,pairs,chrom=Counter(),Counter(),Counter()
        record_fields(r,counts,pairs,chrom)
        self.assertEqual(counts["source_explicit_QC_FAIL"],1)
        self.assertEqual(counts["QC_PASS_valid_N_beta_SE_P"],0)
        self.assertEqual(counts["finite_beta_positive_se_valid_p"],1)
    def test_source_qc_pass_can_be_numerically_eligible(self):
        r=synthetic_row()
        counts,pairs,chrom=Counter(),Counter(),Counter()
        record_fields(r,counts,pairs,chrom)
        self.assertEqual(counts["source_explicit_QC_PASS"],1)
        self.assertEqual(counts["QC_PASS_valid_N_beta_SE_P"],1)
    def test_unknown_source_qc_not_auto_pass(self):
        r=synthetic_row()
        r[19]="CUSTOM"
        counts,pairs,chrom=Counter(),Counter(),Counter()
        record_fields(r,counts,pairs,chrom)
        self.assertEqual(counts["source_QC_UNRECOGNIZED"],1)
        self.assertEqual(counts["QC_PASS_valid_N_beta_SE_P"],0)
    def test_invalid_p_is_counted_not_approved(self):
        r=synthetic_row()
        r[12]="1.7"
        counts,pairs,chrom=Counter(),Counter(),Counter()
        record_fields(r,counts,pairs,chrom)
        self.assertEqual(counts["invalid_beta_se_p"],1)
        self.assertEqual(counts["finite_beta_positive_se_valid_p"],0)
    def test_invalid_se_is_counted_not_approved(self):
        r=synthetic_row()
        r[9]="0"
        counts,pairs,chrom=Counter(),Counter(),Counter()
        record_fields(r,counts,pairs,chrom)
        self.assertEqual(counts["invalid_beta_se_p"],1)
    def test_invalid_casecontrol_does_not_pass(self):
        r=synthetic_row()
        r[17]="NA"
        counts,pairs,chrom=Counter(),Counter(),Counter()
        record_fields(r,counts,pairs,chrom)
        self.assertEqual(counts["missing_or_noninteger_N_case_ctrl"],1)
        self.assertEqual(len(pairs),0)

if __name__=="__main__":
    unittest.main()
