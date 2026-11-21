"""Fail-closed pinned source metadata N reconciliation tests; no network."""
import unittest
from scripts.ckm_stroke import verify_tpmi_43321_source_obsct_discrepancy as m

def fake():
    n=m.EXPECTED_ROWS
    chrom={str(i):1 for i in range(1,23)}
    chrom["1"]=n-21
    return {
      "schema":"TPMI43321-original-gwas-complete-row-QC-v1",
      "technical_status":"PASS_FULL_SOURCE_SHA256_AND_ALL_ROWS_PARSED_WITH_QC_FLAGS",
      "full_original_compressed_sha256":m.PINNED_SHA,
      "pinned_prior_original_sha256_match":True,
      "original_compressed_bytes":58_302_393,
      "rows_total_in_original":n,
      "pinned_original_columns":["c"]*17,
      "scientific_status":"HOLD","ready_for_ldsc":False,
      "MR_coloc_LDSC_executed":False,
      "by_chr_all_original_rows":chrom,
      "row_level_QC_counts":{key:n for key in
        ("records","fully_numeric_ADD_no_error","P_valid",
         "OR_positive","LOGOR_SE_positive","OBS_CT_valid",
         "freq_numeric","P_vs_Z_comparable")} | {"P_le_5e_8":1},
      "duplicate_variant_test_128bit_hashes":0,
      "original_additive_test_categories":{"ADD":n},
      "original_ERRORCODE_categories":{".":n},
      "source_rows_out_of_coordinate_order":0,
      "original_MAF_bins_from_A1_FREQ":{"MAF_ge_0.05":n},
      "value_min_max":{"min_OBS_CT":204_099,"max_OBS_CT":204_099},
      "original_full_file_deleted_after_scan":True,
    }
class TestPinnedMismatch(unittest.TestCase):
    def test_correct_denominator_mismatch_is_hard_HOLD(self):
        r=m.validate(fake())
        self.assertEqual(r["observed_total_difference"],109810)
        self.assertEqual(r["observed_total_difference_percent_of_page_N"],34.98148)
        self.assertFalse(r["binary_case_control_Neff_approved"])
        self.assertFalse(r["phenotype_N_vs_regression_N_reconciled_by_individual_records"])
        self.assertFalse(r["ready_for_ldsc"])
    def test_mutated_compressed_source_digest_refused(self):
        r=fake();r["full_original_compressed_sha256"]="0"*64
        with self.assertRaisesRegex(ValueError,"integrity"):
            m.validate(r)
    def test_mismatched_22_chr_sum_refused(self):
        r=fake();r["by_chr_all_original_rows"]["1"]-=1
        with self.assertRaisesRegex(ValueError,"22-chromosome"):
            m.validate(r)
    def test_gwas_numeric_invalid_count_refused(self):
        r=fake();r["row_level_QC_counts"]["P_valid"]-=1
        with self.assertRaisesRegex(ValueError,"complete-case"):
            m.validate(r)
    def test_bad_MAF_ledger_refused(self):
        r=fake();r["original_MAF_bins_from_A1_FREQ"]["MAF_ge_0.05"]-=10
        with self.assertRaisesRegex(ValueError,"MAF bins"):
            m.validate(r)
    def test_gwas_technical_pass_not_sufficient_for_Neff(self):
        r=m.validate(fake())
        self.assertFalse(r["SNP_case_specific_N_available"])
        self.assertFalse(r["SNP_control_specific_N_available"])
        self.assertFalse(r["study_overlap_with_GIGASTROKE_excluded_by_individual_ID"])
    def test_sample_n_change_must_be_investigated(self):
        r=fake();r["value_min_max"]["max_OBS_CT"]+=1
        with self.assertRaisesRegex(ValueError,"OBS_CT"):
            m.validate(r)
    def test_fake_publisher_original_test_error_refused(self):
        r=fake();r["original_ERRORCODE_categories"]={"FAIL":m.EXPECTED_ROWS}
        with self.assertRaisesRegex(ValueError,"mutation"):
            m.validate(r)
    def test_not_deleting_original_source_refused(self):
        r=fake();r["original_full_file_deleted_after_scan"]=False
        with self.assertRaisesRegex(ValueError,"not removed"):
            m.validate(r)

if __name__=="__main__":
    unittest.main()
