"""Offline synthetic regression of full original file chr10 logOR/CI unit checks.

Each synthetic fixture can use reduced test-specific denominator thresholds;
production calls never override the original complete file / chr10 guards.
"""
from __future__ import annotations
import csv
import gzip
import importlib
import json
import math
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

M=importlib.import_module("scripts.is.audit_p3a_gigastroke_entire_chr10_logOR_scale")

FIELDS=("chromosome","base_pair_location","effect_allele","other_allele",
        "effect_allele_frequency","beta","standard_error","p_value",
        "odds_ratio","ci_lower","ci_upper","rsid")

def rows_for_test(n=15, or_offset=0, ci_offset=0, missing_or=False,
                  missing_ci=False, malformed_beta=False):
    output=[]
    for i in range(n):
        beta=.1 if i%2 else -.03
        se=.02
        row={"chromosome":"10","base_pair_location":str(104_500_000+i),
             "effect_allele":"A","other_allele":"G",
             "effect_allele_frequency":"0.31",
             "beta":"nan" if malformed_beta and i==1 else str(beta),
             "standard_error":str(se),
             "p_value":str(math.erfc(abs(beta/se)/math.sqrt(2))),
             "odds_ratio":"" if missing_or and i==1 else
                         str(math.exp(beta+or_offset)),
             "ci_lower":"" if missing_ci and i==1 else
                        str(math.exp(beta-1.96*se+ci_offset)),
             "ci_upper":"" if missing_ci and i==1 else
                        str(math.exp(beta+1.96*se+ci_offset)),
             "rsid":f"rsSynthetic{i}"}
        output.append(row)
    return output


def write_original_fixture(path,rows):
    with gzip.open(path,"wt",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def scan(rows, ancestry="EAS"):
    with tempfile.TemporaryDirectory(prefix="synthetic_gwas_") as directory:
        p=Path(directory)/"gwas.tsv.gz"
        write_original_fixture(p,rows)
        return M.audit_all_chr10_rows(
            p,ancestry,min_chr10=1,expected_total=False,min_total_rows=1)


class FullGIGASTROKEOriginalScaleControls(unittest.TestCase):
    def test_all_synthetic_chr10_OR_and_CI_valid(self):
        x=scan(rows_for_test())
        self.assertEqual(x["original_total_chr10_rows"],15)
        self.assertTrue(x["all_source_chr10_OR_and_CI_logscale_consistent"])
        self.assertEqual(x["aggregates"]["chr10_rows_OR_and_CI_populated"],15)
        self.assertLess(x["max_abs_beta_minus_ln_original_OR_all_chr10"],1e-12)
        self.assertFalse(x["source_beta_units_approved_for_automated_formal_coloc"])

    def test_EUR_source_synthetic_valid(self):
        x=scan(rows_for_test(),"EUR")
        self.assertTrue(x["all_source_chr10_OR_and_CI_logscale_consistent"])
        self.assertFalse(x["exact_source_file_cohort_case_fraction_verified"])

    def test_or_scale_mismatch_fails_closed(self):
        x=scan(rows_for_test(or_offset=.5))
        self.assertFalse(x["all_source_chr10_OR_and_CI_logscale_consistent"])
        self.assertEqual(x["aggregates"]["chr10_rows_OR_not_exp_beta"],15)
        self.assertEqual(x["source_beta_lnOR_and_CI_consistency_chr10"],
                         "HOLD_SOME_CHR10_SOURCE_ROWS")

    def test_CI_wrong_scale_fails_closed(self):
        x=scan(rows_for_test(ci_offset=.15))
        self.assertFalse(x["all_source_chr10_OR_and_CI_logscale_consistent"])
        self.assertEqual(x["aggregates"][
            "chr10_rows_CI_not_beta_plusminus_1_96_SE"],15)

    def test_optional_missing_OR_is_HOLD_not_imputation(self):
        x=scan(rows_for_test(missing_or=True))
        self.assertEqual(x["aggregates"]["chr10_rows_OR_missing_or_invalid"],1)
        self.assertFalse(x["all_source_chr10_OR_and_CI_logscale_consistent"])

    def test_optional_missing_CI_is_HOLD_not_imputation(self):
        x=scan(rows_for_test(missing_ci=True))
        self.assertEqual(x["aggregates"]["chr10_rows_CI_missing_or_invalid"],1)
        self.assertFalse(x["all_source_chr10_OR_and_CI_logscale_consistent"])

    def test_nonfinite_beta_is_invalid_not_auto_fixed(self):
        x=scan(rows_for_test(malformed_beta=True))
        self.assertEqual(x["aggregates"]["chr10_rows_invalid_beta_SE_p_or_alleles"],1)
        self.assertFalse(x["all_source_chr10_OR_and_CI_logscale_consistent"])

    def test_non_chr10_rows_ignored_without_flipping_origin(self):
        rows=rows_for_test(15)
        rows[0]["chromosome"]="11"
        x=scan(rows)
        self.assertEqual(x["original_total_chr10_rows"],14)
        self.assertTrue(x["all_source_chr10_OR_and_CI_logscale_consistent"])

    def test_numeric_beta_se_vs_p_original_checked(self):
        x=scan(rows_for_test())
        self.assertEqual(x["aggregates"]["chr10_rows_z_p_comparable"],15)
        self.assertEqual(x["aggregates"][
            "chr10_rows_z_p_discrepant_log10_gt_0_1"],0)

    def test_source_outcome_scientific_holds_independent_BBJ_replication(self):
        x=scan(rows_for_test())
        self.assertFalse(x["BBJ_independent_replication_approved"])
        self.assertFalse(x["GWAS_per_SNP_effective_N_verified"])
        self.assertFalse(x["exact_source_file_cohort_case_fraction_verified"])

    def test_does_not_export_per_variant_beta_or_rsid(self):
        x=scan(rows_for_test())
        txt=json.dumps(x)
        self.assertNotIn("rsSynthetic0",txt)
        self.assertNotIn('"beta":',txt)
        self.assertNotIn('"odds_ratio":',txt)
        self.assertFalse(x["all_original_SNP_effect_rows_or_IDs_emitted"])

    def test_strict_publisher_full_source_row_count_on_by_default(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"gwas.tsv.gz"
            write_original_fixture(p,rows_for_test(15))
            with self.assertRaises(ValueError):
                M.audit_all_chr10_rows(p,"EAS")

    def test_missing_original_publisher_scale_header_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"gwas.tsv.gz"
            rows=rows_for_test()
            with gzip.open(p,"wt",newline="") as f:
                fields=[x for x in FIELDS if x!="odds_ratio"]
                w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",
                                 extrasaction="ignore")
                w.writeheader()
                w.writerows(rows)
            with self.assertRaises(ValueError):
                M.audit_all_chr10_rows(p,"EAS",min_total_rows=1,min_chr10=1,
                                       expected_total=False)

    def test_original_assay_manifest_identity_preserved(self):
        p=Path("docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json")
        groups=M.SOURCE.checked_manifest(p)
        self.assertEqual(len(groups),36)
        self.assertEqual(len({x["gene_id"] for x in groups}),20)
        self.assertTrue(all(g["approved_for_new_coloc"] is False for g in groups))

    def test_invalid_ancestry_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"gwas.tsv.gz"
            write_original_fixture(p,rows_for_test())
            with self.assertRaises(ValueError):
                M.audit_all_chr10_rows(p,"AFR",min_chr10=1,
                                       expected_total=False,min_total_rows=1)


if __name__=="__main__":
    unittest.main()
