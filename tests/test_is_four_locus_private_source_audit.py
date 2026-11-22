"""Negative scientific-claim/source-drift tests for legacy BBJ four locus audit."""
import copy
import json
from pathlib import Path
import unittest

from importlib import import_module

module = import_module("scripts.is.audit_four_locus_private_inputs")
GENES = module.GENES
profile = module.profile
validate_snapshot = module.validate_snapshot

SNAPSHOT = Path("docs/is/evidence/IS_FOUR_LOCUS_PRIVATE_INPUT_QC_20261122.json")


def frozen():
    return json.loads(SNAPSHOT.read_text(encoding="utf-8"))


def source_rows():
    return (
        "snp\tgwas_beta\tgwas_se\tgwas_maf\teqtl_beta\teqtl_se\teqtl_maf\tqtl_n\n"
        "4:101:A:T\t0.1\t0.02\t0.1\t-0.2\t0.05\t0.3\t175.0\n"
        "4:102:C:T\t-0.1\t0.03\t0.4\t0.2\t0.05\t0.42\t175.0\n"
    ).encode("ascii")


class FourLocusInputTests(unittest.TestCase):
    def test_real_inspected_snapshot_is_technical_only(self):
        x = validate_snapshot(frozen())
        self.assertEqual(x["technical_status"], "PASS_FROZEN_SOURCE_AGGREGATE_CONTRACT_ONLY")
        self.assertEqual(x["historical_GWAS_eQTL_input_rows_verified_in_prior_local_read"], 7969)
        self.assertFalse(x["cohort_matched_LD_verified"])
        self.assertEqual(x["replicated_MR_estimates"], 0)
        self.assertEqual(x["independently_recomputed_coloc_results"], 0)
        self.assertEqual(x["scientific_status"], "HOLD")

    def test_source_hash_tamper_rejected(self):
        s = frozen()
        s["genes"][0]["source_sha256"] = "0"*64
        with self.assertRaisesRegex(ValueError, "source row/hash"):
            validate_snapshot(s)

    def test_raw_n_and_palindromes_are_not_scientific_approval(self):
        s = frozen()
        s["genes"][1]["per_variant_effect_n_independently_verified"] = True
        with self.assertRaisesRegex(ValueError, "unjustified MR"):
            validate_snapshot(s)
        s = frozen()
        s["genes"][2]["palindromic_variant_count"] = 0
        with self.assertRaisesRegex(ValueError, "source row/hash"):
            validate_snapshot(s)

    def test_matrix_size_is_not_signed_ld_validation(self):
        s = frozen()
        s["genes"][0]["reference_ld_file_size_only"]["binary_content_digest_and_snp_order_verified"] = True
        with self.assertRaisesRegex(ValueError, "LD size-only"):
            validate_snapshot(s)

    def test_genome_build_cannot_be_guessed(self):
        s = frozen()
        s["genes"][0]["source_genome_build_independently_verified"] = True
        with self.assertRaisesRegex(ValueError, "unjustified MR"):
            validate_snapshot(s)

    def test_separate_ebi_abf_not_independent_replication(self):
        s = frozen()
        s["separate_source_model"]["historical_GTEx_model_independent_replication"] = True
        with self.assertRaisesRegex(ValueError, "replication"):
            validate_snapshot(s)

    def test_unexpected_scientific_pass_rejected(self):
        s = frozen()
        s["scientific_status"] = "PASS"
        with self.assertRaisesRegex(ValueError, "scientific causal"):
            validate_snapshot(s)

    def test_input_synthetic_valid_rows(self):
        result = profile(source_rows(), 4)
        self.assertEqual(result["source_rows"], 2)
        self.assertEqual(result["unique_variant_keys"], 2)
        self.assertEqual(result["palindromic_variant_count"], 1)
        self.assertEqual(result["abs_gwas_eqtl_maf_gap_gt_0_2_count"], 0)
        self.assertEqual(result["distinct_qtl_n_values"], [175.0])

    def test_duplicate_snp_rejected(self):
        row = source_rows().decode().splitlines()
        data = ("\n".join([row[0], row[1], row[1]])+"\n").encode()
        with self.assertRaisesRegex(ValueError, "duplicate SNP"):
            profile(data, 4)

    def test_negative_se_rejected(self):
        data = source_rows().replace(b"0.02", b"-0.02", 1)
        with self.assertRaisesRegex(ValueError, "standard error"):
            profile(data, 4)

    def test_wrong_chromosome_rejected(self):
        data = source_rows().replace(b"4:101", b"12:101", 1)
        with self.assertRaisesRegex(ValueError, "coordinate"):
            profile(data, 4)

    def test_maf_not_silently_flipped(self):
        data = source_rows().replace(b"0.42", b"0.71")
        with self.assertRaisesRegex(ValueError, "MAF out of bounds"):
            profile(data, 4)

    def test_four_locus_counts_and_ld_shape_pinned(self):
        self.assertEqual(sum(x[1] for x in GENES.values()), 7969)
        self.assertTrue(all(x[1]**2*8 == x[7] for x in GENES.values()))


if __name__ == "__main__":
    unittest.main()
