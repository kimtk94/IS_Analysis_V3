"""Fail-closed original 36 x two-ancestry whole-cis overlap reconciliation tests.

Fixtures construct all 36 actual candidate gene-tissue IDs with entirely
synthetic counts. Never use simulated SNP effects as original observations.
"""
import copy
import importlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

M=importlib.import_module("scripts.is.audit_p3a_36x2_original_fullcis_gwas_coverage")


def original_manifest():
    return json.loads(Path(
        "docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json"
    ).read_text(encoding="utf8"))


def fixture():
    manifest=original_manifest()
    s={"rows":[],"n_real_assay_groups":36,
       "n_formal_coloc_abf_executed":0,
       "matched_nonpal_summary":{"EAS_match":108,"EAS_no_GWAS_pos":108,
                                 "EUR_match":172,"EUR_no_GWAS_pos":36},
       "GWAS_ancestries":{
           a:{"original_publisher_accession":accession,
              "original_publisher_full_file_MD5_SHA256_bytes_pass":True,
              "GWAS_target_invalid_numeric_or_allele_rows":0,
              "GWAS_target_duplicate_locus_rows":0}
           for a,accession in M.ANCESTRIES.items()}}
    c={"groups":[],"n_gene_tissue_groups_matching_full_cis_tested_variant_COUNT":36,
       "n_formal_coloc_abf_executed":0}
    n={"groups":[],"n_real_gene_tissue_assays":36,"n_formal_coloc_abf_executed":0}
    discrepancy_genes=[g for g in manifest["groups"] if g["QTL_dataset"]=="QTD000434"][:8]
    conflicts={g["gene_id"] for g in discrepancy_genes}
    for row in manifest["groups"]:
        basic={k:row[k] for k in ("QTL_dataset","gene_id","gene_symbol")}
        s["rows"].append({
            **basic,
            "QTL_source_indexed_nonpal_valid":6,
            "QTL_uniquely_reciprocally_lifted_REF_valid":6,
            "QTL_liftover_blockers":{},
            "per_ancestry":{
                "EAS":{
                    "dataset":"GCST90104545",
                    "counts_by_signed_GWAS_allele_match":{
                        "MATCH_ALT_EFFECT":2,
                        "MATCH_REF_EFFECT_SIGN_FLIPPED":1,
                        "NO_GWAS_POSITION":3,
                        "AF_ABS_DELTA_GT_0_2_DESCRIPTIVE":1},
                    "n_nonpal_QTL_GWAS_recipient_allele_matches":3,
                    "n_valid_QTL_ALT_AF_vs_GWAS_ALT_EAF_pairs":3},
                "EUR":{
                    "dataset":"GCST90104540",
                    "counts_by_signed_GWAS_allele_match":{
                        "MATCH_ALT_EFFECT":3,
                        "MATCH_REF_EFFECT_SIGN_FLIPPED":1 if row["gene_id"] in conflicts and row["QTL_dataset"]=="QTD000434" else 2,
                        "NO_GWAS_POSITION":1,
                        **({"HOLD_ALLELE_PAIR_MISMATCH":1}
                           if row["gene_id"] in conflicts and row["QTL_dataset"]=="QTD000434" else {})},
                    "n_nonpal_QTL_GWAS_recipient_allele_matches":4 if row["gene_id"] in conflicts and row["QTL_dataset"]=="QTD000434" else 5,
                    "n_valid_QTL_ALT_AF_vs_GWAS_ALT_EAF_pairs":4 if row["gene_id"] in conflicts and row["QTL_dataset"]=="QTD000434" else 5}
            }
        })
        c["groups"].append({
            **basic,
            "nominal_distinct_all_variant_records":10,
            "nominal_strict_nonpal_SNP_source_records":6,
            "original_variant_class_counts":{"SNP":8,"INDEL":2},
            "published_permutation_vs_recovered_nominal_unique_agree":True,
            "source_original_cis_n_variants_from_permutation":10})
        n["groups"].append({
            **basic,
            "QTL_source_variant_N_and_residual_proxy":{
                "original_strict_nonpal_SNPs":6,
                "positive_even_genotyped_AN_count":6,
                "invalid_AN_AC_count":0,
                "AN_over_2_higher_than_catalogue_study_N":0,
                "min_observed_SNP_N":200,
                "max_observed_SNP_N":205}})
    datasets={key:{**val,"all_cis":10*val["n"],"strict":6*val["n"]}
              for key,val in M.DATASETS.items()}
    expected={"cis":360,"strict_nonpal":216,
              "EAS_signed":108,"EUR_signed":172,
              "EAS_absent":108,"EUR_absent":36,"EUR_conflict":8}
    return manifest,s,c,n,datasets,expected


def run_synthetic(change=None):
    a,s,c,n,ds,e=fixture()
    if change:
        change(a,s,c,n)
    with patch.object(M,"DATASETS",ds),patch.object(M,"EXPECTED_TOTAL",e):
        return M.compute_36x2(a,s,c,n)


class FullCisOriginal36x2CoverageTests(unittest.TestCase):
    def test_all_real_36_gene_assay_keys_x_2_in_fixture(self):
        a,_,_,_,_,_=fixture()
        self.assertEqual(len(a["groups"]),36)
        self.assertEqual(len({(x["QTL_dataset"],x["gene_id"]) for x in a["groups"]}),36)
        self.assertEqual(len({x["gene_id"] for x in a["groups"]}),20)

    def test_original_pinned_sha256_is_not_mutable_reference(self):
        self.assertEqual(M.INPUTS["reciprocal"]["sha256"],
          "53b6d71229ed98ab79f404194fa80e94e95bd838fa119ba2961ef234f1308ba8")
        self.assertEqual(M.INPUTS["cis"]["run_id"],45036525615)
        self.assertEqual(M.INPUTS["qtl_n"]["run_id"],45038801844)

    def test_original_source_artifact_digest_conflict_is_fatal(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"source.json"
            p.write_text('{"schema":"fake"}',encoding="utf8")
            with self.assertRaises(ValueError):
                M.load_verified_artifact(p,"reciprocal")

    def test_full_36x2_synthetic_partition_exact(self):
        d=run_synthetic()
        self.assertEqual(d["n_assay_ancestry_rows"],72)
        self.assertEqual(d["aggregate_by_ancestry"]["EAS"]["matched"],108)
        self.assertEqual(d["aggregate_by_ancestry"]["EUR"]["matched"],172)
        self.assertEqual(d["aggregate_by_ancestry"]["EUR"]["discordant"],8)
        self.assertEqual(d["aggregate_by_ancestry"]["EAS"]["original_all_cis"],360)
        self.assertEqual(d["formal_coloc_abf_executed"],0)
        self.assertEqual(d["formal_coloc_abf_eligible"],0)

    def test_matched_coverage_two_denominators_are_not_equated(self):
        r=run_synthetic()
        x=next(a for a in r["rows"] if a["ancestry"]=="EAS")
        self.assertEqual(x["fraction_QTL_strict_nonpal_with_GWAS_match"],0.5)
        self.assertEqual(x["fraction_QTL_ALL_original_cis_with_nonpal_GWAS_match"],0.3)
        self.assertEqual(
            x["source_cis_variants_excluded_before_GWAS_matching"][
                "source_INDEL_all_cis_not_in_nonpal_subset"],2)

    def test_AF_difference_not_added_to_exclusive_missing_counts(self):
        d=run_synthetic()
        x=next(a for a in d["rows"] if a["ancestry"]=="EAS")
        self.assertEqual(x["original_QTL_vs_GWAS_ALT_AF_delta_gt_0_2_DESCRIPTIVE_nonexclusive"],1)
        self.assertEqual(d["aggregate_by_ancestry"]["EAS"]["matched"]+
                         d["aggregate_by_ancestry"]["EAS"]["absent"],216)

    def test_EUR_8_original_gene_assay_mismatches_retained(self):
        r=run_synthetic()
        self.assertEqual(len(r["EUR_only_eight_discordant_gene_assay_rows"]),8)
        self.assertTrue(all(x["QTL_dataset"]=="QTD000434"
                            for x in r["EUR_only_eight_discordant_gene_assay_rows"]))

    def test_wrong_GWAS_allele_direction_breaks_partition(self):
        def bad(a,s,c,n):
            s["rows"][0]["per_ancestry"]["EAS"][
                "counts_by_signed_GWAS_allele_match"]["MATCH_ALT_EFFECT"]+=1
        with self.assertRaises(ValueError):
            run_synthetic(bad)

    def test_changing_source_cis_nvariants_breaks_accounting(self):
        def bad(a,s,c,n):
            c["groups"][0]["nominal_distinct_all_variant_records"]+=1
        with self.assertRaises(ValueError):
            run_synthetic(bad)

    def test_nonpal_QTL_liftover_REF_failure_not_ignored(self):
        def bad(a,s,c,n):
            s["rows"][0]["QTL_liftover_blockers"]={"HOLD_GRCh37_REFERENCE_MISMATCH":1}
        with self.assertRaises(ValueError):
            run_synthetic(bad)

    def test_QTL_sample_AN_not_accepted_as_model_N(self):
        result=run_synthetic()
        self.assertFalse(result["original_gene_specific_expression_model_N_verified"])
        self.assertFalse(result["original_transformed_expression_phenotype_sdY_directly_verified"])
        self.assertTrue(all(x["original_SNP_genotype_N_is_not_approved_regression_N"]
                            for x in result["rows"]))

    def test_CC_case_fraction_s_not_false_hard_blocker(self):
        d=run_synthetic()
        self.assertIs(d["original_binary_cc_beta_varbeta_s_required"],False)
        self.assertTrue(d["original_gwas_chr10_beta_log_OR_source_verified_prior_725707_rows"])

    def test_duplicate_assay_key_fails(self):
        def bad(a,s,c,n):
            a["groups"][1]=copy.deepcopy(a["groups"][0])
        with self.assertRaises(ValueError):
            run_synthetic(bad)

    def test_discordant_EUR_source_record_not_silently_flipped(self):
        d=run_synthetic()
        flagged=[x for x in d["rows"]
                 if x["ancestry"]=="EUR" and x["source_GWAS_allele_discordance_HOLD"]]
        self.assertEqual(len(flagged),8)
        self.assertTrue(all(x["original_allele_mismatch_categories"]==
                            {"HOLD_ALLELE_PAIR_MISMATCH":1} for x in flagged))

    def test_original_permutation_group_count_not_automatically_eligible(self):
        d=run_synthetic()
        self.assertEqual(d["n_permutation_vs_nominal_variant_count_exact_groups"],36)
        self.assertTrue(all(x["formal_coloc_abf_eligible"] is False for x in d["rows"]))

    def test_summary_has_no_per_SNP_effects(self):
        d=run_synthetic()
        txt=json.dumps(d)
        self.assertNotIn('"beta":',txt)
        self.assertNotIn('"se":',txt)
        self.assertNotIn('"pvalue":',txt)
        self.assertTrue(d["no_per_SNP_effect_values_or_patient_ids_exported"])

    def test_wrong_original_file_ancestry_prevents_match(self):
        def bad(a,s,c,n):
            s["GWAS_ancestries"]["EAS"]["original_publisher_accession"]="GCST90104540"
        with self.assertRaises(ValueError):
            run_synthetic(bad)


if __name__=="__main__":
    unittest.main()
