"""Synthetic negative and numerical tests for source-QTL coloc sdY.est OLS.

No network or genuine SNP effect records are loaded into fixtures.
All original original-qtl 15 IDs derived from versioned 36-source manifest.
"""
from __future__ import annotations

import importlib
import json
import math
import unittest
from pathlib import Path
from unittest.mock import patch

M=importlib.import_module("scripts.is.audit_p3a_rosmap_real_sdY_est_15gene")
Q=M.QTL


def manifest():
    return json.loads(Path(
        "docs/is/evidence/IS_P3A_RESTORED_REAL_36_ASSAY_KEYS_20261122.json"
    ).read_text(encoding="utf-8"))


def candidate_genes():
    return [g for g in manifest()["groups"] if g["QTL_dataset"]==Q.DATASET]


def metadata():
    return {
        g["gene_id"]:{"chromosome":"10",
                      "phenotype_pos":102_000_000+i*1000,
                      "gene_interval_start":102_000_000+i*1000,
                      "gene_interval_end":102_000_200+i*1000,
                      "strand":"+"}
        for i,g in enumerate(candidate_genes())
    }


def region():
    return Q.one_combined_region(metadata(),set(metadata()))


def one_gene_variant(gid,pos,maf,sd=2.4,an="1120",
                     beta="0.01",rsid="rsTEST"):
    maf=float(maf)
    N=560
    se=math.sqrt(sd*sd/(2*N*maf*(1-maf)))
    fields={
        "molecular_trait_id":gid,"chromosome":"10","position":str(pos),
        "ref":"C","alt":"G","variant":f"chr10_{pos}_C_G",
        "ma_samples":"80","maf":str(maf),"pvalue":"0.5",
        "beta":beta,"se":str(se),"type":"SNP","ac":"120","an":an,"r2":"0.99",
        "molecular_trait_object_id":gid,"gene_id":gid,
        "median_tpm":"2","rsid":rsid,
    }
    return "\t".join(fields[c] for c in Q.COLUMNS)+"\n"


def synthetic_lines(with_duplicates=True,an="1120",bad_duplicate=False,
                    leave_out_last=False):
    ret=[]
    genes=candidate_genes()
    if leave_out_last:
        genes=genes[:-1]
    for g in genes:
        gid=g["gene_id"]
        pos=metadata()[gid]["phenotype_pos"]
        maf=[.02,.04,.2,.35]
        for i,v in enumerate(maf):
            ret.append(one_gene_variant(gid,pos+100+i*2,v,an=an,
                                        rsid=f"rsTEST_{i}"))
        if with_duplicates:
            ret.append(one_gene_variant(
                gid,pos+100,maf[0],an=an,
                beta="0.03" if bad_duplicate else "0.01",rsid="rsALIAS"))
    return ret


def mocked_source_test(lines,extra_patches=None):
    expected_counts={g["gene_id"]:4 for g in candidate_genes()}
    with patch.object(M,"EXACT_15_RECORD_COUNTS",expected_counts), \
         patch.object(M,"EXPECTED_SOURCE_ROWS_IN_UNION",len(lines)), \
         patch.object(M,"EXPECTED_TARGET_SOURCE_ROWS",len(lines)), \
         patch.object(M,"MIN_GENE_INPUT_ROWS",2), \
         patch.object(M,"MIN_STRATUM_ROWS",2):
        return M.source_record_qc(lines,manifest(),region(),metadata())


class OriginalColocSDY15GeneSourceTests(unittest.TestCase):
    def test_coloc_R_zero_intercept_OLS_equivalence(self):
        n=560
        sigma=2.4
        data=[(m,math.sqrt(sigma*sigma/(2*n*m*(1-m))))
              for m in [.02,.04,.1,.2,.35,.48]]
        out=M.estimate_coloc_sdY(data,n)
        self.assertAlmostEqual(out["sdY_est_coloc_zero_intercept_OLS"],sigma,places=11)
        self.assertAlmostEqual(out["zero_intercept_OLS_R_squared_diagnostic"],1,places=10)
        self.assertFalse(out["this_is_original_normalized_expression_sdY_direct_measurement"])

    def test_OLS_not_median_per_variant_for_variable_SE(self):
        n=560
        data=[(.02,.05),(.04,.08),(.2,.09),(.35,.03)]
        a=M.estimate_coloc_sdY(data,n)
        self.assertGreater(a["sdY_est_coloc_zero_intercept_OLS"],0)
        self.assertNotAlmostEqual(a["sdY_est_coloc_zero_intercept_OLS"],
                                  a["median_single_SNP_implied_sdY_diagnostic"],
                                  delta=.0001)

    def test_source_AN560_stability_hypothetical_N_scaling(self):
        d=[(m,.2) for m in [.02,.04,.20,.35]]
        s=M.estimate_coloc_sdY(d,560)["sdY_est_coloc_zero_intercept_OLS"]
        for n in (540,580):
            x=M.estimate_coloc_sdY(d,n)["sdY_est_coloc_zero_intercept_OLS"]
            self.assertAlmostEqual(x/s,math.sqrt(n/560),places=12)

    def test_invalid_n_and_nonpositive_maf_and_se_rejected(self):
        for n in (0,1,float("nan"),float("inf"),100_000):
            with self.subTest(n=n),self.assertRaises(ValueError):
                M.estimate_coloc_sdY([(.1,.02),(.2,.03)],n)
        for d in ([(0,.02),(.2,.03)],[(.1,0),(.2,.03)],
                  [(.1,float("nan")),(.2,.02)],
                  [(.6,.02),(.2,.03)]):
            with self.subTest(d=d),self.assertRaises(ValueError):
                M.estimate_coloc_sdY(d,560)

    def test_single_source_variant_not_enough_for_estimator(self):
        with self.assertRaises(ValueError):
            M.estimate_coloc_sdY([(.3,.04)],560)

    def test_synthetic_original_15_genes_recovered_and_estimated(self):
        lines=synthetic_lines()
        a=mocked_source_test(lines)
        self.assertEqual(a["n_genes_real_sdY_OLS_estimated_N560_scenario"],15)
        self.assertEqual(a["real_original_gene_wise_unique_SNP_sum"],60)
        self.assertEqual(a["real_original_target_gene_QTL_source_rows"],75)
        self.assertEqual(a["n_genes_with_both_MAF_strata_estimable"],15)
        for g in a["groups"]:
            self.assertAlmostEqual(g["official_coloc_sdY_estimator_N560_hypothesis"][
                "sdY_est_coloc_zero_intercept_OLS"],2.4,places=10)
            self.assertEqual(g["source_AN_div2_nominal_N_hypothesis"],560)
            self.assertFalse(g["sdY_estimate_approved_for_formal_coloc"])
            self.assertFalse(g["original_expression_matrix_sdY_directly_measured"])

    def test_duplicate_rsid_equivalent_effect_dedup(self):
        a=mocked_source_test(synthetic_lines())
        self.assertEqual(a["real_original_target_gene_QTL_source_rows"],75)
        self.assertEqual(a["real_original_gene_wise_unique_SNP_sum"],60)

    def test_duplicate_rsid_conflicting_effect_fails_closed(self):
        with self.assertRaises(ValueError):
            mocked_source_test(synthetic_lines(bad_duplicate=True))

    def test_qtl_one_gene_N1118_fails_not_silently_assumed_560(self):
        with self.assertRaises(ValueError):
            mocked_source_test(synthetic_lines(an="1118"))

    def test_missing_a_target_gene_source_rows_fails(self):
        with self.assertRaises(ValueError):
            mocked_source_test(synthetic_lines(leave_out_last=True))

    def test_original_total_source_row_drift_fails(self):
        lines=synthetic_lines()
        with patch.object(M,"EXPECTED_SOURCE_ROWS_IN_UNION",len(lines)+1), \
             patch.object(M,"EXPECTED_TARGET_SOURCE_ROWS",len(lines)+1), \
             patch.object(M,"EXACT_15_RECORD_COUNTS",
                          {g["gene_id"]:4 for g in candidate_genes()}):
            with self.assertRaises(ValueError):
                M.source_record_qc(lines,manifest(),region(),metadata())

    def test_original_per_gene_variant_count_drift_fails(self):
        with patch.object(M,"EXACT_15_RECORD_COUNTS",
                          {g["gene_id"]:5 for g in candidate_genes()}), \
             patch.object(M,"EXPECTED_SOURCE_ROWS_IN_UNION",75), \
             patch.object(M,"EXPECTED_TARGET_SOURCE_ROWS",75), \
             patch.object(M,"MIN_GENE_INPUT_ROWS",2):
            with self.assertRaises(ValueError):
                M.source_record_qc(synthetic_lines(),manifest(),region(),metadata())

    def test_privacy_original_effects_not_exported(self):
        data=mocked_source_test(synthetic_lines())
        txt=json.dumps(data)
        self.assertNotIn("rsTEST",txt)
        self.assertNotIn("chr10_102000100_C_G",txt)
        self.assertNotIn('"beta":',txt)
        self.assertTrue(data["no_variant_effect_rows_or_sample_IDs_exported"])

    def test_R_estimator_lineage_pinned_original_blob(self):
        self.assertEqual(M.AN_DIV2_PROVISIONAL_N,560)
        self.assertEqual(M.PINNED_SOURCE_INDEX_SHA256,
          "6cbf12826e19a61594662972a1205ba38cdd3071984e51c24261bff8ace9806a")
        self.assertEqual(M.PINNED_SOURCE_RANGE_SHA256,
          "5a54c1faadb0530ecdca37497f82f2e28f1003c000d7a24bb9fa8422de3fa2f1")

    def test_historical_source_window_not_reclassified_original_R7(self):
        self.assertFalse(region()["R7_actual_qtlmap_run_cis_parameter_independently_verified"])
        self.assertTrue(region()["all_15_historical_phenotype_windows_contained"])

    def test_sort_odd_even_split_not_independent_replication(self):
        data=mocked_source_test(synthetic_lines())
        self.assertTrue(all(x["odd_even_split_is_NOT_independent_due_to_LD"]
                            for x in data["groups"]))

    def test_remaining_21_GTEx_groups_are_not_miscounted(self):
        x=manifest()
        self.assertEqual(len(x["groups"]),36)
        self.assertEqual(len([g for g in x["groups"]
                              if g["QTL_dataset"]!="QTD000434"]),21)


if __name__=="__main__":
    unittest.main()
