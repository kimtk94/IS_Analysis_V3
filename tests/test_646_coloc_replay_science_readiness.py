"""Defensive classification of full-source SNP replay vs biological inference."""
import importlib.util
import sys
import unittest
from pathlib import Path

SOURCE=Path(__file__).resolve().parents[1]/"scripts/is/summarize_646_coloc_replay_science_readiness.py"
spec=importlib.util.spec_from_file_location("replay_science",SOURCE)
mod=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=mod
spec.loader.exec_module(mod)

class Test646ScienceQC(unittest.TestCase):
    def source(self):
        r={"locus":"BBJ_IS_L001","dataset_key":"GTEx_V8__Brain_Cerebellum",
            "gene_base":"ENSG00000138675","status":"PASS",
            "maf_diff_fraction":"0.41","maf_diff_gt_0p1":"41",
            "n_qtl_sample_N_distinct":"8","delta_max":"1e-15",
            "new_PP_H4":"0.7","source_file":"test.tsv"}
        master={"locus":r["locus"],"dataset_key":r["dataset_key"],
                "gene_base":r["gene_base"],"nsnps":"100","PP.H4":"0.7"}
        return r,master

    def test_agreement_not_proof(self):
        r,m=self.source()
        result=mod.assess([r],[m])
        self.assertEqual(result[0]["new_replay_status"],"PASS")
        self.assertEqual(result[0]["maf_QC_review"],"MAF_DIFF_GT_0P1_IN_OVER_QUARTER_SNPS")
        self.assertTrue(result[0]["QTL_sample_N_variable_across_SNPs"])
        self.assertEqual(result[0]["scientific_status"],"NO_MATCHED_QTL_LD_NOT_CAUSAL")

    def test_nonmatching_posterior_reject_pass(self):
        r,m=self.source()
        r["new_PP_H4"]="0.9"
        with self.assertRaises(ValueError):
            mod.assess([r],[m])

    def test_unverified_missing_input_is_not_negative(self):
        r,m=self.source()
        r.update({"status":"MISSING_INPUT","new_PP_H4":"","delta_max":"",
                  "maf_diff_fraction":"","maf_diff_gt_0p1":"",
                  "n_qtl_sample_N_distinct":""})
        result=mod.assess([r],[m])
        self.assertEqual(result[0]["maf_QC_review"],"NOT_RECOMPUTED")
        self.assertEqual(result[0]["new_replay_status"],"MISSING_INPUT")

if __name__=="__main__":
    unittest.main()
