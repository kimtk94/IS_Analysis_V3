import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_is_priority4_ld_proxy_ancestry_stability import classify,run


class AncestryProxyGateTests(unittest.TestCase):
    def test_four_way_tagging_and_missing_not_negative(self):
        self.assertEqual(classify(".9",".95"),"HIGH_LD_BOTH")
        self.assertEqual(classify(".8",".5"),"HIGH_LD_EAS_ONLY")
        self.assertEqual(classify(".5",".8"),"HIGH_LD_EUR_ONLY")
        self.assertEqual(classify(".1",".2"),"BELOW_HIGH_LD_BOTH")
        self.assertEqual(classify("nan",".9"),"UNINFORMATIVE_BOTH_ANCESTRIES")
        with self.assertRaises(ValueError):
            classify("1.1",".4")

    def test_full_locus_counts_and_noncausal_fields_if_mounted(self):
        base=Path("/srv/is-analysis/results/is/audits")
        src=base/"is_priority4_eur_eas_signed_ld_20261011_v1/IS_PRIORITY4_EUR_EAS_GWAS_LEAD_TAGGING_BY_SNP.tsv"
        original=base/"legacy_coloc_drive_recovery_20261011_v1/input"
        if not src.exists() or not original.is_dir():
            self.skipTest("Original IS research server source unavailable")
        with tempfile.TemporaryDirectory() as tmp:
            x=run(src,original,Path(tmp))
            self.assertEqual(x["reference_snp_rows"],3277)
            self.assertEqual(x["unique_loci"],2)
            self.assertEqual(x["gene_tissue_assays"],4)
            self.assertFalse(x["PP_H4_recomputed"])
            self.assertFalse(x["LD_matrix_directly_used_by_coloc_abf"])
            self.assertEqual(x["original_candidate_gene_universe"],2225)
            regions={r["locus"]:r for r in x["locus_summaries"]}
            self.assertEqual(regions["BBJ_IS_L001"]["material_abs_delta_r2_gt_0p2"],22)
            self.assertEqual(regions["BBJ_IS_L002"]["material_abs_delta_r2_gt_0p2"],72)
            gs={r["symbol"]:r for r in x["gene_summaries"]}
            self.assertEqual(gs["CALHM2"]["qtl_top_proxy_category"],"HIGH_LD_BOTH")
            self.assertEqual(gs["CALHM2"]["n_nominal_QTL_p_lt_1e4"],0)
            self.assertEqual(gs["FGF5"]["n_nominal_QTL_p_lt_1e4"],19)
            self.assertEqual(gs["C4orf22"]["n_nominal_QTL_p_lt_1e4"],1)
            self.assertEqual(gs["NEURL"]["n_nominal_QTL_p_lt_1e4"],16)

if __name__=="__main__":
    unittest.main()
