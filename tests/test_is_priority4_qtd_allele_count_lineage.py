import csv
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_is_priority4_qtd_allele_count_lineage import classify,run

class TestAlleleCountLineage(unittest.TestCase):
    def test_count_classification_not_beta_sign_inference(self):
        self.assertEqual(classify(5,5,20),"AC_SAME")
        self.assertEqual(classify(5,15,20),"AC_COMPLEMENT")
        self.assertEqual(classify(5,7,20),"AC_OTHER")
        self.assertEqual(classify(10,10,20),"AC_SAME")
    def test_real_four_gene_reanalysis_gate_when_data_mounted(self):
        root=Path("/srv/is-analysis/results/is/audits")
        source=root/"legacy_coloc_drive_recovery_20261011_v1/input"
        external=root/"is_priority4_qtd000161_full_cis_universe_20261011_v1"
        if not source.exists() or not all((external/f"BBJ_IS_L00{i}.QTD000161.region.priority4.tsv").exists() for i in (1,2)):
            self.skipTest("Original archived and EBI QTL inputs unavailable")
        with tempfile.TemporaryDirectory() as td:
            m=run(source,external,Path(td))
            self.assertEqual(m["genes"],4)
            self.assertEqual(m["snp_gene_allele_exact_matches"],6242)
            self.assertFalse(m["can_recalculate_archived_coloc_pp_h4"])
            self.assertTrue(m["no_automatic_allele_beta_inversion"])
            self.assertEqual(m["original_candidate_universe"],2225)
            a={x["symbol"]:x for x in m["summaries"]}
            self.assertEqual(a["FGF5"]["exact_variant_rows"],1568)
            self.assertEqual(a["FGF5"]["ac_complement"],423)
            self.assertEqual(a["CALHM2"]["ac_complement"],270)
            for row in a.values():
                self.assertEqual(row["maf_equal_2e6"],row["exact_variant_rows"])
                self.assertEqual(row["an_identical"],row["exact_variant_rows"])
                self.assertEqual(row["ac_other"],0)
                self.assertEqual(row["an_different"],0)
                self.assertFalse(row["beta_effect_direction_automatically_flipped"])
if __name__=="__main__":unittest.main()
