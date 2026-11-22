"""Negative and reference consistency tests for archived 15-source ABF batch."""
from __future__ import annotations
import copy
import importlib
import json
from pathlib import Path
import unittest

m = importlib.import_module("scripts.is.replay_legacy_15_abf_sufficient")
SOURCE = Path("docs/is/evidence/IS_LEGACY_15_GROUP_COLABF_SUFFICIENT_20261122.json")

def raw():
    return json.loads(SOURCE.read_text(encoding="utf-8"))

class OriginalSourceAggregateTests(unittest.TestCase):
    def test_baseline_all_fifteen_replayed(self):
        full=m.reconcile(raw(),0,1)
        self.assertEqual(full["panel_count"],15)
        self.assertEqual(full["scenario_count"],225)
        self.assertEqual(sum(x["original_source_SNP_rows"] for x in full["panels"]),29661)
        self.assertEqual(len({x["gene"] for x in full["panels"]}),11)
        self.assertTrue(all(x["scientific_status"]=="EXPLORATORY_HOLD" for x in full["panels"]))

    def test_parallel_three_shards_reconcile(self):
        shards=[m.reconcile(raw(),i,3) for i in range(3)]
        self.assertEqual([x["panel_count"] for x in shards],[5,5,5])
        self.assertEqual(sum(x["scenario_count"] for x in shards),225)
        self.assertEqual(len({x["panel"] for y in shards for x in y["panels"]}),15)

    def test_sensitivity_both_priors_and_sample_scalar(self):
        full=m.reconcile(raw(),0,1)
        for panel in full["panels"]:
            self.assertEqual(len(panel["sensitivity_grid"]),15)
            self.assertTrue(panel["min_H4"]<=panel["baseline"]["H4"]<=panel["max_H4"])
            self.assertGreater(panel["max_H4"]-panel["min_H4"],0)
            self.assertTrue(all(abs(sum(p[x] for x in ("H0","H1","H2","H3","H4"))-1)<1e-10
                                for p in panel["sensitivity_grid"]))

    def test_reference_FGF5_and_CALHM2_replay(self):
        panels={x["panel"]:x for x in m.reconcile(raw())["panels"]}
        self.assertAlmostEqual(panels["P01"]["baseline"]["H4"],0.77900427,places=6)
        self.assertAlmostEqual(panels["P06"]["baseline"]["H4"],0.7979356,places=6)
        self.assertGreater(panels["P06"]["baseline"]["H4"],panels["P01"]["baseline"]["H4"])

    def test_ALDH2_multiple_tissues_not_earlier_single_maximum(self):
        panels=m.reconcile(raw())["panels"]
        aldh=[x for x in panels if x["gene"]=="ALDH2"]
        self.assertEqual(len(aldh),3)
        self.assertTrue(all(x["baseline"]["H4"]<0.001 for x in aldh))
        self.assertTrue(all(abs(x["baseline"]["H4"]-0.0790399)>0.05 for x in aldh))

    def test_no_source_model_approval(self):
        r=raw()
        r["study_N_and_alleles_approved"]=True
        with self.assertRaisesRegex(ValueError,"provenance"):
            m.reconcile(r)

    def test_coloc_h3_computation_must_exclude_same_causal_variant(self):
        with self.assertRaisesRegex(ValueError,"H3"):
            m.posterior([1,1,3],1e-5)
        pp=m.posterior([1,1,0],1e-5)
        self.assertAlmostEqual(sum(pp),1,places=12)

    def test_source_SNP_row_hash_drift_rejected(self):
        r=raw()
        r["panels"][0][5]="not_a_sha256"
        with self.assertRaisesRegex(ValueError,"hash"):
            m.reconcile(r)

    def test_missing_original_panels_are_fail_closed(self):
        r=raw()
        r["panels"].pop()
        with self.assertRaisesRegex(ValueError,"15"):
            m.reconcile(r)

    def test_gene_group_distinctness_preserved(self):
        r=raw()
        r["panels"][2][0]="P02"
        with self.assertRaisesRegex(ValueError,"identity"):
            m.reconcile(r)

    def test_genomewide_cohort_approval_not_created_by_math(self):
        r=m.reconcile(raw())
        self.assertFalse(r["raw_source_files_re_read_in_ci"])
        self.assertFalse(r["R_coloc_executed"])
        self.assertFalse(r["P3A_36_group_executed"])
        self.assertEqual(r["scientific_status"],"EXPLORATORY_HOLD")

    def test_all_original_qtl_n_values_positive(self):
        r=raw()
        r["panels"][1][6]=0
        with self.assertRaisesRegex(ValueError,"hash or original"):
            m.reconcile(r)

    def test_invalid_prior_grid_rejected(self):
        r=raw()
        r["p12_values"][0]=0.0
        with self.assertRaisesRegex(ValueError,"prior/N"):
            m.reconcile(r)

if __name__ == "__main__":
    unittest.main()
