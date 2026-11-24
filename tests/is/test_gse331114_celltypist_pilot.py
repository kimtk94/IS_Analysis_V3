#!/usr/bin/env python3
"""Unit tests for reference-pilot confidence semantics and cohort provenance."""
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"is"))
from run_gse331114_celltypist_pilot import (
    select_gsm, summarize_predictions, MODEL_NAME, MIN_MODEL_OVERLAP_FRACTION,
)


class ReferencePilotUnitTests(unittest.TestCase):
    def test_fixed_broad_immune_model(self):
        self.assertEqual(MODEL_NAME,"Immune_All_High.pkl")
        self.assertGreaterEqual(MIN_MODEL_OVERLAP_FRACTION,0.6)

    def test_gsm_by_group_not_invented_donor(self):
        samples={
            "GSM0001":{"accession":"GSM0001","title":"AIS_brain1","characteristics":[]},
            "GSM0002":{"accession":"GSM0002","title":"AIS_brain2","characteristics":[]},
            "GSM0003":{"accession":"GSM0003","title":"AIS_brain3","characteristics":[]},
            "GSM0004":{"accession":"GSM0004","title":"AIS_brain4","characteristics":[]},
            "GSM0005":{"accession":"GSM0005","title":"Control_brain1","characteristics":[]},
            "GSM0006":{"accession":"GSM0006","title":"Control_brain2","characteristics":[]},
            "GSM0007":{"accession":"GSM0007","title":"Control_brain3","characteristics":[]},
            "GSM0008":{"accession":"GSM0008","title":"Control_brain4","characteristics":[]},
        }
        self.assertEqual(select_gsm("AIS",0,samples)["accession"],"GSM0001")
        self.assertEqual(select_gsm("Control",0,samples)["accession"],"GSM0005")

    def test_low_score_not_claimed_accuracy(self):
        labels=["T cell"]*60+["Monocyte"]*40
        scores=[0.9]*60+[0.3]*40
        output=summarize_predictions(labels,scores,100)
        self.assertEqual(output["n_predictions"],100)
        self.assertEqual(output["fraction_score_ge_threshold"],0.6)
        self.assertTrue(output["celltypist_raw_top_class_score_is_not_calibrated_probability"])
        self.assertEqual(output["provisional_labels_top_20"][0]["n_cells"],60)

    def test_predictions_missing_cells_raise(self):
        with self.assertRaisesRegex(ValueError,"Prediction rows"):
            summarize_predictions(["a"]*100,[0.5]*99,100)

    def test_invalid_score_raise(self):
        with self.assertRaisesRegex(ValueError,"outside"):
            summarize_predictions(["a"]*100,[1.5]*100,100)

    def test_low_n_not_a_valid_pilot(self):
        with self.assertRaisesRegex(ValueError,"Not enough"):
            summarize_predictions(["a"]*99,[0.5]*99,99)


if __name__=="__main__":
    unittest.main()
