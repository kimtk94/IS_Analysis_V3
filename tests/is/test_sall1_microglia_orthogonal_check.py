#!/usr/bin/env python3
"""Marker evidence contracts; all tests use synthetic matrices, not stroke data."""
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"is"))
from sall1_microglia_orthogonal_check import (
    check_binary_markers, MARKERS, MODEL_LABELS,
)


def fixture():
    labels=["Monocytes","Monocytes","Macrophages","Macrophages","T cells","Mono-mac"]
    genes={
        "SALL1":[1,0,1,0,0,1],
        "TMEM119":[1,1,0,0,0,1],
        "P2RY12":[1,1,1,0,0,0],
        "FCN1":[0,1,1,1,0,0],
        "CD14":[0,1,1,1,0,0],
        "S100A8":[0,0,0,1,0,0],
    }
    return genes,labels


class MarkerCrosscheckTests(unittest.TestCase):
    def test_required_panel(self):
        self.assertEqual(set(MARKERS),{"SALL1","TMEM119","P2RY12","FCN1","CD14","S100A8"})

    def test_cross_detect_not_call_celltype(self):
        genes,labels=fixture()
        result=check_binary_markers(genes,labels)
        self.assertFalse(result["celltypes_validated"])
        self.assertFalse(result["paper_SALL1_high_classifier_reproduced"])
        self.assertEqual(result["n_total_barcodes"],len(labels))
        rows={r["model_provisional_label"]:r for r in result["model_categories"]}
        self.assertEqual(rows["Monocytes"]["n_barcodes"],2)
        self.assertEqual(rows["Monocytes"]["counts"]["SALL1_TMEM119_P2RY12_triple"],1)
        self.assertEqual(rows["Monocytes"]["counts"]["at_least_two_FCN1_CD14_S100A8"],1)
        self.assertEqual(rows["Macrophages"]["counts"]["both_SALL1_canonical_and_peripheral"],1)
        self.assertEqual(rows["ALL_PREDICTED_MY_MARKERS"]["n_barcodes"],5)
        self.assertEqual(rows["ALL_PREDICTED_MY_MARKERS"]["counts"]["SALL1_any"],3)
        self.assertEqual(rows["ALL_PREDICTED_MY_MARKERS"]["rates"]["SALL1_any"],0.6)

    def test_no_false_negative_identity(self):
        genes,labels=fixture()
        for gene in MARKERS:
            genes[gene]=[0]*len(labels)
        result=check_binary_markers(genes,labels)
        self.assertEqual(next(r for r in result["model_categories"]
                              if r["model_provisional_label"]=="ALL_PREDICTED_MY_MARKERS")
                         ["counts"]["SALL1_any"],0)
        self.assertFalse(result["myeloid_resident_versus_infiltrating_resolved"])

    def test_missing_marker_fails(self):
        genes,labels=fixture();genes.pop("SALL1")
        with self.assertRaisesRegex(ValueError,"all six"):
            check_binary_markers(genes,labels)

    def test_misaligned_marker_fails(self):
        genes,labels=fixture();genes["SALL1"]=genes["SALL1"][:-1]
        with self.assertRaisesRegex(ValueError,"length"):
            check_binary_markers(genes,labels)

    def test_invalid_binary_fails(self):
        genes,labels=fixture();genes["SALL1"][0]=2
        with self.assertRaisesRegex(ValueError,"not binary"):
            check_binary_markers(genes,labels)

    def test_zero_predicted_labels_is_a_defined_nonresult(self):
        genes,labels=fixture()
        result=check_binary_markers(genes,labels)
        row=next(r for r in result["model_categories"] if r["model_provisional_label"]=="DC")
        self.assertEqual(row["status"],"NO_PREDICTED_CELLS")
        self.assertIsNone(row["counts"])


if __name__=="__main__":
    unittest.main()
