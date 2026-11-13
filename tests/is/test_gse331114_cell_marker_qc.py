#!/usr/bin/env python3
"""Synthetic marker QC tests only. No claims of validated cell identity."""
from array import array
import gzip
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"is"))
from audit_gse331114_cell_marker_qc import (  # noqa: E402
    PANELS, feature_map, stream_cell_qc, classify_marker_support, aggregate_cell_qc,
    QC_MIN_GENES, QC_MAX_MITO_PCT,
)


class MarkerReadinessTest(unittest.TestCase):
    def test_marker_panels_unique(self):
        allgenes=[x for genes in PANELS.values() for x in genes]
        self.assertEqual(len(allgenes),len(set(allgenes)))
        self.assertTrue({"CD3D","TRAC"}.issubset(allgenes))
        self.assertTrue({"P2RY12","TMEM119"}.issubset(allgenes))

    def test_coexpression_requires_two_distinct_markers(self):
        marks={"CD3D":bytearray([1,1,1]),
               "CD3E":bytearray([1,0,1]),
               "TRAC":bytearray([0,0,1]),
               "LYZ":bytearray([0,0,1]),
               "LST1":bytearray([0,0,1])}
        self.assertEqual(classify_marker_support(0, marks),"T_cell_marker_support")
        self.assertEqual(classify_marker_support(1, marks),"UNASSIGNED")
        self.assertEqual(classify_marker_support(2, marks),"MULTIPLE_PANELS")

    def test_matrix_market_cell_count_qc(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"matrix.gz"
            with gzip.open(p,"wt") as f:
                f.write("%%MatrixMarket matrix coordinate integer general\n%\n4 2 5\n")
                f.write("1 1 5\n2 1 3\n3 1 1\n1 2 7\n4 2 2\n")
            x=stream_cell_qc(p,4,{3},{1:"CD3D",2:"CD3E"})
            self.assertEqual(x["cells"],2)
            self.assertEqual(x["total_umi"],18)
            self.assertEqual(list(x["detected_features_by_cell"]),[3,2])
            self.assertEqual(x["mitochondrial_umi_by_cell"][0],1)
            self.assertEqual(classify_marker_support(0,x["marker_detected_by_cell"]),"T_cell_marker_support")

    def test_fail_closed_if_header_rows_mismatch(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"bad.gz"
            with gzip.open(p,"wt") as f:
                f.write("%%MatrixMarket matrix coordinate integer general\n2 1 3\n1 1 5\n")
            with self.assertRaisesRegex(ValueError,"truncated"):
                stream_cell_qc(p,2,{1},{1:"CD3D"})

    def test_counts_filtered_and_labels_dont_imply_true_celltypes(self):
        marks={"CD3D":bytearray([1,1,1]),
               "CD3E":bytearray([1,0,1]),
               "TRAC":bytearray([0,0,1]),
               "LYZ":bytearray([0,0,1]),
               "LST1":bytearray([0,0,1])}
        result=aggregate_cell_qc({
            "cells":3,
            "total_umi_by_cell":array("Q",[1000,400,1000]),
            "detected_features_by_cell":array("I",[300,100,450]),
            "mitochondrial_umi_by_cell":array("Q",[50,20,300]),
            "marker_detected_by_cell":marks,
        },3)
        self.assertEqual(result["cells_pass_qc"],1)
        self.assertEqual(result["marker_panel_support_unvalidated_counts_after_QC"]["T_cell_marker_support"],1)
        self.assertTrue(result["requires_reference_validation_before_pseudobulk"])
        self.assertTrue(result["no_p_values_or_degs"])
        self.assertLess(QC_MAX_MITO_PCT,30)

    def test_feature_map_reports_mito_gene_and_marker(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"features.gz"
            markers=["CD3D","CD3E","TRAC","P2RY12","TMEM119","SALL1"]
            with gzip.open(p,"wt") as f:
                for i in range(36601):
                    g=(f"MT-{i}" if i<12 else
                       markers[i-12] if 12<=i<12+len(markers) else f"GENE{i}")
                    f.write(f"ENSG{i:011d}\t{g}\tGene Expression\n")
            n,mt,found,counts=feature_map(p)
            self.assertEqual(n,36601)
            self.assertEqual(len(mt),12)
            self.assertEqual(set(found.values()),set(markers))

    def test_bad_barcode_number_rejected(self):
        with self.assertRaisesRegex(ValueError,"Barcode count"):
            aggregate_cell_qc({
                "cells":2,
                "total_umi_by_cell":array("Q",[10,20]),
                "detected_features_by_cell":array("I",[300,300]),
                "mitochondrial_umi_by_cell":array("Q",[1,1]),
                "marker_detected_by_cell":{},
            },3)


if __name__=="__main__":
    unittest.main()
