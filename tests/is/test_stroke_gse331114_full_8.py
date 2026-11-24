#!/usr/bin/env python3
import gzip
import sys
import tempfile
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"is"))
from audit_stroke_gse331114_full_8 import GENES, parse_features, stream_gene_counts

class FullStrokeGeneCountsTests(unittest.TestCase):
    def test_preset_is_fixed(self):
        self.assertEqual(len(GENES), 15)
        self.assertEqual(len(GENES),len(set(GENES)))
    def test_sparse_aggregate_and_cells(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d)
            feat=d/"features.gz"
            mat=d/"matrix.gz"
            with gzip.open(feat,"wt") as f:
                f.write("ENSG1\tABO\tGene Expression\nENSG2\tGRK5\tGene Expression\nENSG3\tMS4A1\tGene Expression\n")
            with gzip.open(mat,"wt") as f:
                f.write("%%MatrixMarket matrix coordinate integer general\n3 2 3\n1 1 5\n2 2 3\n3 2 2\n")
            n,symbols,duplicates=parse_features(feat)
            self.assertEqual((n,duplicates),(3,[]))
            x=stream_gene_counts(mat,n,symbols)
            self.assertEqual(x["cell_barcodes"],2)
            self.assertEqual(x["selected_gene_umis"]["ABO"],5.0)
            self.assertEqual(x["selected_genes_detected_cells"]["GRK5"],1)
            self.assertEqual(x["total_umi_count"],10.0)
    def test_matrix_truncation_is_not_success(self):
        with tempfile.TemporaryDirectory() as d:
            mat=Path(d)/"matrix.gz"
            with gzip.open(mat,"wt") as f:
                f.write("%%MatrixMarket matrix coordinate integer general\n3 2 4\n1 1 5\n")
            with self.assertRaises(ValueError):
                stream_gene_counts(mat,3,{})
    def test_duplicate_gene_symbol_is_review(self):
        with tempfile.TemporaryDirectory() as d:
            feat=Path(d)/"features.gz"
            with gzip.open(feat,"wt") as f:
                f.write("ENSG1\tABO\nENSG2\tABO\n")
            _,_,duplicates=parse_features(feat)
            self.assertEqual(duplicates,["ABO"])
if __name__=="__main__":
    unittest.main()
