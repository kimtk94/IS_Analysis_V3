import csv
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "is"))
from build_public_nominal_qtl_region_queue import build

class TestISNominalRegionQueue(unittest.TestCase):
    def fixture(self, d, ready=True):
        path = Path(d) / "input.tsv"
        fields = ["gene_id_stable", "region_id", "new_nominal_QTL_gene_source_status",
                  "source_nominal_QTL_ids", "coloc_recomputed",
                  "GWAS_and_QTL_alleles_coharmonized"]
        row = dict(gene_id_stable="ENSG00000000001", region_id="IS_XDATA_G0001",
                   new_nominal_QTL_gene_source_status=(
                       "RAW_NOMINAL_QTL_RECORDS_AVAILABLE_NO_COLOC" if ready else "NOT_ASSESSED"),
                   source_nominal_QTL_ids="QTD000131" if ready else "",
                   coloc_recomputed="False", GWAS_and_QTL_alleles_coharmonized="False")
        with path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            w.writeheader()
            w.writerow(row)
        return path

    def test_synthetic_queue(self):
        with tempfile.TemporaryDirectory() as d:
            m = build(self.fixture(d), Path(d) / "out", expected_genes=1)
            self.assertEqual(m["qtl_ready_genes"], 1)
            self.assertEqual(m["qtl_ready_regions"], 1)
            self.assertEqual(m["status"], "EXPLORATORY_QUEUE_ONLY")

    def test_denominator_mismatch(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError, "Expected"):
                build(self.fixture(d), Path(d) / "out")

    def test_no_ready_genes(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError, "No nominal"):
                build(self.fixture(d, ready=False), Path(d) / "out", expected_genes=1)

if __name__ == "__main__":
    unittest.main()
