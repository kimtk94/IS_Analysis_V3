import csv
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "is"))
from validate_public_nominal_qtl_readiness import validate, REQUIRED

class TestReadinessValidator(unittest.TestCase):
    def write(self, directory, **changes):
        row = dict.fromkeys(REQUIRED, "")
        row.update(gene_id_stable="ENSG00000000001", canonical_gene_symbol="TEST",
                   archived_646_coloc_assays="0", n_new_QTL_tissues_with_any_nominal_rows="0",
                   n_new_nominal_QTL_gene_variant_rows="0",
                   new_nominal_QTL_gene_source_status="NOT_ASSESSED",
                   GWAS_and_QTL_alleles_coharmonized="False",
                   matched_QTL_LD_verified="False", coloc_recomputed="False",
                   ischemic_stroke_causal_gene="NOT_ESTABLISHED")
        row.update(changes)
        p = Path(directory) / "fixture.tsv"
        with p.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=REQUIRED, delimiter="\t")
            writer.writeheader()
            writer.writerow(row)
        return p

    def test_valid_synthetic(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(validate(self.write(d), 1)["unique_genes"], 1)

    def test_unreviewed_causal_claim_fails(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError, "causal"):
                validate(self.write(d, ischemic_stroke_causal_gene="YES"), 1)

    def test_missing_pvalue_fails(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError, "p-value"):
                validate(self.write(d, n_new_nominal_QTL_gene_variant_rows="2"), 1)

    def test_bad_boolean_fails(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError, "boolean"):
                validate(self.write(d, matched_QTL_LD_verified="NA"), 1)

    def test_denominator_fails(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError, "Expected"):
                validate(self.write(d), 2225)

if __name__ == "__main__":
    unittest.main()
