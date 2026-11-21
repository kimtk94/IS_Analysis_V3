import gzip
import tempfile
import unittest
from pathlib import Path

from scripts.ckm_stroke import variant_qc

HEADER="CHR\tPOS\tEA\tNEA\tBETA\tSE\tP\n"

def write_data(path, lines):
    with gzip.open(path,"wt",encoding="utf-8") as out:
        out.write(HEADER+"".join(lines))


class TestVariantQC(unittest.TestCase):
    def test_unique_autosomal_variant_qc_pass(self):
        with tempfile.TemporaryDirectory() as d:
            target=Path(d)/"gwas.tsv.gz"
            write_data(target,[
                "1\t123\tA\tG\t0.1\t0.1\t0.31731050786291415\n",
                "2\t456\tT\tC\t-0.1\t0.1\t0.31731050786291415\n",
            ])
            v=variant_qc.variant_qc(target,work_dir=d,stride=1)
            self.assertEqual(v["status"],"PASS")
            self.assertEqual(v["rows"],2)
            self.assertEqual(v["distinct_variant_keys"],2)
            self.assertEqual(v["p_consistency_mismatches"],0)

    def test_swapped_allele_variants_are_flagged_not_silently_collapsed(self):
        with tempfile.TemporaryDirectory() as d:
            target=Path(d)/"gwas.tsv.gz"
            write_data(target,[
                "1\t123\tA\tT\t0.1\t0.1\t0.31731050786291415\n",
                "1\t123\tT\tA\t-0.1\t0.1\t0.31731050786291415\n",
                "1\t124\tG\tC\t0.1\t0.1\t0.31731050786291415\n",
            ])
            v=variant_qc.variant_qc(target,work_dir=d,stride=1)
            self.assertEqual(v["status"],"HOLD")
            self.assertEqual(v["duplicate_excess_rows"],1)
            self.assertEqual(v["duplicate_key_groups"],1)
            self.assertEqual(v["duplicate_opposing_effect_allele_groups"],1)
            self.assertEqual(v["palindromic_rows"],3)

    def test_discordant_reported_p_is_hold(self):
        with tempfile.TemporaryDirectory() as d:
            target=Path(d)/"gwas.tsv.gz"
            write_data(target,["1\t123\tA\tG\t0.1\t0.1\t1e-15\n"])
            v=variant_qc.variant_qc(target,work_dir=d,stride=1)
            self.assertEqual(v["status"],"HOLD")
            self.assertEqual(v["p_consistency_mismatches"],1)

    def test_bad_allele_is_fail(self):
        with tempfile.TemporaryDirectory() as d:
            target=Path(d)/"gwas.tsv.gz"
            write_data(target,["1\t123\tA\tA\t0.1\t0.1\t0.3\n"])
            with self.assertRaisesRegex(ValueError,"invalid standardized"):
                variant_qc.variant_qc(target,work_dir=d,stride=1)

if __name__=="__main__":
    unittest.main()
