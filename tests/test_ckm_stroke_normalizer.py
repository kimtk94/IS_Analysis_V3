import gzip
import tempfile
import unittest
from pathlib import Path
from scripts.ckm_stroke import normalize_sumstats as norm

class TestCKMNormalizer(unittest.TestCase):
    def test_ckdgen_effect_allele_and_build37(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            raw=root/"source.tsv.gz"
            with gzip.open(raw,"wt") as fh:
                fh.write("Chr Pos_b37 Allele1 Allele2 Effect StdErr P-value\n")
                fh.write("1   100 A G 0.22 0.07 0.001\n")
            output=root/"out.tsv.gz"
            report=norm.normalize(raw,output,"ckdgen2019")
            self.assertEqual(report["status"],"PASS")
            with gzip.open(output,"rt") as fh:
                lines=fh.read()
            self.assertIn("CHR\tPOS\tEA\tNEA\tBETA\tSE\tP",lines)
            self.assertIn("1\t100\tA\tG\t0.22",lines)

    def test_glgc_log_tg_inverse_normal_effect_allele(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            raw=root/"glgc.tsv"
            raw.write_text("CHROM\tPOS_b37\tREF\tALT\tEFFECT_SIZE\tSE\tpvalue\n"
                           "2\t500\tT\tC\t-0.42\t0.06\t0.002\n")
            output=root/"out.tsv.gz"
            report=norm.normalize(raw,output,"glgc2021")
            self.assertEqual(report["status"],"PASS")
            with gzip.open(output,"rt") as fh:
                self.assertIn("2\t500\tC\tT\t-0.42",fh.read())

    def test_exclusions_are_visible_and_hold(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            raw=root/"glgc.tsv"
            raw.write_text("CHROM\tPOS_b37\tREF\tALT\tEFFECT_SIZE\tSE\tpvalue\n"
                           "2\t500\tT\tC\t-0.42\t0.06\t0.002\n"
                           "X\t700\tG\tA\t0.3\t0.1\t0.01\n"
                           "2\t800\tT\tG\t0.4\t0.1\t0\n")
            result=norm.normalize(raw,root/"out.tsv.gz","glgc2021")
            self.assertEqual(result["status"],"HOLD")
            self.assertEqual(result["excluded"]["non_autosomal"],1)
            self.assertEqual(result["excluded"]["p_underflow"],1)
            self.assertEqual(result["output_rows"],1)

    def test_bad_source_sha_rejected_without_output(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            raw=root/"glgc.tsv"
            raw.write_text("CHROM\tPOS_b37\tREF\tALT\tEFFECT_SIZE\tSE\tpvalue\n"
                           "2\t500\tT\tC\t-0.42\t0.06\t0.002\n")
            target=root/"out.tsv.gz"
            with self.assertRaisesRegex(ValueError,"SHA256 mismatch"):
                norm.normalize(raw,target,"glgc2021","not-the-expected-hash")
            self.assertFalse(target.exists())

if __name__=="__main__":
    unittest.main()
