"""Regression checks for allele orientation and direct vs absent variant."""
import importlib.util
import tempfile
import unittest
from pathlib import Path

MODULE = Path(__file__).resolve().parents[1] / "scripts/is/audit_alcohol_cs_ais_direct_overlap.py"
spec = importlib.util.spec_from_file_location("audit_alcohol_cs_ais_direct_overlap", MODULE)
tool = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tool)

class TestAltHarmonization(unittest.TestCase):
    def test_effect_alt(self):
        x=tool.harmonize_alt({"effect_allele":"A","other_allele":"G","beta":"0.2",
             "se":"0.1","eaf":"0.3"},"A","G")
        self.assertEqual(x["orientation"],"EFFECT_ALT")
        self.assertAlmostEqual(x["beta_alt"],0.2)
        self.assertAlmostEqual(x["eaf_alt"],0.3)

    def test_effect_ref_reversed(self):
        x=tool.harmonize_alt({"effect_allele":"G","other_allele":"A","beta":"0.2",
             "se":"0.1","eaf":"0.3"},"A","G")
        self.assertEqual(x["orientation"],"EFFECT_REF_FLIPPED")
        self.assertAlmostEqual(x["beta_alt"],-0.2)
        self.assertAlmostEqual(x["eaf_alt"],0.7)

    def test_invalid_effect_alleles(self):
        x=tool.harmonize_alt({"effect_allele":"C","other_allele":"A",
             "beta":"0.2","se":"0.1","eaf":"0.3"},"A","G")
        self.assertEqual(x["orientation"],"INVALID_EFFECT_ALLELES")

    def test_scan_exact_vs_other_alleles(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/"sample.tsv"
            path.write_text("dataset\tbuild\tchr\tpos\tref\talt\tvariant_id\teffect_allele\tother_allele\tbeta\tse\tp\n"
                 "BBJ\tGRCh37\t12\t10\tG\tA\t12:10:G:A\tA\tG\t-0.3\t0.1\t0.01\n"
                 "BBJ\tGRCh37\t12\t20\tG\tT\t12:20:G:T\tT\tG\t0.3\t0.1\t0.01\n")
            targets=[{"chrom":"12","pos":"10","ref":"G","alt":"A","variant_id":"12:10:G:A"},
                     {"chrom":"12","pos":"20","ref":"G","alt":"A","variant_id":"12:20:G:A"}]
            good,other,n=tool.scan(path,targets)
            self.assertEqual(n,2)
            self.assertEqual(len(good["12:10:G:A"]),1)
            self.assertEqual(len(good["12:20:G:A"]),0)
            self.assertEqual(len(other["12:20:G:A"]),1)

if __name__=="__main__":
    unittest.main()
