import csv
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "is"))
from audit_gwas_qtl_snp_harmonization import audit, classify

class TestHarmonization(unittest.TestCase):
    def test_matched(self):
        self.assertEqual(classify(("A","C","C","A",1,0.01),("A","C","C","A",2,0.02)), ("MATCHED",1))
    def test_swapped(self):
        self.assertEqual(classify(("A","C","C","A",1,0.01),("A","C","A","C",2,0.02)), ("SWAPPED",-1))
    def test_palindromic_fail_closed(self):
        self.assertEqual(classify(("A","T","A","T",1,0.01),("A","T","A","T",2,0.02))[0], "PALINDROMIC_UNRESOLVED")
    def test_strand_fail_closed(self):
        self.assertEqual(classify(("A","C","A","C",1,0.01),("T","G","T","G",2,0.02))[0], "STRAND_UNRESOLVED")
    def test_build_mismatch(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError, "Genome builds"):
                audit("unused", "unused", "GRCh37", "GRCh38", d)
    def test_end_to_end(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            fields=["chrom","pos","ref","alt","effect_allele","other_allele","beta","pvalue"]
            for filename, allele in (("gwas.tsv","C"),("qtl.tsv","A")):
                with (p/filename).open("w",newline="") as f:
                    w=csv.DictWriter(f,fieldnames=fields,delimiter="\t");w.writeheader()
                    w.writerow(dict(chrom="1",pos="123",ref="A",alt="C",
                                    effect_allele=allele,other_allele="A" if allele=="C" else "C",
                                    beta="0.2",pvalue="0.01"))
            report=audit(p/"gwas.tsv",p/"qtl.tsv","GRCh37","GRCh37",p/"out")
            self.assertEqual(report["counts"]["SWAPPED"],1)
            self.assertFalse(report["ld_verified"])
if __name__=="__main__":unittest.main()
