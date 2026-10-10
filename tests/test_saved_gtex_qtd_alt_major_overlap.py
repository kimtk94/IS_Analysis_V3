"""Protect ALT-major identity across saved eQTL Catalogue QTD and coloc input."""
import csv,gzip,importlib.util,sys,tempfile,unittest
from pathlib import Path
P=Path(__file__).resolve().parents[1]/"scripts/is/audit_saved_gtex_qtd_alt_major_overlap.py"
spec=importlib.util.spec_from_file_location("qtd_saved",P)
m=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=m
spec.loader.exec_module(m)

class TestOriginalSavedQTD(unittest.TestCase):
    def row(self,ac="640",an="774",ref="C",alt="A"):
        row=["ENSG00000089009","12","111652218",ref,alt,
          "chr12_111652218_"+ref+"_"+alt,"15","0.173127",
          "0.001","0.02","0.002","SNP",ac,an,"NA",
          "ENSG00000089009","ENSG00000089009","12.5","rs123"]
        return row
    def test_alt_major_ac_an_with_minor_maf(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"qtd.gz"
            with gzip.open(p,"wt") as f:
                f.write("\t".join(self.row())+"\n")
            r=list(m.raw_qtd_iter(p))
            self.assertEqual(len(r),1)
            self.assertTrue(r[0]["ALT_major_GTEx_ac_over_an"])
            self.assertAlmostEqual(r[0]["GTEx_ALT_ac_an"],640/774)
            self.assertTrue(r[0]["is_snv"])
    def test_indel_is_not_snp_and_kept_separately(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"qtd.gz"
            with gzip.open(p,"wt") as f:
                f.write("\t".join(self.row(ref="CT",alt="C"))+"\n")
            r=list(m.raw_qtd_iter(p))
            self.assertFalse(r[0]["is_snv"])
    def test_impossible_raw_allele_count_fails(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"qtd.gz"
            with gzip.open(p,"wt") as f:
                f.write("\t".join(self.row(ac="775"))+"\n")
            with self.assertRaises(ValueError):
                list(m.raw_qtd_iter(p))
if __name__=="__main__":
    unittest.main()
