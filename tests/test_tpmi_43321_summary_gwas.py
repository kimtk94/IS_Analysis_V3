import gzip
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

s=Path(__file__).resolve().parents[1]/"scripts/is/audit_tpmi_43321_summary_gwas.py"
spec=importlib.util.spec_from_file_location("tpmi_gwas_parser",s)
m=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=m
spec.loader.exec_module(m)

class TestTPMIGwasParser(unittest.TestCase):
    def setUp(self):
        self.targets=[
            {"gene":"ALDH2","rsid":"rs671","GRCh37":"12:112241766:G:A",
             "GRCh38":"12:111803962:G:A"},
            {"gene":"ADH1B","rsid":"rs1229984","GRCh37":"4:100239319:T:C",
             "GRCh38":"4:99318162:T:C"}]
        self.header=["chrom","pos","ref","alt","rsids","pval",
                     "beta","sebeta","af","case_af","control_af"]

    def run_gwas(self,records):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"433.21.tsv.gz"
            with gzip.open(p,"wt") as f:
                f.write("\t".join(self.header)+"\n")
                for r in records:
                    f.write("\t".join(str(r.get(x,"")) for x in self.header)+"\n")
            return m.extract(p,self.targets)

    def test_exact_ref_alt_and_missing(self):
        row={"chrom":"12","pos":111803962,"ref":"G","alt":"A",
             "rsids":"rs671","pval":"1e-6","beta":"-0.2",
             "sebeta":"0.04","af":"0.25"}
        records,n=self.run_gwas([row])
        self.assertEqual(n,1)
        self.assertEqual(records[0]["status"],"EXACT_GRCH38_REF_ALT_MATCH")
        self.assertAlmostEqual(records[0]["tpmi_beta_ALT"],-0.2)
        self.assertEqual(records[1]["status"],"NOT_REPORTED")

    def test_wrong_allele_pair_review(self):
        row={"chrom":"12","pos":111803962,"ref":"A","alt":"G",
             "rsids":"rs671","pval":"0.5","beta":"0.02",
             "sebeta":"0.04","af":"0.25"}
        rows,_=self.run_gwas([row])
        self.assertEqual(rows[0]["status"],"SAME_POSITION_DIFFERENT_ALLELE_PAIR_REVIEW")

    def test_invalid_probability(self):
        row={"chrom":"12","pos":111803962,"ref":"G","alt":"A",
             "rsids":"rs671","pval":"3","beta":"-0.2",
             "sebeta":"0.04","af":"0.25"}
        with self.assertRaises(ValueError):self.run_gwas([row])

    def test_duplicate_exact_is_error(self):
        row={"chrom":"12","pos":111803962,"ref":"G","alt":"A",
             "rsids":"rs671","pval":"1e-6","beta":"-0.2",
             "sebeta":"0.04","af":"0.25"}
        with self.assertRaises(ValueError):self.run_gwas([row,row])

if __name__=="__main__":unittest.main()
