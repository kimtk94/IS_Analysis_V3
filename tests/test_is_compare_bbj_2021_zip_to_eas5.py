import csv,gzip,importlib.util,io,tempfile,unittest,zipfile
from pathlib import Path

SCRIPT=Path(__file__).resolve().parents[1]/"scripts/is/compare_bbj_2021_zip_to_eas5.py"
spec=importlib.util.spec_from_file_location("bbjzip",SCRIPT)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class SourceGWASComparisonTests(unittest.TestCase):
    def test_matches_both_allele_orders_but_keeps_source_unapproved(self):
        with tempfile.TemporaryDirectory() as d:
            base=Path(d)
            archive=base/"eas5.tsv.gz"
            official=base/"bbj.zip"
            with gzip.open(archive,"wt",newline="") as f:
                cols=["mapped_genomic_variant_GRCh37","BBJ_beta_ALT","BBJ_p","BBJ_allele_exact_Fasta37"]
                w=csv.DictWriter(f,fieldnames=cols,delimiter="\t");w.writeheader()
                w.writerow(dict(zip(cols,["1:100:A:G","0.2","0.001","True"])))
                w.writerow(dict(zip(cols,["1:101:C:T","-0.3","0.02","True"])))
            buf=io.StringIO()
            cols=["CHR","POS","Allele1","Allele2","BETA","SE","p.value"]
            w=csv.DictWriter(buf,fieldnames=cols,delimiter="\t");w.writeheader()
            w.writerow({"CHR":"1","POS":"100","Allele1":"A","Allele2":"G","BETA":"0.2","SE":"0.03","p.value":"0.001"})
            w.writerow({"CHR":"1","POS":"101","Allele1":"T","Allele2":"C","BETA":"0.3","SE":"0.04","p.value":"0.02"})
            with zipfile.ZipFile(official,"w") as z:
                z.writestr("study.auto.txt.gz",gzip.compress(buf.getvalue().encode()))
            result=m.audit(archive,official)
            self.assertEqual(result["summary"]["same_coordinate_allele_pair"],2)
            self.assertEqual(result["summary"]["p_and_beta_match_either_sign"],2)
            self.assertEqual(result["summary"]["hypothesis_effect_A2_agree"],2)
            self.assertFalse(result["gwas_source_accepted"])
            self.assertFalse(result["sample_size_adopted"])
            self.assertFalse(result["coloc_executed"])
    def test_invalid_variant_and_incomplete_gwas(self):
        self.assertEqual(m.parsed_variant("chr1:123:a:g"),("1",123,"A","G"))
        self.assertIsNone(m.parsed_variant("rs123"))
        self.assertFalse(m.close_p(1e-10,1e-3))
if __name__=="__main__":unittest.main()
