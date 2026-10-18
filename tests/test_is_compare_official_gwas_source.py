import importlib.util
from pathlib import Path
import csv
import tempfile
import unittest

SCRIPT=Path(__file__).resolve().parents[1]/"scripts/is/compare_official_gwas_source.py"
spec=importlib.util.spec_from_file_location("comparison",SCRIPT)
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

def table(p,fields,rows):
    with p.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t")
        w.writeheader()
        w.writerows(rows)

class MatchingTest(unittest.TestCase):
    def test_source_matching_and_sign_flip_are_separate(self):
        with tempfile.TemporaryDirectory() as td:
            a=Path(td)/"archive.tsv"
            o=Path(td)/"official.tsv"
            ar=[]
            off=[]
            for i in range(1,26):
                ref,alt=("A","C") if i%2 else ("G","T")
                v=f"chr1:{i}:{ref}:{alt}"
                ar.append({"mapped_genomic_variant_GRCh37":v,
                           "BBJ_beta_ALT":"0.2","BBJ_p":"0.004",
                           "EAS_AIS_beta_ALT":"-0.3","EAS_AIS_p":"0.003"})
                off.append({"variant_GRCh37":v,"beta_ALT":"0.2" if i!=25 else "-0.2",
                            "p_value":"0.004"})
            table(a,["mapped_genomic_variant_GRCh37","BBJ_beta_ALT","BBJ_p",
                     "EAS_AIS_beta_ALT","EAS_AIS_p"],ar)
            table(o,["variant_GRCh37","beta_ALT","p_value"],off)
            res=module.match(a,o,"BBJ")
            self.assertEqual(res["counts"]["same_ALT_beta_and_p"],24)
            self.assertEqual(res["counts"]["opposite_beta_orientation"],1)
            self.assertEqual(res["source_statistics_screen"],"PASS_REVIEW_REQUIRED")
            self.assertFalse(res["approved_source"])
            self.assertFalse(res["sample_size_adopted"])
            self.assertFalse(res["coloc_executed"])

    def test_wrong_source_and_missing_columns_fail(self):
        with tempfile.TemporaryDirectory() as td:
            a=Path(td)/"a.tsv"
            o=Path(td)/"o.tsv"
            table(a,["mapped_genomic_variant_GRCh37","BBJ_beta_ALT","BBJ_p"],[])
            table(o,["random"],[])
            with self.assertRaises(ValueError):
                module.match(a,o,"NOT_BBJ")
            with self.assertRaises(ValueError):
                module.match(a,o,"BBJ")

    def test_variant_string_and_p_mismatch(self):
        self.assertEqual(module.variant("chr1:123:a:g"),"1:123:A:G")
        self.assertIsNone(module.variant("rs1"))
        self.assertIsNone(module.variant("1:123:A:A"))
        self.assertFalse(module.same_p(0.0001,0.1))
        self.assertFalse(module.same_p(0,1e-50))

if __name__=="__main__":unittest.main()
