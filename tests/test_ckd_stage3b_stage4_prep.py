#!/usr/bin/env python3
import importlib.util, unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

s3=load("s3",ROOT/"scripts/refine_ckd_stage3b_specificity.py")
s4a=load("s4a",ROOT/"scripts/audit_ckd_stage4_koges_inputs.py")

class Tests(unittest.TestCase):
    def test_compartments(self):
        self.assertEqual(s3.compartment("proximal tubule cells"),"renal_epithelial")
        self.assertEqual(s3.compartment("macrophages"),"immune")
        self.assertEqual(s3.compartment("vascular endothelial cells"),"vascular")
    def test_koges_audit_helpers(self):
        self.assertFalse(s4a.is_palindromic("C","T"))
        self.assertTrue(s4a.is_palindromic("C","G"))
        anchors=[{
            "gene_symbol":"ACP1","rsid":"rs1","chrom_hg19":"2",
            "pos_hg19":"100","ref_allele":"C","effect_allele_alt":"T",
            "palindromic":0
        }]
        records=[{
            "source":"x.pvar","source_kind":"plink2_pvar","line":1,
            "chrom":"2","id":"rs1","pos":"100","ref":"C","alt":"T",
            "a1":"","a2":""
        }]
        out=s4a.match_anchors(records,anchors,False)
        self.assertEqual(out[0]["found"],1)
        self.assertEqual(out[0]["match_type"],"rsid")
        self.assertEqual(out[0]["allele_status"],"ref_alt_match")

    def test_summary_patterns(self):
        rows=[
          {"gene_symbol":"A","cell_type":"x proximal tubule","expression":"10"},
          {"gene_symbol":"A","cell_type":"macrophages","expression":"1"},
          {"gene_symbol":"B","cell_type":"loop of henle epithelial cells","expression":"0.2"},
          {"gene_symbol":"B","cell_type":"collecting duct cells","expression":"0.2"},
        ]
        out={x["gene_symbol"]:x for x in s3.summarize(rows)}
        self.assertEqual(out["A"]["localization_pattern"],"focused")
        self.assertEqual(out["A"]["top_compartment"],"renal_epithelial")
        self.assertEqual(out["B"]["localization_pattern"],"very_low_expression")

if __name__=="__main__": unittest.main()
