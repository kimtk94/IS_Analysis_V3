import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from ocular_frozen_panel_archives import compare

class ArchiveTest(unittest.TestCase):
    def test_missing_is_not_negative(self):
        panel={"a":{"acute_log2fc":2,"ecm_member":True},"b":{"acute_log2fc":-2,"ecm_member":True}}
        models={"s":{"a":(1,.01),"b":(-2,.01)},"b":{"a":(1,.01),"b":(-1,.01)},"f":{"a":(.5,.01)}}
        chronic={"a":{"chronic_log2ratio":.25,"n_strains_same_sign":3,"source_subset":"P9"}}
        rows,r=compare(panel,models,chronic)
        self.assertEqual(r["all_241"]["all_three_alternate_models_tested"],1)
        self.assertEqual(r["all_241"]["all_alternate_signed_fdr05"],1)
        self.assertEqual(r["all_241"]["all_alternate_signed_fdr05_abs1"],0)
        self.assertEqual(r["all_241"]["chronic_covered"],1)
        self.assertEqual(r["all_241"]["both_annotation_and_chronic_3of3"],1)
        self.assertIsNone(rows[1]["f_fdr"])
        self.assertIsNone(rows[1]["n_strains_same_sign"])
if __name__=="__main__":unittest.main()
