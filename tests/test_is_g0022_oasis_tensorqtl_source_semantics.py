"""Offline fail-closed regressions for original DDBJ Bokeh data lineage."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts/is"))
from audit_g0022_oasis_tensorqtl_source_semantics import (
    compare_source_browser, source_rows_from_prefix
)

def make_pair(n=600):
    raw=[];browser=[]
    for i in range(n):
        p=0.8 if i%2 else 0.06
        af=0.8 if i%3 else 0.2
        k=f"chr1_{55326+i}_T_C"
        raw.append({"phenotype_id":"ENSG00000237491","gene":"LINC01409",
                    "variant_id":k,"pval_nominal":str(p),
                    "slope":"-0.125","slope_se":"0.086","af":str(af)})
        browser.append({"variant_id_hg38":k.replace("_",":"),
                        "pval_nominal":p,"effect_size":-0.125,
                        "effect_size_SE":0.086,"maf":af})
    return raw,browser

class TestOriginalBokehLineage(unittest.TestCase):
    def test_600_exact_rows_with_high_af(self):
        a,b=make_pair()
        rows=compare_source_browser(a,b)
        self.assertEqual(len(rows),600)
        self.assertEqual(sum(float(r["source_af"])>.5 for r in rows),400)

    def test_precision_rounding_allowed(self):
        a,b=make_pair()
        b[0]["pval_nominal"]+=0.00049
        b[0]["effect_size"]+=0.00049
        b[0]["effect_size_SE"]+=0.00049
        b[0]["maf"]+=0.00049
        self.assertEqual(len(compare_source_browser(a,b)),600)

    def test_beta_sign_error_rejected(self):
        a,b=make_pair()
        b[200]["effect_size"]=+0.125
        with self.assertRaisesRegex(ValueError,"mismatch"):
            compare_source_browser(a,b)

    def test_allele_frequency_mislabel_rejected_when_values_differ(self):
        a,b=make_pair()
        b[4]["maf"]=0.2
        with self.assertRaisesRegex(ValueError,"mismatch"):
            compare_source_browser(a,b)

    def test_missing_variant_not_promoted_as_empty(self):
        a,b=make_pair()
        with self.assertRaisesRegex(ValueError,"absent"):
            compare_source_browser(a,b[:-1])

    def test_wrong_original_zip_rejected(self):
        with self.assertRaisesRegex(ValueError,"bounded ZIP"):
            source_rows_from_prefix(b"hello")

if __name__=="__main__":
    unittest.main()
