import io
import math
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import publication_sprint_followup as f

class FollowupTests(unittest.TestCase):
    def test_bh_includes_all_terms(self):
        self.assertEqual(f.bh([.01,1,.02]),[.03,1,.03])
    def test_hypergeom_exact(self):
        self.assertAlmostEqual(f.hypergeom_tail(4,2,2,2),1/6)
        self.assertEqual(f.hypergeom_tail(4,2,2,0),1)
    def test_skin_scan(self):
        text="CHR\tBP\tEA\tNEA\tBETA\tSE\tP\n1\t123\tA\tC\t0.1\t0.01\t1e-9\n2\t456\tA\tG\t0.1\t0.01\t0\n1\t2\tA\tA\t0.1\t0.01\t1e-12\n"
        r=f.scan_gwas(io.StringIO(text))
        self.assertEqual(r["valid_autosomal_snv_rows"],2)
        self.assertEqual(r["positive_p_lt_5e8_rows"],1)
        self.assertEqual(r["zero_p_rows_not_counted_significant"],1)
        self.assertEqual(r["excluded_rows"],{"not_biallelic_snv":1})
    def test_missing_schema(self):
        with self.assertRaises(ValueError):
            f.skin_schema(["CHR","BP","P"])
    def test_compare_not_independent_model(self):
        a={"ENSG1":(2,.01,"A"),"ENSG2":(2,.01,"B")}
        b={"ENSG1":(1,.01,"A"),"ENSG2":(1,.01,"B")}
        c={"ENSG1":{"fdr":.01,"log2fc":3,"symbol":"A"},"ENSG2":{"fdr":.01,"log2fc":-3,"symbol":"B"}}
        n,rows=f.compare(a,b,c)
        self.assertEqual(n,2)
        self.assertTrue(rows[0]["strong_concordant"])
        self.assertTrue(rows[1]["sign_discordant"])

if __name__=="__main__":
    unittest.main()
