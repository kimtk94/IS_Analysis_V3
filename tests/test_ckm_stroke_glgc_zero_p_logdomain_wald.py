import unittest
from scripts.ckm_stroke.glgc_zero_p_logdomain_wald import compare,normal_two_sided_neglog10p


def source(z=40,logp=None,p="0",pos=12345):
    if logp is None:
        logp=normal_two_sided_neglog10p(z)
    return {
        "CHROM":"19","POS_b37":str(pos),"rsID":"rs1",
        "REF":"A","ALT":"G",
        "N":"1000","POOLED_ALT_AF":"0.3",
        "EFFECT_SIZE":str(z),"SE":"1",
        "pvalue":p,"pvalue_neg_log10":str(logp)
    }


class TestLogDomainZeroP(unittest.TestCase):
    def test_large_z_pzero_remains_finite_and_agrees(self):
        report,rows=compare([source(z=50)],expected_count=1)
        self.assertEqual(report["status"],"PASS")
        self.assertGreater(rows[0]["source_negative_log10_p"],500)
        self.assertAlmostEqual(rows[0]["delta_log10_source_minus_wald"],0,places=8)
        self.assertEqual(report["scientific_status"],"HOLD")

    def test_sign_of_effect_does_not_change_two_sided_p(self):
        p1=normal_two_sided_neglog10p(60)
        p2=normal_two_sided_neglog10p(-60)
        self.assertAlmostEqual(p1,p2,places=10)

    def test_differing_p_values_reported_not_misleading_go(self):
        report,rows=compare([source(z=40,logp=400)],expected_count=1)
        self.assertEqual(report["status"],"PASS")
        self.assertEqual(report["abs_delta_gt_5"],1)
        self.assertEqual(report["scientific_status"],"HOLD")

    def test_duplicate_canonical_key_blocks(self):
        row1=source()
        row2=source()
        row2["REF"]="G";row2["ALT"]="A"
        report,rows=compare([row1,row2],expected_count=2)
        self.assertEqual(report["status"],"HOLD")
        self.assertTrue(report["errors"])

    def test_mismatched_chromosome_fails_closed(self):
        r=source()
        report,rows=compare([r],expected_count=1,expected_chromosome="8")
        self.assertEqual(report["status"],"HOLD")
        self.assertEqual(len(rows),0)

    def test_nonzero_source_p_fails_closed(self):
        report,rows=compare([source(p="0.05")],expected_count=1)
        self.assertEqual(report["status"],"HOLD")
        self.assertEqual(len(rows),0)

    def test_incomplete_expected_source_count_blocks(self):
        report,rows=compare([source()],expected_count=332)
        self.assertEqual(report["status"],"HOLD")

if __name__=="__main__":
    unittest.main()
