import importlib.util
import unittest
from pathlib import Path

SCRIPT=Path(__file__).resolve().parents[1]/"scripts"/"prepare_ckd_stage2b_coloc.py"
SPEC=importlib.util.spec_from_file_location("stage2b",SCRIPT)
M=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)

class TestStage2B(unittest.TestCase):
    def test_chr(self):
        self.assertEqual(M.norm_chr("chr16"),"16")
        self.assertEqual(M.norm_chr("23"),"X")
    def test_complement(self):
        self.assertEqual(M.complement("A"),"T")
        self.assertEqual(M.complement("CG"),"GC")
    def test_snv(self):
        self.assertTrue(M.is_snv("A","G"))
        self.assertFalse(M.is_snv("A","AT"))
    def test_p(self):
        self.assertAlmostEqual(M.p_from_log10(2),0.01)
    def test_match_direct(self):
        p={"allele0_pqtl_hg37":"A","allele1_pqtl_hg37":"G"}
        o={"effect_allele_outcome":"G","other_allele_outcome":"A","p_outcome":0.1}
        m=M.match_outcome(p,[o])
        self.assertEqual(m[1],"direct"); self.assertEqual(m[2],1)
    def test_match_swapped(self):
        p={"allele0_pqtl_hg37":"A","allele1_pqtl_hg37":"G"}
        o={"effect_allele_outcome":"A","other_allele_outcome":"G","p_outcome":0.1}
        m=M.match_outcome(p,[o])
        self.assertEqual(m[1],"swapped"); self.assertEqual(m[2],-1)

if __name__=="__main__":
    unittest.main()
