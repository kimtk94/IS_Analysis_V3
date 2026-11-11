import unittest
from unittest.mock import patch
from scripts.ckm_stroke import glgc_region_probe as probe

HEADER="rsID\tCHROM\tPOS_b37\tREF\tALT\tN\tEFFECT_SIZE\tSE\tpvalue\n"
DATA=(HEADER+"rs1\t22\t20000010\tA\tG\t100000\t0.10\t0.03\t0.000858120666\n"+
           "rs2\t22\t20000020\tA\tAT\t100000\t0.10\t0.03\t0.000858120666\n")

class TestGLGCIndexedRegion(unittest.TestCase):
    def test_one_valid_snp_one_indel_region_pass(self):
        r=probe.validate_rows(DATA,"22",20000000,20050000)
        self.assertEqual(r["status"],"PASS")
        self.assertEqual(r["rows_region"],2)
        self.assertEqual(r["non_snp_rows"],1)
        self.assertEqual(r["sample_variants"][0]["alt"],"G")

    def test_wrong_region_fail_closed(self):
        r=probe.validate_rows(HEADER+"rs1\t1\t20000010\tA\tG\t1000\t0.1\t0.03\t0.001\n",
                               "22",20000000,20050000)
        self.assertEqual(r["status"],"HOLD")
        self.assertEqual(r["invalid_rows"],1)

    def test_unavailable_index_is_not_false_pass(self):
        class Completed:
            returncode=1
            stdout=""
            stderr="failed to load .tbi"
        with patch.object(probe.subprocess,"run",return_value=Completed()) as stub:
            x=probe.probe("TG")
        self.assertEqual(x["status"],"SOURCE_UNAVAILABLE")
        self.assertIn("index",x["advice"])
        self.assertEqual(stub.call_args.args[0][0],"tabix")

    def test_successful_probe_does_not_claim_full_gwas(self):
        class Completed:
            returncode=0
            stdout=DATA
            stderr=""
        with patch.object(probe.subprocess,"run",return_value=Completed()):
            x=probe.probe("LDL_C")
        self.assertEqual(x["status"],"PASS")
        self.assertEqual(x["scope"],"INDEXED_REGION_ONLY_NOT_FULL_GWAS")
        self.assertEqual(x["scientific_status"],"NOT_RUN")
        self.assertIn("LDL_INV_EUR",x["url"])

    def test_unknown_trait_rejected(self):
        with self.assertRaises(ValueError):
            probe.probe("unknown")

if __name__=="__main__":
    unittest.main()
