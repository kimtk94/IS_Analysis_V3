"""Offline safety checks: no fake GIGASTROKE replication on mere coordinate mapping."""
import copy
import unittest
from scripts.ckm_stroke.prepare_tpmi_43321_gigastroke_lookup import prepare, EXPECTED_SHA

def fixture():
    return {
      "schema":"ckm-TPMI43321-lead-variant-Ensembl-GRCh38-GRCh37-v1",
      "technical_status":"PASS_BIDIRECTIONAL_ENSEMBL_BUILD_MAP_AND_REF_ALLELE",
      "original_source_SHA256_verified":EXPECTED_SHA,
      "liftover_reverse_verified":True,
      "scientific_status":"HOLD","ready_for_ldsc":False,
      "GIGASTROKE_EUR_source_effect_at_target_observed":False,
      "lead_GRCh38":{"original_OBS_CT":204099,"P_original":1e-8,"reported_rsid":"rs123",
                      "logOR_original_A1":0.18,"LOGOR_SE_original":0.03},
      "mapped_GRCh37":{
         "chrom_GRCh37":"11","pos_GRCh37":6590000,
         "REF_GRCh37":"A","ALT_GRCh37":"G","A1_effect_allele_GRCh37":"G",
         "mapped_GRCh37_reference_verified":True,
         "allele_source_GRCh38_reference_verified":True
      }
    }

class LookupTests(unittest.TestCase):
    def test_request_is_a_lookup_not_a_replication(self):
        r=prepare(fixture())
        self.assertEqual(r["target_variant_key"],"11:6590000:A:G")
        self.assertEqual(r["technical_status"],"PASS_VERIFIED_TARGET_ONLY")
        self.assertFalse(r["GIGASTROKE_original_row_matched"])
        self.assertFalse(r["binary_case_control_Neff_approved"])
        self.assertEqual(r["scientific_status"],"HOLD")
        self.assertEqual(r["TPMI_A1_logOR"],0.18)
        self.assertEqual(r["TPMI_source_SHA256"],EXPECTED_SHA)
    def test_reject_nonfinite_effect(self):
        x=fixture();x["lead_GRCh38"]["logOR_original_A1"]=float("nan")
        with self.assertRaises(ValueError):prepare(x)
    def test_reject_missing_logor(self):
        x=fixture();del x["lead_GRCh38"]["logOR_original_A1"]
        with self.assertRaises(ValueError):prepare(x)
    def test_reject_wrong_pinned_hash(self):
        x=fixture();x["original_source_SHA256_verified"]="0"*64
        with self.assertRaises(ValueError):prepare(x)
    def test_reject_unverified_reference(self):
        x=fixture();x["mapped_GRCh37"]["mapped_GRCh37_reference_verified"]=False
        with self.assertRaises(ValueError):prepare(x)
    def test_reject_missing_or_invalid_alleles(self):
        for alt in ("A","N",""):
            with self.subTest(alt=alt):
                x=fixture();x["mapped_GRCh37"]["ALT_GRCh37"]=alt
                with self.assertRaises(ValueError):prepare(x)
    def test_reject_invalid_sample_n(self):
        x=fixture();x["lead_GRCh38"]["original_OBS_CT"]=313909
        with self.assertRaises(ValueError):prepare(x)
    def test_reject_fake_go(self):
        x=fixture();x["scientific_status"]="GO"
        with self.assertRaises(ValueError):prepare(x)
    def test_reject_unsupported_gigastroke_claim(self):
        x=fixture();x["GIGASTROKE_EUR_source_effect_at_target_observed"]=True
        with self.assertRaises(ValueError):prepare(x)

if __name__=="__main__":
    unittest.main()
