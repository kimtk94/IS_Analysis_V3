"""Offline tests for publisher-MD5-pinned single outcome locus parsing."""
import gzip,hashlib,io,unittest
from unittest.mock import patch
from scripts.ckm_stroke import match_tpmi_43321_gigastroke_ais as m

def target():
    return {"schema":"ckm-TPMI43321-GIGASTROKE-GRCh37-target-request-v1",
      "technical_status":"PASS_VERIFIED_TARGET_ONLY",
      "target_build":"GRCh37","scientific_status":"HOLD",
      "GIGASTROKE_original_row_matched":False,
      "binary_case_control_Neff_approved":False,
      "target_chr":"11","target_pos":10000,"target_REF":"A",
      "target_ALT":"G","target_A1":"G","target_variant_key":"11:10000:A:G",
      "TPMI_A1_logOR":0.18,"TPMI_A1_logOR_SE":0.03,"TPMI_original_P":1e-8}
def mockgzip(effect="G",other="A",duplicate=False):
    hdr="\t".join(m.COLUMNS)+"\n"
    row="\t".join(["11","10000","0.3","0.12","0.05","0.016","1.1275","1.02","1.24",effect,other])+"\n"
    otherrow="\t".join(["2","1234","0.2","0.02","0.05","0.7","1.02","0.95","1.08","T","C"])+"\n"
    return gzip.compress((hdr+otherrow+row+(row if duplicate else "")).encode())
class TestAIS(unittest.TestCase):
    def run_source(self, b,t=None):
        with patch.object(m,"SIZE",len(b)),patch.object(m,"PUBLISHER_MD5",hashlib.md5(b).hexdigest()):
            return m.scan(io.BytesIO(b),target() if t is None else t)
    def test_exact_locus_effect_orientation(self):
        r=self.run_source(mockgzip())
        self.assertEqual(r["technical_status"],"PASS_EXACT_LOCUS_AND_ALLELE_MATCH")
        self.assertEqual(r["signed_beta_to_TPMI_A1"],0.12)
        self.assertFalse(r["independent_replication_claimed"])
        self.assertFalse(r["GIGASTROKE_Neff_approved"])
        self.assertTrue(r["observed_effect_direction_concordant"])
    def test_swapped_effect_sign(self):
        r=self.run_source(mockgzip("A","G"))
        self.assertEqual(r["signed_beta_to_TPMI_A1"],-0.12)
        self.assertFalse(r["observed_effect_direction_concordant"])
    def test_duplicate_locus_is_rejected(self):
        with self.assertRaisesRegex(ValueError,"duplicate"):
            self.run_source(mockgzip(duplicate=True))
    def test_publisher_md5_is_required(self):
        b=mockgzip()
        with patch.object(m,"SIZE",len(b)):
            with self.assertRaisesRegex(ValueError,"MD5"):
                m.scan(io.BytesIO(b),target())
    def test_missing_original_tpmi_effect_is_rejected(self):
        x=target();del x["TPMI_A1_logOR"]
        with self.assertRaises(ValueError):m.target_gate(x)
    def test_invalid_target_is_rejected(self):
        x=target();x["scientific_status"]="GO"
        with self.assertRaises(ValueError):m.target_gate(x)
    def test_missing_coordinate_is_hold(self):
        b=mockgzip()
        x=target();x["target_pos"]=12345;x["target_variant_key"]="11:12345:A:G"
        self.assertEqual(self.run_source(b,x)["technical_status"],"HOLD_NO_EXACT_GIGASTROKE_COORDINATE")
    def test_palindromic_orientation_remains_hold(self):
        x=target();x.update(target_REF="A",target_ALT="T",target_A1="T",target_variant_key="11:10000:A:T")
        r=self.run_source(mockgzip("T","A"),x)
        self.assertEqual(r["technical_status"],"HOLD_ALLELE_ORIENTATION")
if __name__=="__main__":unittest.main()
