"""Pure offline tests for exact duplicate and zero-P original BMI reconciliation."""
import gzip
import hashlib
import io
import math
import unittest
from unittest.mock import patch

from scripts.ckm_stroke.audit_bmi_duplicate_zero_p import (
    analyze, minus_log10_two_sided_p_from_z, signed_beta_se_z, iter_original_rows,
)
from scripts.ckm_stroke.audit_bmi_full_stream import EXPECTED, SOURCE


def row(ch, pos, snp, p="0.2", beta="0.1", se="0.05", n="600000"):
    return "\t".join([str(ch),str(pos),snp,"A","C","0.33",beta,se,p,n])


class Response(io.BytesIO):
    def __init__(self, raw):
        super().__init__(raw)
        self.headers={"Content-Length": str(len(raw))}
        self.status=200

    def geturl(self):
        return SOURCE


def mock_open(rows):
    raw = gzip.compress(("\t".join(EXPECTED)+"\n" + "\n".join(rows)+"\n").encode())
    return lambda req, timeout: Response(raw), hashlib.sha256(raw).hexdigest()


class DuplicateAndZeroPTests(unittest.TestCase):
    def test_z_sign_and_normal_p(self):
        self.assertEqual(signed_beta_se_z("-0.4","0.1"),-4.0)
        self.assertAlmostEqual(minus_log10_two_sided_p_from_z(0),0)
        self.assertAlmostEqual(minus_log10_two_sided_p_from_z(3.0),
                               -math.log10(math.erfc(3.0/math.sqrt(2))),places=10)

    def test_large_z_avoids_numeric_underflow(self):
        large=minus_log10_two_sided_p_from_z(90)
        self.assertTrue(math.isfinite(large))
        self.assertGreater(large,1700)
        self.assertLess(large,1800)

    def test_reject_zero_or_negative_standard_error(self):
        with self.assertRaisesRegex(ValueError, "invalid"):
            signed_beta_se_z("1","0")

    def test_two_passes_identical_duplicate(self):
        op,sha=mock_open([row(1,100,"rs1"),row(1,100,"rs1"),row(2,300,"rs2",p="0",beta="1",se=".02")])
        with patch("scripts.ckm_stroke.audit_bmi_duplicate_zero_p.EXPECTED_BMI_SHA256",sha):
            r=analyze(op)
        self.assertEqual(r["counts"]["total_rows"],3)
        self.assertEqual(r["counts"]["duplicate_rsid_extra_rows"],1)
        self.assertEqual(r["counts"]["duplicate_coordinate_extra_rows"],1)
        self.assertEqual(r["counts"]["p_underflow_zero_rows"],1)
        self.assertEqual(r["counts"]["two_complete_sha256_matched_passes"],2)
        self.assertEqual(r["duplicate_groups"][0]["dedup_decision"],"SAFE_TO_DROP_IDENTICAL_REPEAT")
        self.assertTrue(r["dedup_ready_for_downstream"])
        self.assertFalse(r["ready_for_ldsc"])
        self.assertEqual(r["zero_p_records"][0]["published_P"],"0")
        self.assertFalse(r["zero_p_records"][0]["derived_P_is_original"])

    def test_discordant_stats_are_not_silently_deduplicated(self):
        op,sha=mock_open([row(1,100,"rs1"),row(1,100,"rs1",beta="0.2"),row(1,300,"rs2")])
        with patch("scripts.ckm_stroke.audit_bmi_duplicate_zero_p.EXPECTED_BMI_SHA256",sha):
            r=analyze(op)
        self.assertEqual(r["discordant_duplicate_rsids"],["rs1"])
        self.assertFalse(r["dedup_ready_for_downstream"])
        self.assertFalse(r["duplicate_groups"][0]["identical_beta_se_p_n"])

    def test_cross_id_coordinate_collision_holds(self):
        op,sha=mock_open([row(1,100,"rs1"),row(1,100,"rs2"),row(2,500,"rs3")])
        with patch("scripts.ckm_stroke.audit_bmi_duplicate_zero_p.EXPECTED_BMI_SHA256",sha):
            r=analyze(op)
        self.assertFalse(r["dedup_ready_for_downstream"])
        self.assertEqual(len(r["duplicate_different_rsids_at_same_coordinate"]),1)

    def test_original_hash_pin_fail_closed(self):
        op,_=mock_open([row(1,100,"rs1")])
        with self.assertRaisesRegex(ValueError,"SHA256"):
            list(iter_original_rows(op))

    def test_no_duplicate_still_no_ldsc_approval(self):
        op,sha=mock_open([row(1,100,"rs1"),row(2,500,"rs2")])
        with patch("scripts.ckm_stroke.audit_bmi_duplicate_zero_p.EXPECTED_BMI_SHA256",sha):
            r=analyze(op)
        self.assertEqual(r["counts"]["two_complete_sha256_matched_passes"],1)
        self.assertFalse(r["ready_for_ldsc"])
        self.assertEqual(r["scientific_status"],"HOLD")


if __name__=="__main__":
    unittest.main()
