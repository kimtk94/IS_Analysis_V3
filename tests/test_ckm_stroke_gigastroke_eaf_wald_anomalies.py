"""Tests for original-source four-phenotype EAF/Wald anomaly summary."""
import unittest
from scripts.ckm_stroke.audit_gigastroke_eaf_gc import review_diagnostics
from scripts.ckm_stroke.review_gigastroke_eaf_wald_anomalies import analyze,to_markdown
from tests.test_ckm_stroke_gigastroke_eaf_gc import full

def source_fixture():
    src=full()
    for row in src["outcomes"]:
        e=row["eaf_gc_source_descriptive"]
        e["categories"]["hrs_absdiff_lt_002"]-=1
        e["categories"]["hrs_absdiff_ge_020"]=1
        ch=e["by_chr"][10]["categories"]
        ch["hrs_absdiff_lt_002"]-=1
        ch["hrs_absdiff_ge_020"]=1
        if row["phenotype"] in ("AIS","CE"):
            e["categories"]["abs_log10p_delta_lt_001"]-=14
            e["categories"]["abs_log10p_delta_010_to_1"]=14
            ch4=e["by_chr"][3]["categories"]
            ch4["abs_log10p_delta_lt_001"]-=14
            ch4["abs_log10p_delta_010_to_1"]=14
    return src

class TestAnomalyReview(unittest.TestCase):
    def test_chr11_is_shared_chromosome_not_shared_variant(self):
        src=source_fixture()
        r=analyze(src,review_diagnostics(src))
        self.assertEqual(r["extreme_HRS_reference_delta_shared_chromosomes"],[11])
        self.assertFalse(r["per_variant_shared_identity_across_subtypes_established"])
        self.assertEqual([x["P_vs_Wald_delta_ge_010"] for x in r["outcomes"]],[14,0,14,0])
        self.assertFalse(r["GC_correction_existence_or_factor_inferred"])
        self.assertFalse(r["ready_for_ldsc"])
        self.assertIn("same rsID",to_markdown(r))
    def test_review_mutation_rejected(self):
        src=source_fixture()
        rev=review_diagnostics(src)
        rev["outcomes"][0]["HRS_reference_freq_absdiff"]["hrs_absdiff_ge_020"]=0
        with self.assertRaisesRegex(ValueError,"disagrees"):
            analyze(src,rev)
    def test_full_publisher_md5_tamper_rejected(self):
        src=source_fixture();rev=review_diagnostics(src)
        src["outcomes"][0]["raw_md5_verified"]=False
        with self.assertRaisesRegex(ValueError,"MD5"):
            analyze(src,rev)
    def test_chr_count_tamper_rejected(self):
        src=source_fixture();rev=review_diagnostics(src)
        src["outcomes"][0]["eaf_gc_source_descriptive"]["by_chr"][10]["categories"]["hrs_absdiff_ge_020"]+=1
        with self.assertRaises(ValueError):
            analyze(src,rev)
    def test_no_extreme_hrs_yields_no_shared_chr(self):
        src=full()
        r=analyze(src,review_diagnostics(src))
        self.assertEqual(r["extreme_HRS_reference_delta_shared_chromosomes"],[])
        self.assertEqual(r["scientific_status"],"HOLD")
if __name__=="__main__":unittest.main()
