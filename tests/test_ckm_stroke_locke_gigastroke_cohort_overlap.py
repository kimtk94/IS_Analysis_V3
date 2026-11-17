"""Offline strict-source cohort label/family and identity overlap safety tests."""
from __future__ import annotations
from collections import Counter
from unittest.mock import patch
import unittest

from scripts.ckm_stroke.audit_locke_gigastroke_cohort_overlap import (
    ALIASES,SOURCE,assert_publisher_file,compare,extract_gig,extract_locke,audit,
    EXPECTED_EUR,
)


def sample_sheets():
    locke={
        2:{"A":"Study","B":"Full name",
           "H":"Samples in analyses (N) (Trait - Cases/Controls)"},
        4:{"A":"Previous genome-wide association studies"},
    }
    # Fills 110 independent named study rows so original Table 1 bounds can be tested.
    for idx in range(5,125):
        locke[idx]={"A":f"StudyMock{idx}","B":f"Study title {idx}",
                    "H":"1000","D":"white European"}
    core=[("ARIC",8108),("CHS",3228),("FRAM",8094),("AGES",3207)]
    labels=[a[0] for a in ALIASES.values()]
    for i,name in enumerate(dict.fromkeys(labels)):
        idx=5+i
        expected=dict(core).get(name,1001)
        locke[idx]={"A":name,"B":name,"H":str(expected),"D":"white European"}
    gig={4:{"K":"AIS","L":"CES","M":"SVS","N":"LAS","P":"Controls"},
         49:{"A":"TOTAL EUR","K":"62100","L":"10804","M":"6811","N":"6399","P":"1234808"}}
    for idx,name in enumerate([*ALIASES.keys(),*["OtherMock"+str(i) for i in range(27)],"UKBiobank"],start=5):
        gig[idx]={"B":name,"K":"10","L":"NA","M":"NA","N":"NA","P":"100"}
    return locke,gig


class AliasTests(unittest.TestCase):
    def test_firstparty_cohort_label_and_family_tiers(self):
        original,stroke=sample_sheets()
        a=extract_locke(original);b=extract_gig(stroke)
        r=compare(a,b)
        self.assertEqual(r["technical_status"],
                         "PASS_FIRSTPARTY_TWO_SUPPLEMENT_COHORT_NAME_RECONCILIATION")
        self.assertEqual(r["GIGASTROKE_EUR_cohorts"],44)
        self.assertEqual(r["direct_cohort_label_matches"],10)
        self.assertEqual(r["study_family_only_matches"],6)
        self.assertEqual(r["unmatched_GIGASTROKE_study_labels"],28)
        self.assertIn("UKBiobank",r["unmatched_GIGASTROKE_labels"])
        self.assertTrue(all(v["same_individuals_verified"] is False for v in r["matched_study_rows"]))
        self.assertTrue(all(v["individual_overlap_count"] is None for v in r["matched_study_rows"]))
        self.assertFalse(r["ready_for_ldsc"])

    def test_expected_source_SHA_not_fake_approved(self):
        with self.assertRaisesRegex(ValueError,"SHA256"):
            assert_publisher_file("LOCKE",b"not the XLSX")
        with self.assertRaisesRegex(ValueError,"SHA256"):
            assert_publisher_file("GIGASTROKE",b"not the XLSX")

    def test_allele_like_duplicate_cohort_name_preserves_multiple_rows(self):
        original,_=sample_sheets()
        original[125]={"A":"ARIC","B":"ARIC Second Batch","H":"200","D":"white European"}
        data=extract_locke(original)
        self.assertEqual(data["ARIC"]["reported_analysis_n_leading_numeric"],8108)
        self.assertEqual(len(data["ARIC"]["additional_source_rows"]),1)
        self.assertEqual(data["ARIC"]["additional_source_rows"][0]["row"],125)

    def test_GIGASTROKE_original_CE_LAA_swap_rejected(self):
        _,stroke=sample_sheets()
        stroke[4]["L"]="LAS";stroke[4]["N"]="CES"
        with self.assertRaisesRegex(ValueError,"phenotype column order"):
            extract_gig(stroke)

    def test_GIGASTROKE_total_drift_rejected(self):
        _,stroke=sample_sheets()
        stroke[49]["N"]="6400"
        with self.assertRaisesRegex(ValueError,"totals changed"):
            extract_gig(stroke)

    def test_locke_anchor_ancestry_discordance_rejected(self):
        original,_=sample_sheets()
        for row in original.values():
            if row.get("A")=="ARIC":row["D"]="African American"
        with self.assertRaisesRegex(ValueError,"European ancestry"):
            extract_locke(original)

    def test_stale_alias_rejected_not_silently_reclassified(self):
        original,stroke=sample_sheets()
        a=extract_locke(original);b=extract_gig(stroke)
        b.pop("ARIC_EA*")
        with self.assertRaisesRegex(ValueError,"unsupported or stale alias"):
            compare(a,b)

    def test_unavailable_publisher_fail_closed(self):
        def broken(*args,**kwargs):
            raise ConnectionError("official source unavailable")
        report=audit(opener=broken)
        self.assertEqual(report["technical_status"],"HOLD_FIRSTPARTY_COHORT_SOURCE_OR_ALIAS_QC")
        self.assertEqual(report["scientific_status"],"HOLD")
        self.assertFalse(report["ready_for_ldsc"])
        self.assertEqual(report["error_type"],"ConnectionError")


if __name__=="__main__":
    unittest.main()
