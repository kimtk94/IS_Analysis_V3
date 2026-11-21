"""No-network source-backed LDSC reference-weight and cohortwise Neff gates."""
from __future__ import annotations
from unittest.mock import patch
import unittest
from scripts.ckm_stroke.audit_ldsc_weight_reference_merge import (
    reconcile_reference,verify_neff_source,EXPECTED_SUPPLEMENT_SHA256,
)

def fake_cohort():
    case={"AIS":62100,"LAA":6399,"CE":10804,"SVD":6811}
    names={"AIS":"GCST90104540","LAA":"GCST90104542","CE":"GCST90104541","SVD":"GCST90104543"}
    d={"technical_status":"PASS_PRIMARY_CORRECTED_COHORT_AND_TABLE3_ONLY",
       "publisher_source_sha256":EXPECTED_SUPPLEMENT_SHA256,
       "publisher_original_table1_eur_cohort_count":44,
       "ready_for_ldsc":False,"LDSC_h2_rg_executed":False,"outcomes":[]}
    for name,n in case.items():
        controls=1234808
        theoretical=round(4*n*controls/(n+controls),3)
        c={"cohortwise_sum_harmonic_neff_illustrative":theoretical*.8,
           "not_per_variant_N_and_not_approved_for_LDSC":True,
           "cases_missing_cohort_control_count":0,
           "cohort_count_with_cases_and_reported_controls":10}
        d["outcomes"].append({
          "phenotype":name,"accession":names[name],
          "EUR_case_count_from_original_table1":n,
          "EUR_control_count_from_original_table1":controls,
          "EUR_total_N_from_original_table1":n+controls,
          "illustrative_harmonic_neff":theoretical,
          "illustrative_cohortwise_neff":c,
          "illustrative_neff_not_per_snp_or_approved":True,
        })
    return d

def fake_original():
    base=[{"candidate_reference_bmi_snps":1 if i<3 else 0} for i in range(22)]
    return {"outcomes":[{"chr_attrition_pre_munge":base}]}

def fake_review():
    return {"source_BMI_candidate_count":3,
            "outcomes":[{"phenotype":ph,"totals":{"allele_aligned_numeric_valid_pre_munge":2}}
                        for ph in ("AIS","LAA","CE","SVD")]}

class ReferenceMergeTests(unittest.TestCase):
    def test_cohortwise_and_pooled_not_approved(self):
        result=verify_neff_source(fake_cohort())
        self.assertEqual(len(result),4)
        self.assertTrue(all(x["difference_percent"]>0 for x in result))
        self.assertTrue(all(x["N_SNP_specific_observed"] is False for x in result))
        self.assertTrue(all(x["scientific_status"]=="HOLD" for x in result))

    def test_original_cohort_sha_tamper_fails_closed(self):
        data=fake_cohort()
        data["publisher_source_sha256"]="0"*64
        with self.assertRaisesRegex(ValueError,"44-cohort"):
            verify_neff_source(data)

    def test_missing_control_metadata_fails_closed(self):
        data=fake_cohort()
        data["outcomes"][0]["illustrative_cohortwise_neff"]["cases_missing_cohort_control_count"]=-3
        with self.assertRaisesRegex(ValueError,"evidence mismatch"):
            verify_neff_source(data)

    @patch("scripts.ckm_stroke.audit_ldsc_weight_reference_merge.audit_attrition")
    @patch("scripts.ckm_stroke.audit_ldsc_weight_reference_merge.verify_neff_source")
    @patch("scripts.ckm_stroke.audit_ldsc_weight_reference_merge.EXPECTED_SNP_COUNT",3)
    def test_real_22chr_ids_intersection_still_scientific_hold(self,_expected,neff,source_review):
        source_review.return_value=fake_review()
        neff.return_value=[{"phenotype":"AIS","difference_percent":17.3}]
        ids={"rs1","rs2","rs3"}
        ld={i:set() for i in range(1,23)}
        ld[1]={"rs1"};ld[2]={"rs2"};ld[3]={"rs3"}
        report=reconcile_reference(ids,ld,ids,fake_original(),fake_cohort())
        self.assertEqual(report["reference_BMI_candidate_22chr_merge_count"],3)
        self.assertEqual(len(report["chr_reference"]),22)
        self.assertFalse(report["case_control_N_eff_per_variant_known"])
        self.assertTrue(report["reference_must_NOT_be_called_final_regression_SNP_count"])
        self.assertFalse(report["ready_for_ldsc"])

    @patch("scripts.ckm_stroke.audit_ldsc_weight_reference_merge.audit_attrition")
    @patch("scripts.ckm_stroke.audit_ldsc_weight_reference_merge.verify_neff_source")
    @patch("scripts.ckm_stroke.audit_ldsc_weight_reference_merge.EXPECTED_SNP_COUNT",3)
    def test_missing_reference_chr_snp_fails(self,_expected,neff,source_review):
        source_review.return_value=fake_review()
        ids={"rs1","rs2","rs3"}
        ref={i:set() for i in range(1,23)}
        ref[1]={"rs1"};ref[2]={"rs2"}
        with self.assertRaisesRegex(ValueError,"absent"):
            reconcile_reference(ids,ref,ids,fake_original(),fake_cohort())

    @patch("scripts.ckm_stroke.audit_ldsc_weight_reference_merge.audit_attrition")
    @patch("scripts.ckm_stroke.audit_ldsc_weight_reference_merge.verify_neff_source")
    @patch("scripts.ckm_stroke.audit_ldsc_weight_reference_merge.EXPECTED_SNP_COUNT",3)
    def test_different_reference_chromosome_bins_fail(self,_expected,neff,source_review):
        source_review.return_value=fake_review()
        ids={"rs1","rs2","rs3"}
        ref={i:set() for i in range(1,23)}
        ref[1]={"rs1","rs2"};ref[3]={"rs3"}
        with self.assertRaisesRegex(ValueError,"chr attrition"):
            reconcile_reference(ids,ref,ids,fake_original(),fake_cohort())

    @patch("scripts.ckm_stroke.audit_ldsc_weight_reference_merge.audit_attrition")
    @patch("scripts.ckm_stroke.audit_ldsc_weight_reference_merge.verify_neff_source")
    @patch("scripts.ckm_stroke.audit_ldsc_weight_reference_merge.EXPECTED_SNP_COUNT",3)
    def test_reference_SNP_in_multiple_chr_is_not_same_as_unique_merge(self,_expected,neff,source_review):
        source_review.return_value=fake_review()
        ids={"rs1","rs2","rs3"}
        ref={i:set() for i in range(1,23)}
        ref[1]={"rs1"};ref[2]={"rs1","rs2"};ref[3]={"rs3"}
        with self.assertRaisesRegex(ValueError,"does not reconcile"):
            reconcile_reference(ids,ref,ids,fake_original(),fake_cohort())

if __name__=="__main__":
    unittest.main()
