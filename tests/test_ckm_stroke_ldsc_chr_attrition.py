"""Fail-closed aggregate chromosome accounting tests; no network or LDSC."""
import unittest
from scripts.ckm_stroke.review_gigastroke_full_evidence import EXPECTED_BMI_CANDIDATE
from scripts.ckm_stroke.audit_ldsc_chr_attrition import audit_attrition
from tests.test_ckm_stroke_review_gigastroke_full_evidence import fixture

def with_chromosomes():
    raw=fixture()
    # Synthetic test only: distribute pinned candidate sites and outcome missing
    # across all 22 chromosomes, preserving independent full-source total counts.
    chrom_candidates=[46000]*21+[EXPECTED_BMI_CANDIDATE-46000*21]
    for item in raw["outcomes"]:
        missing=item["counts"]["candidate_coordinates_not_in_outcome"]
        rows=[]
        for index,total in enumerate(chrom_candidates,1):
            excluded=min(total,missing)
            missing-=excluded
            passed=total-excluded
            rows.append({
                "chromosome":index,
                "candidate_reference_bmi_snps":total,
                "unique_outcome_coordinate_snps":passed,
                "outcome_duplicate_coordinate_exclusions":0,
                "candidate_not_in_outcome":excluded,
                "allele_aligned_numeric_valid_pre_munge":passed,
                "nonaligned_or_numeric_invalid":0,
            })
        assert missing==0
        item["chr_attrition_pre_munge"]=rows
    return raw

class ChrAttritionTests(unittest.TestCase):
    def test_valid_22_chromosome_counts_are_not_ldsc_go(self):
        result=audit_attrition(with_chromosomes())
        self.assertEqual(len(result["outcomes"]),4)
        self.assertTrue(all(len(x["chromosomes"])==22 for x in result["outcomes"]))
        self.assertFalse(result["ready_for_ldsc"])
        self.assertFalse(result["approved_case_control_N_eff"])
        self.assertTrue(result["not_actual_ldsc_regression_counts"])

    def test_missing_chromosome_rejected(self):
        raw=with_chromosomes()
        raw["outcomes"][0]["chr_attrition_pre_munge"].pop()
        with self.assertRaisesRegex(ValueError,"22 per-chromosome"):
            audit_attrition(raw)

    def test_wrong_chromosome_or_duplicate_rejected(self):
        raw=with_chromosomes()
        raw["outcomes"][0]["chr_attrition_pre_munge"][1]["chromosome"]=1
        with self.assertRaisesRegex(ValueError,"order"):
            audit_attrition(raw)

    def test_per_chromosome_partition_mismatch(self):
        raw=with_chromosomes()
        raw["outcomes"][0]["chr_attrition_pre_munge"][0]["allele_aligned_numeric_valid_pre_munge"]-=1
        with self.assertRaisesRegex(ValueError,"partition"):
            audit_attrition(raw)

    def test_cross_subtype_reference_drift(self):
        raw=with_chromosomes()
        a=raw["outcomes"][1]["chr_attrition_pre_munge"]
        a[0]["candidate_reference_bmi_snps"]-=1
        a[0]["candidate_not_in_outcome"]-=1
        a[1]["candidate_reference_bmi_snps"]+=1
        a[1]["candidate_not_in_outcome"]+=1
        with self.assertRaisesRegex(ValueError,"inconsistent reference chromosome"):
            audit_attrition(raw)

    def test_boolean_or_negative_counts_rejected(self):
        raw=with_chromosomes()
        raw["outcomes"][0]["chr_attrition_pre_munge"][0]["candidate_reference_bmi_snps"]=True
        with self.assertRaisesRegex(ValueError,"invalid per-chromosome"):
            audit_attrition(raw)
        raw=with_chromosomes()
        raw["outcomes"][0]["chr_attrition_pre_munge"][0]["candidate_reference_bmi_snps"]=-1
        with self.assertRaisesRegex(ValueError,"invalid per-chromosome"):
            audit_attrition(raw)

    def test_false_md5_rejected(self):
        raw=with_chromosomes()
        raw["outcomes"][0]["raw_md5_verified"]=False
        with self.assertRaisesRegex(ValueError,"MD5"):
            audit_attrition(raw)

if __name__=="__main__":
    unittest.main()
