"""No-network independent original allele order reviewer and matched-coordinate tests."""
import gzip
import io
import unittest
from scripts.ckm_stroke.audit_gigastroke_bmi_allele_harmonization import examine_outcome_rows
from scripts.ckm_stroke.audit_gigastroke_independent_full_direction import direction_counts
from tests.test_ckm_stroke_gigastroke_bmi_allele_harmonization import fake_data,panel
from tests.test_ckm_stroke_review_gigastroke_full_evidence import fixture

def synthetic_full():
    data=fixture()
    data["independent_locke_panel_evidence"]={
      "candidate_rows":1070753,
      "plain_sha256":"a1cc5ec06cc81d2bb0a583eb76c0f7e81784cd5dc3f622270de61023f95de4bf",
      "fractional_N_rows":11506,"selected_without_Yengo_2018":True}
    data["independent_2015_vs_2018_order_evidence"]={
      "total_Yengo_candidate":1007899,
      "both_original_BMI_order_comparable":1007858,
      "2015_order_reversed_vs_2018":491428,
      "individual_participant_independence_confirmed":False,
      "categories":{"same":516430,"swapped":491428,"not_in_independent_2015":41}}
    for source in data["outcomes"]:
        aligned=source["counts"]["unique_outcome_coordinate_matches"]
        values={"2015_same":aligned//2-20,
                "2015_swapped":aligned-aligned//2,
                "not_in_locke_2015":20}
        by_chr=[{"chr":i+1,"candidate_positions":0,"categories":{}} for i in range(22)]
        for key,value in values.items():
            for i in range(22):
                n=value//22+(i<value%22)
                if n:by_chr[i]["categories"][key]=n
                by_chr[i]["candidate_positions"]+=n
        source["locke_2015_independent_full_orientation"]={
          "unique_original_outcome_candidate_positions":aligned,
          "categories":values,
          "aligned_to_both_2015_and_2018":aligned-20,
          "2015_effect_sign_reversals":values["2015_swapped"],
          "by_chr":by_chr,"complete_original_source_scanned":True,
          "source_MD5_validated_by_caller":True,
          "approved_LDSC_analysis":False}
        remaining=source["counts"]["candidate_coordinates_not_in_outcome"]
        chrlist=[]
        for i in range(22):
            capacity=1007899//22+(i<1007899%22)
            excluded=min(remaining,capacity)
            remaining-=excluded
            kept=capacity-excluded
            chrlist.append({"chromosome":i+1,"candidate_reference_bmi_snps":capacity,
                            "unique_outcome_coordinate_snps":kept,
                            "outcome_duplicate_coordinate_exclusions":0,
                            "candidate_not_in_outcome":excluded,
                            "allele_aligned_numeric_valid_pre_munge":kept,
                            "nonaligned_or_numeric_invalid":0})
        source["chr_attrition_pre_munge"]=chrlist
    return data

class TestFullIndependentDirection(unittest.TestCase):
    def test_variant_flip_and_duplicate_quarantine(self):
        source={"rs1":("C","A"),"rs2":("A","C"),"rs3":("A","C"),
                "rs4":("A","C"),"rs5":("A","C"),"rs6":("A","C")}
        raw,_=fake_data()
        with gzip.GzipFile(fileobj=io.BytesIO(raw),mode="rb") as file:
            r=examine_outcome_rows((line.decode() for line in file),panel(),locke_alleles=source)
        a=r["locke_2015_independent_full_orientation"]
        self.assertEqual(a["2015_effect_sign_reversals"],3)
        self.assertEqual(a["categories"]["2015_swapped"],2)
        self.assertEqual(a["categories"]["2015_complement_same"],1)
        self.assertEqual(a["categories"]["2015_complement_swapped"],1)
        self.assertEqual(a["categories"]["yengo_numeric_or_allele_invalid"],1)
        self.assertEqual(sum(x["candidate_positions"] for x in a["by_chr"]),5)

    def test_missing_locke_is_not_assumed_same(self):
        raw,_=fake_data()
        with gzip.GzipFile(fileobj=io.BytesIO(raw),mode="rb") as file:
            r=examine_outcome_rows((line.decode() for line in file),panel(),
                                   locke_alleles={"rs1":("A","C")})
        self.assertEqual(r["locke_2015_independent_full_orientation"]["categories"]["not_in_locke_2015"],3)

    def test_full_source_qc_pass_but_science_hold(self):
        r=direction_counts(synthetic_full())
        self.assertEqual(len(r["outcomes"]),4)
        self.assertFalse(r["approved_outcome_N_eff"])
        self.assertFalse(r["ready_for_ldsc"])
        self.assertFalse(r["LDSC_h2_rg_MR_coloc_executed"])

    def test_md5_failure_cannot_pass(self):
        r=synthetic_full();r["outcomes"][0]["raw_md5_verified"]=False
        with self.assertRaisesRegex(ValueError,"MD5"):
            direction_counts(r)

    def test_bad_locke_sha_cannot_pass(self):
        r=synthetic_full();r["independent_locke_panel_evidence"]["plain_sha256"]="0"*64
        with self.assertRaisesRegex(ValueError,"independent Locke"):
            direction_counts(r)

    def test_modified_source_study_direction_rejected(self):
        r=synthetic_full();r["independent_2015_vs_2018_order_evidence"]["2015_order_reversed_vs_2018"]-=1
        with self.assertRaisesRegex(ValueError,"independent 2015/2018"):
            direction_counts(r)

    def test_missing_chr_rejected(self):
        r=synthetic_full();r["outcomes"][0]["locke_2015_independent_full_orientation"]["by_chr"].pop()
        with self.assertRaisesRegex(ValueError,"22-chromosome"):
            direction_counts(r)

    def test_global_direction_count_drift_rejected(self):
        r=synthetic_full();r["outcomes"][0]["locke_2015_independent_full_orientation"]["categories"]["2015_swapped"]-=1
        with self.assertRaisesRegex(ValueError,"not all"):
            direction_counts(r)

    def test_boolean_counts_rejected(self):
        r=synthetic_full();r["outcomes"][0]["locke_2015_independent_full_orientation"]["categories"]["2015_swapped"]=True
        with self.assertRaisesRegex(ValueError,"invalid full-source"):
            direction_counts(r)

if __name__=="__main__":
    unittest.main()
