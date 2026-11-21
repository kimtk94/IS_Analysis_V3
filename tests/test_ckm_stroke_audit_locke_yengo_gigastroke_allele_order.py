"""Offline study-independent allele order fixtures, no original data downloaded."""
import math
import unittest
from scripts.ckm_stroke.audit_locke_yengo_gigastroke_allele_order import (
    compare_source_orders,compare_rows
)

def effect(ch,pos,ea="A",oa="C"):
    return ["1" if ch==1 else str(ch),str(pos),"0.2","0.1","0.05","0.04",
            str(math.exp(.1)),"1.01","1.3",ea,oa]

class IndependentOrderTests(unittest.TestCase):
    def setUp(self):
        self.yengo={(1,100):("rs1","A","C"),
                    (1,200):("rs2","A","C"),
                    (1,300):("rs3","A","C"),
                    (1,400):("rs4","A","C")}
        self.locke={"rs1":("A","C"),"rs2":("C","A"),
                    "rs3":("T","G"),"rs4":("A","C")}

    def test_original_bmi_allele_orders_are_not_assumed_identical(self):
        r=compare_source_orders(self.locke,self.yengo,strict=False)
        self.assertEqual(r["both_original_BMI_order_comparable"],4)
        self.assertEqual(r["categories"]["same"],2)
        self.assertEqual(r["categories"]["swapped"],1)
        self.assertEqual(r["categories"]["complement_same"],1)
        self.assertEqual(r["2015_order_reversed_vs_2018"],1)
        self.assertFalse(r["individual_participant_independence_confirmed"])

    def test_2015_absence_never_imputed_as_same(self):
        self.locke.pop("rs1")
        r=compare_source_orders(self.locke,self.yengo,strict=False)
        self.assertEqual(r["categories"]["not_in_independent_2015"],1)

    def test_outcome_same_yengo_but_2015_sign_reversal_detected(self):
        records=[effect(1,100),effect(1,200),effect(1,300),effect(1,400)]
        r=compare_rows(records,self.locke,self.yengo)
        self.assertEqual(r["counts"]["against_2018_same"],4)
        self.assertEqual(r["counts"]["against_2015_swapped"],1)
        self.assertEqual(r["GIGASTROKE_vs_Locke_sign_reversals"],1)
        self.assertEqual(r["counts"]["against_2015_complement_same"],1)

    def test_true_outcome_reversal_not_silently_counted_same(self):
        r=compare_rows([effect(1,100,"C","A")],self.locke,self.yengo)
        self.assertEqual(r["counts"]["against_2018_swapped"],1)
        self.assertEqual(r["GIGASTROKE_vs_Yengo_sign_reversals"],1)

    def test_duplicate_sampled_original_coordinate_rejected(self):
        with self.assertRaisesRegex(ValueError,"duplicate GWAS"):
            compare_rows([effect(1,100),effect(1,100)],self.locke,self.yengo)

    def test_nonmatching_origins_are_excluded_not_reoriented(self):
        r=compare_rows([effect(1,999),effect(1,100),effect(1,200)],self.locke,self.yengo)
        self.assertEqual(r["counts"]["no_2018_candidate"],1)
        self.assertEqual(r["counts"]["triple_intersection"],2)

    def test_bad_row_width_is_blocked(self):
        with self.assertRaisesRegex(ValueError,"column count"):
            compare_rows([["1","100"]],self.locke,self.yengo)

    def test_missing_independent_2015_does_not_validate_direction(self):
        self.locke.pop("rs1")
        r=compare_rows([effect(1,100)],self.locke,self.yengo)
        self.assertEqual(r["counts"]["no_2015_independent_candidate"],1)
        self.assertEqual(r["counts"]["triple_intersection"],0)

if __name__=="__main__":
    unittest.main()
