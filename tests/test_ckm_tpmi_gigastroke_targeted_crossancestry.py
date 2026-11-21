"""Original GIGASTROKE rs187759247 exact GRCh37 effect direction gates."""
import math
import unittest
from scripts.ckm_stroke import probe_gigastroke_tpmi_43321_lead as p

def source_map():
    return {
       "schema":"ckm-TPMI43321-lead-variant-Ensembl-GRCh38-GRCh37-v1",
       "technical_status":"PASS_BIDIRECTIONAL_ENSEMBL_BUILD_MAP_AND_REF_ALLELE",
       "original_source_SHA256_verified":p.TPMI_SHA,
       "liftover_reverse_verified":True,
       "lead_GRCh38":{"reported_rsid":"rs187759247","chrom_GRCh38":"19",
                      "position_GRCh38":15462568},
       "mapped_GRCh37":{"chrom_GRCh37":"19","pos_GRCh37":15573379,
                        "REF_GRCh37":"T","A1_effect_allele_GRCh37":"C"},
       "scientific_status":"HOLD","ready_for_ldsc":False}
def original(ea="C",oa="T",beta="0.2",se="0.1",pval="0.05"):
    return ["19","15573379","0.03",beta,se,pval,"1.22","1.1","1.4",ea,oa]

class TestFourOutcomes(unittest.TestCase):
    def test_pinned_mapping_pass(self):
        self.assertTrue(p.gate_mapped(source_map()))
    def test_mutated_liftover_coordinate_is_hard_hold(self):
        x=source_map();x["mapped_GRCh37"]["pos_GRCh37"]+=1
        with self.assertRaisesRegex(ValueError,"drift"):
            p.gate_mapped(x)
    def test_identical_allele_uses_original_beta(self):
        r=p.effect_source(original())
        self.assertEqual(r["technical_harmonization_status"],
                         "PASS_ORIGINAL_ALLELE_DIRECTION_AND_NUMERIC")
        self.assertAlmostEqual(r["original_GIGASTROKE_beta_logOR_aligned_C"],0.2)
        self.assertTrue(r["effect_sign_concordant"])
        self.assertFalse(r["GIGASTROKE_original_samples_case_control_per_variant_known"])
    def test_swapped_effect_allele_negates_beta(self):
        r=p.effect_source(original(ea="T",oa="C",beta="0.2"))
        self.assertEqual(r["allele_orientation_relative_to_TPMI_A1_C"],"swapped")
        self.assertAlmostEqual(r["original_GIGASTROKE_beta_logOR_aligned_C"],-0.2)
        self.assertFalse(r["effect_sign_concordant"])
    def test_reverse_complement_same_orientation(self):
        r=p.effect_source(original(ea="G",oa="A"))
        self.assertEqual(r["allele_orientation_relative_to_TPMI_A1_C"],"complement_same")
        self.assertTrue(r["allele_aligned_effects_verified"])
    def test_incompatible_alleles_are_not_imputed(self):
        r=p.effect_source(original(ea="A",oa="C"))
        self.assertEqual(r["technical_harmonization_status"],"HOLD_ALLELE_ORIENTATION")
        self.assertFalse(r["allele_aligned_effects_verified"])
    def test_invalid_numeric_cannot_claim_direction(self):
        r=p.effect_source(original(beta="NA"))
        self.assertEqual(r["technical_harmonization_status"],"HOLD_EFFECT_NUMERIC_INVALID")
        self.assertFalse(r["allele_aligned_effects_verified"])
    def test_variant_position_wrong_fail_closed(self):
        x=original();x[1]="15573380"
        with self.assertRaisesRegex(ValueError,"wrong original"):
            p.effect_source(x)
    def test_no_effective_N_assumed_from_EAF(self):
        r=p.effect_source(original())
        self.assertEqual(r["original_effect_frequency"],0.03)
        self.assertFalse(r["GIGASTROKE_original_samples_case_control_per_variant_known"])
    def test_pinned_stroke_four_trait_accessions_exist(self):
        self.assertEqual(set(p.EXPECTED_TRAITS),set(p.PINNED))
if __name__=="__main__":
    unittest.main()
