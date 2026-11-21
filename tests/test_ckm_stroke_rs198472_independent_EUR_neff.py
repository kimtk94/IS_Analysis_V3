"""Network-free population allele orientation and case-control Neff HOLD tests."""
from __future__ import annotations
import copy
import unittest
from scripts.ckm_stroke.verify_rs198472_independent_EUR_and_neff import (
    verify_giant_tested_allele,verified_target_artifact,parse_ensembl_variant,
    effective_n_scenario,audit,POPULATION,EXPECTED_EAF,VARIANT,BP,OUTCOMES
)
def giant(a1="A"):
    return {"rsid":VARIANT,"tested_allele":a1,
            "other_allele":"G" if a1=="A" else "A",
            "HRS_tested_allele_frequency":0.5341}
def ensembl():
    return {"name":"rs198472",
            "mappings":[{"assembly_name":"GRCh37","seq_region_name":"11",
                         "start":61527266,"end":61527266,"allele_string":"A/G"}],
            "populations":[{"population":POPULATION,"allele":"A","frequency":0.58},
                           {"population":POPULATION,"allele":"G","frequency":0.42},
                           {"population":"1000GENOMES:phase_3:EAS","allele":"A","frequency":0.71}]}
def source():
    recs=[]
    for ph in OUTCOMES:
        recs.append({"chromosome":11,"position_GRCh37":61527266,
                     "rsid_from_pinned_GIANT":"rs198472",
                     "source_reference_alleles_sorted":"AG",
                     "GIANT_HRS_tested_allele_frequency":0.5341,
                     "GIGASTROKE_frequency_oriented_to_GIANT":EXPECTED_EAF[ph]})
    return {
        "schema":"ckm-GIGASTROKE-EUR-chr11-variant-identity-chr4-P-Wald-v1",
        "technical_status":"PASS_PINNED_FULL_MD5_TARGET_CHR4_11_VARIANT_IDENTITY",
        "all_four_chr11_outlier_same_genomic_variant":True,
        "full_source_reverified_four_original_MD5_and_SHA256":True,
        "outcome_per_variant_N_eff_approved":False,
        "scientific_status":"HOLD",
        "chr11_source_outlier_frequency_cross_trait_details":recs,
        "chr11_outlier_variant_identity":recs[0],
    }
class TestIndependentEUR(unittest.TestCase):
    def test_orientation_allele_A_and_G_1000G_EUR(self):
        r=parse_ensembl_variant(ensembl(),giant("A"))
        self.assertEqual(r["eur_AF_of_GIANT_tested_allele"],0.58)
        r=parse_ensembl_variant(ensembl(),giant("G"))
        self.assertEqual(r["eur_AF_of_GIANT_tested_allele"],0.42)
        self.assertTrue(r["GRCh37_alleles_and_coordinate_verified"])
    def test_only_phase3_EUR_not_EAS_used(self):
        r=parse_ensembl_variant(ensembl(),giant())
        self.assertEqual(r["population"],POPULATION)
        self.assertNotEqual(r["eur_AF_of_GIANT_tested_allele"],0.71)
    def test_grch38_coordinate_or_allele_confusion_rejected(self):
        e=ensembl();e["mappings"][0]["assembly_name"]="GRCh38"
        with self.assertRaisesRegex(ValueError,"GRCh37"):
            parse_ensembl_variant(e,giant())
        e=ensembl();e["mappings"][0]["allele_string"]="A/T"
        with self.assertRaisesRegex(ValueError,"GRCh37"):
            parse_ensembl_variant(e,giant())
    def test_missing_european_or_frequency_sum_drift_rejected(self):
        e=ensembl();e["populations"]=[e["populations"][2]]
        with self.assertRaisesRegex(ValueError,"biallelic"):
            parse_ensembl_variant(e,giant())
        e=ensembl();e["populations"][1]["frequency"]=0.90
        with self.assertRaisesRegex(ValueError,"biallelic"):
            parse_ensembl_variant(e,giant())
    def test_publisher_2018_tested_allele_from_completed_source(self):
        def fake():
            yield 1,["11",str(BP),"rs198472","A","G","0.5341","0.2","0.1","0.1","200"]
            yield 2,["11","777","rs123","A","C","0.1","0.1","0.1","0.5","200"]
            yield 0,["__SOURCE_COMPLETE_SHA256_PASS__","46754937"]
        r=verify_giant_tested_allele(fake,expected_source_rows=2)
        self.assertEqual(r["tested_allele"],"A")
        self.assertEqual(r["HRS_tested_allele_frequency"],0.5341)
    def test_missing_publisher_sha_completion_rejected(self):
        def fake():
            yield 1,["11",str(BP),"rs198472","A","G","0.5341","0.2","0.1","0.1","200"]
        with self.assertRaisesRegex(ValueError,"sentinel"):
            verify_giant_tested_allele(fake,expected_source_rows=1)
    def test_bad_original_source_rs_or_frequency_fails_closed(self):
        s=source();s["chr11_source_outlier_frequency_cross_trait_details"][1]["rsid_from_pinned_GIANT"]="rs0"
        with self.assertRaisesRegex(ValueError,"source allele/frequency drift"):
            verified_target_artifact(s)
    def test_full_comparison_keeps_scientific_hold(self):
        r=audit(source(),ensembl(),giant())
        self.assertEqual(r["technical_status"],"PASS_GRCh37_1000G_EUR_SOURCE_BIALLELIC_VARIANT_MATCH")
        self.assertEqual(len(r["GIGASTROKE_allele_aligned_subtype_comparisons"]),4)
        self.assertFalse(r["LDSC_effective_sample_size_evidence"]["constant_illustrative_N_eff_approvable_for_LDSC"])
        self.assertFalse(r["ready_for_ldsc"])
        self.assertFalse(r["LDSC_h2_rg_MR_coloc_executed"])
    def test_sample_size_total_and_illustration_not_approved(self):
        n=effective_n_scenario()
        self.assertEqual([x["publisher_reported_total_N"] for x in n["illustrative_case_control_harmonic_N_eff_by_trait"]],
                         [1296908,1241207,1245612,1241619])
        self.assertFalse(n["per_variant_N_eff_available"])
        self.assertEqual(n["N_eff_scientific_status"],"HOLD")
if __name__=="__main__":
    unittest.main()
