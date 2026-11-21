"""EAF and GC input diagnostics: synthetic only, no invented N_eff or GC coefficient."""
from __future__ import annotations
from collections import Counter
import gzip
import io
import unittest
from scripts.ckm_stroke.audit_gigastroke_eaf_gc import (
    diagnostic,summarize,review_diagnostics,DIAGNOSTIC_KEYS,MAF_KEYS,HRS_KEYS,
    WALD_KEYS,LOGOR_KEYS,EAF_KEYS
)
from scripts.ckm_stroke.audit_gigastroke_bmi_allele_harmonization import examine_outcome_rows
from tests.test_ckm_stroke_gigastroke_bmi_allele_harmonization import panel,fake_data
from tests.test_ckm_stroke_gigastroke_independent_full_direction import synthetic_full

def row(eaf="0.2",beta="0.1",se="0.05",p="0.0455003",odds="1.1051709180756477"):
    return ["1","100",eaf,beta,se,p,odds,"1.01","1.3","A","C"]

def full():
    r=synthetic_full()
    for outcome in r["outcomes"]:
        chrlist=outcome["chr_attrition_pre_munge"]
        by_chr=[]
        global_counts=Counter()
        for entry in chrlist:
            n=entry["allele_aligned_numeric_valid_pre_munge"]
            c={"maf_ge_005":n,"hrs_absdiff_lt_002":n,
               "abs_log10p_delta_lt_001":n,
               "beta_logOR_absdelta_lt_1e10":n,
               "raw_eaf_strict_between_zero_and_one":n}
            by_chr.append({"chr":entry["chromosome"],"matched_diagnostic_SNPs":n,
                           "categories":c})
            global_counts.update(c)
        outcome["eaf_gc_source_descriptive"]={
           "matched_original_source_SNPs_with_valid_EAF_and_Wald_inputs":
                outcome["counts"]["unique_harmonizable_BMI_positions"],
           "categories":dict(global_counts),"by_chr":by_chr,
           "HRS_frequency_is_reference_only":True,
           "HRS_ancestry_and_cohort_vs_GIGASTROKE_frequency_equivalence_not_assumed":True,
           "Wald_P_discrepancies_are_NOT_proof_of_genomic_control":True,
           "missing_per_variant_N_eff":True,
           "source_original_effect_sizes_uploaded":False}
    return r

class EAFGCInputTests(unittest.TestCase):
    def test_allele_flip_reorients_frequency_only_for_hrs_comparison(self):
        self.assertEqual(diagnostic(row(),0.2,"same")[1],"hrs_absdiff_lt_002")
        self.assertEqual(diagnostic(row(),0.8,"swapped")[1],"hrs_absdiff_lt_002")
        self.assertEqual(diagnostic(row(),0.8,"same")[1],"hrs_absdiff_ge_020")

    def test_reference_missing_is_not_filled_with_outcome_eaf(self):
        self.assertEqual(diagnostic(row(),None,"same")[1],"hrs_reference_not_available")

    def test_wald_and_source_p_coherent_when_from_same_z(self):
        r=diagnostic(row(),0.2,"same")
        self.assertEqual(r[2],"abs_log10p_delta_lt_001")
        self.assertEqual(r[3],"beta_logOR_absdelta_lt_1e10")

    def test_discordant_p_does_not_trigger_inferred_genomic_control(self):
        data=row(p="1e-100")
        r=diagnostic(data,0.2,"same")
        self.assertEqual(r[2],"abs_log10p_delta_ge_1")

    def test_underflowed_original_p_counted_separately(self):
        self.assertEqual(diagnostic(row(p="0"),0.2,"same")[2],"reported_p_is_zero")

    def test_invalid_reference_hrs_and_numeric_allele_rejected(self):
        with self.assertRaisesRegex(ValueError,"HRS"):
            diagnostic(row(),-0.1,"same")
        with self.assertRaisesRegex(ValueError,"frequency"):
            diagnostic(row(eaf="1.3"),0.5,"same")
        with self.assertRaisesRegex(ValueError,"orientation"):
            diagnostic(row(),0.5,"allele_incompatible")

    def test_source_eaf_boundaries_and_low_maf(self):
        r=diagnostic(row(eaf="0.0"),0.0,"same")
        self.assertEqual(r[0],"maf_lt_001")
        self.assertEqual(r[4],"raw_eaf_boundary_zero_or_one")
        self.assertEqual(diagnostic(row(eaf="0.04"),0.04,"same")[0],"maf_001_to_005")

    def test_quarantined_duplicate_not_double_counted_in_full_eaf(self):
        raw,_=fake_data()
        with gzip.GzipFile(fileobj=io.BytesIO(raw),mode="rb") as stream:
            r=examine_outcome_rows((x.decode() for x in stream),panel(),
                hrs_by_coord={(ch,pos):0.2 for ch,pos in panel()})
        e=r["eaf_gc_source_descriptive"]
        self.assertEqual(e["matched_original_source_SNPs_with_valid_EAF_and_Wald_inputs"],4)
        self.assertEqual(sum(x["matched_diagnostic_SNPs"] for x in e["by_chr"]),4)
        self.assertTrue(e["Wald_P_discrepancies_are_NOT_proof_of_genomic_control"])

    def test_independent_full_md5_eaf_qc_pass_but_science_hold(self):
        r=review_diagnostics(full())
        self.assertEqual(len(r["outcomes"]),4)
        self.assertEqual(r["technical_status"],"PASS_ORIGINAL_MD5_FOUR_EAF_AND_WALD_DESCRIPTIVE_ONLY")
        self.assertFalse(r["N_eff_approved"])
        self.assertFalse(r["GC_method_or_correction_factor_INFERRED"])
        self.assertFalse(r["LDSC_h2_intercept_rg_run"])
        self.assertFalse(r["ready_for_ldsc"])

    def test_md5_tamper_rejected(self):
        r=full();r["outcomes"][0]["raw_md5_verified"]=False
        with self.assertRaisesRegex(ValueError,"MD5"):
            review_diagnostics(r)

    def test_hrs_partition_and_per_chromosome_corruption_rejected(self):
        r=full()
        e=r["outcomes"][1]["eaf_gc_source_descriptive"]
        e["categories"]["hrs_absdiff_lt_002"]-=1
        with self.assertRaisesRegex(ValueError,"whole-genome"):
            review_diagnostics(r)
        r=full()
        e=r["outcomes"][1]["eaf_gc_source_descriptive"]
        e["by_chr"][21]["categories"]["maf_ge_005"]-=1
        with self.assertRaisesRegex(ValueError,"chromosome"):
            review_diagnostics(r)

    def test_outcome_snp_ledger_mismatch_rejected(self):
        r=full()
        e=r["outcomes"][0]["eaf_gc_source_descriptive"]
        e["matched_original_source_SNPs_with_valid_EAF_and_Wald_inputs"]-=1
        with self.assertRaisesRegex(ValueError,"valid frequency"):
            review_diagnostics(r)

    def test_gc_source_scope_flag_removed_is_held(self):
        r=full();r["outcomes"][0]["eaf_gc_source_descriptive"]["Wald_P_discrepancies_are_NOT_proof_of_genomic_control"]=False
        with self.assertRaisesRegex(ValueError,"assumptions"):
            review_diagnostics(r)

if __name__=="__main__":
    unittest.main()
