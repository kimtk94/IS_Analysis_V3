"""Offline negative contracts for GIGASTROKE meta-analysis sample-size safety."""
import unittest
from scripts.ckm_stroke.audit_gigastroke_ldsc_neff_method_gate import evaluate
from scripts.ckm_stroke.audit_gigastroke_supplement_cohort_neff import (
    EXPECTED_SHA256, EXPECTED_SUBTYPES, EXPECTED_CONTROLS
)
from scripts.ckm_stroke.audit_gigastroke_bmi_allele_harmonization import COLUMNS

MISS={"AIS":103,"LAA":0,"CE":28,"SVD":39}
COHORT_SUM={"AIS":195570.72,"LAA":23149.57,"CE":39605.59,"SVD":23758.72}

def fixtures():
    outcomes=[]; probes=[]
    for ph,ref in EXPECTED_SUBTYPES.items():
        n=ref["cases"]; m=MISS[ph]
        outcomes.append({
            "phenotype":ph,"accession":ref["accession"],
            "EUR_case_count_from_original_table1":n,
            "EUR_control_count_from_original_table1":EXPECTED_CONTROLS,
            "EUR_total_N_from_original_table1":ref["total"],
            "illustrative_neff_not_per_snp_or_approved":True,
            "this_pipeline_ldsc_executed":False,
            "illustrative_harmonic_neff":round(4*n*EXPECTED_CONTROLS/(n+EXPECTED_CONTROLS),3),
            "illustrative_cohortwise_neff":{
                "reported_case_total_from_rows":n,
                "cases_missing_cohort_control_count":m,
                "cases_with_reported_cohort_controls":n-m,
                "missing_cohort_control_rows":[{"cases":m,"cohort":"HVH2"}] if m else [],
                "not_per_variant_N_and_not_approved_for_LDSC":True,
                "cohortwise_sum_harmonic_neff_illustrative":COHORT_SUM[ph],
            }})
        probes.append({"phenotype":ph,"accession":ref["accession"],
                       "source_per_variant_N_available":False,
                       "source_rsid_available":False,
                       "source_header_columns":list(COLUMNS),
                       "outcome_beta_scale":"LOG_ODDS_ratio_approx_beta",
                       "sample_prefix_compressed_bytes":131072,
                       "bounded_sample_not_genomewide_proof":True,
                       "observed":{"parsed_effect_rows":600},
                       "maximum_abs_beta_minus_lnOR":1e-15})
    c={"schema":"ckm-gigastroke-corrected-primary-cohort-ldsc-gate-v1",
       "publisher_source_sha256":EXPECTED_SHA256,
       "technical_status":"PASS_PRIMARY_CORRECTED_COHORT_AND_TABLE3_ONLY",
       "publisher_original_table1_eur_cohort_count":44,
       "original_full_sumstats_downloaded":False,
       "ready_for_ldsc":False,"LDSC_h2_rg_executed":False,"scientific_status":"HOLD",
       "uk_biobank":{"AIS_case_count":2211,"UKB_reported_control_count":369419,
                     "overlap_in_BMI_actual_count_known":False},
       "outcomes":outcomes}
    e={"schema":"ckm-gigastroke-eur-bounded-original-logor-scale-v1",
       "technical_status":"PASS_FOUR_BOUNDED_LOGOR_EFFECT_SCALE_ONLY",
       "full_source_verified_in_this_run":False,
       "LDSC_h2_intercept_rg_executed":False,
       "ready_for_ldsc":False,"scientific_status":"HOLD","prefixes":probes}
    return c,e

class TestMetaNeffGate(unittest.TestCase):
    def test_original_counts_are_not_Neff(self):
        c,e=fixtures(); r=evaluate(c,e)
        self.assertEqual(len(r["outcomes"]),4)
        self.assertFalse(r["approved_N_eff_available"])
        self.assertFalse(r["ready_for_ldsc"])
        self.assertFalse(r["LDSC_h2_intercept_rg_run"])
        self.assertGreater(r["outcomes"][0]["illustrative_not_approved_for_LDSC"]["relative_difference_pct"],10)
        self.assertTrue(all(x["effective_N_model"]=="UNIDENTIFIED" for x in r["outcomes"]))

    def test_wrong_publisher_digest(self):
        c,e=fixtures(); c["publisher_source_sha256"]="0"*64
        with self.assertRaisesRegex(ValueError,"first-party"):evaluate(c,e)

    def test_wrong_subtype_accession(self):
        c,e=fixtures(); e["prefixes"][1]["accession"]="GCST90104541"
        with self.assertRaisesRegex(ValueError,"phenotype"):evaluate(c,e)

    def test_invented_variant_N(self):
        c,e=fixtures(); e["prefixes"][0]["source_per_variant_N_available"]=True
        with self.assertRaisesRegex(ValueError,"per-variant"):evaluate(c,e)

    def test_unknown_overlap_cannot_be_marked_known(self):
        c,e=fixtures(); c["uk_biobank"]["overlap_in_BMI_actual_count_known"]=True
        with self.assertRaisesRegex(ValueError,"UKB"):evaluate(c,e)

    def test_missing_controls_not_zero(self):
        c,e=fixtures()
        c["outcomes"][0]["illustrative_cohortwise_neff"]["cases_missing_cohort_control_count"]=0
        with self.assertRaisesRegex(ValueError,"missing cohort controls"):evaluate(c,e)

    def test_cohortwise_scenario_not_silently_pooled(self):
        c,e=fixtures()
        c["outcomes"][0]["illustrative_cohortwise_neff"]["cohortwise_sum_harmonic_neff_illustrative"]=999999
        with self.assertRaisesRegex(ValueError,"scenario"):evaluate(c,e)

    def test_pooled_arithmetic_not_approved(self):
        c,e=fixtures()
        c["outcomes"][0]["illustrative_neff_not_per_snp_or_approved"]=False
        with self.assertRaisesRegex(ValueError,"sample metadata"):evaluate(c,e)

    def test_unverified_effect_scale_rejected(self):
        c,e=fixtures(); e["prefixes"][0]["outcome_beta_scale"]="OR"
        with self.assertRaisesRegex(ValueError,"effect-source"):evaluate(c,e)

    def test_nonfinite_scenario_rejected(self):
        c,e=fixtures()
        c["outcomes"][0]["illustrative_cohortwise_neff"]["cohortwise_sum_harmonic_neff_illustrative"]=float("nan")
        with self.assertRaisesRegex(ValueError,"nonfinite"):evaluate(c,e)

if __name__=="__main__":
    unittest.main()
