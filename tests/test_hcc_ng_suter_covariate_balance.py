"""Source-specific synthetic patient-level confounder audit tests."""
import pathlib,sys,unittest
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"/"hcc"))
from probe_ng_suter_covariate_balance import audit,canonical_columns
from test_hcc_ng_suter_clinical_sample_join import fixture,COLS

def synthetic():
    labels,rows=fixture()
    cols={**COLS,"underlyingliverdisease":4,"age":5,
          "sex":6,"bclc":7,"treatment":8}
    for k,row in enumerate(rows):
        row[4]="ALD" if k<12 else "HCV" if k<23 else "HBV"
        # Both duplicated patient biopsies at k>=48 need stable source values:
        if k>=48:row[4]="ALD"
        row[5]="60" if k<12 or k>=48 else "70"
        row[6]="M" if k<12 or k>=48 else "F"
        row[7]="A" if k<12 or k>=48 else "C"
        row[8]="Unknown" # presence not proof of treatment control
    return labels,rows,cols

class TestCovariate(unittest.TestCase):
    def test_patient_based_n_and_covariates(self):
        a,b,c=synthetic()
        report=audit(a,b,c)
        self.assertEqual(report["source_unique_patients_n"],48)
        self.assertEqual(report["exclusive_group_counts"],{"ALD":12,"HCV":11})
        self.assertEqual(report["groups"]["ALD"]["age"]["median"],60.0)
        self.assertEqual(report["groups"]["HCV"]["age"]["median"],70.0)
        self.assertEqual(report["ALD_minus_HCV_unadjusted_age_median_years"],-10.0)
        self.assertFalse(report["full_differential_omics_GO"])
        self.assertIsNone(report["absolute_standardized_mean_difference"]["male_sex"])
        self.assertTrue(report["complete_group_separation_no_finite_SMD"]["male_sex"])
        self.assertTrue(report["complete_group_separation_no_finite_SMD"]["BCLC_early_0_A"])
        self.assertTrue(report["imbalance_attention_threshold_0_2"]["BCLC_early_0_A"])
        self.assertNotIn("patient0",str(report).lower())

    def test_real_child_column_is_recognized_but_no_treatment_assumption(self):
        a,b,c=synthetic()
        c["child"]=9
        for row in b:row[9]="A"
        out=audit(a,b,c)
        self.assertTrue(out["source_child_pugh_header_detected"])
        self.assertEqual(out["groups"]["ALD"]["child_pugh"]["present_n"],12)
        self.assertFalse(out["candidate_baseline_adjustment_ready"])

    def test_child_pugh_low_coverage_prevents_misleading_smd(self):
        a,b,c=synthetic()
        c["child"]=9
        for i,row in enumerate(b):
            row[9]="A" if i<8 or i>=48 else ("B" if 12<=i<15 else "")
        report=audit(a,b,c)
        self.assertTrue(report["ChildPugh_effect_smd_withheld_for_low_coverage"])
        self.assertIsNone(report["absolute_standardized_mean_difference"]["ChildPugh_A_vs_BC"])
        self.assertFalse(report["candidate_baseline_adjustment_ready"])

    def test_missing_age_not_fabricated(self):
        a,b,c=synthetic()
        c.pop("age")
        report=audit(a,b,c)
        self.assertIn("age",report["unavailable_clinical_fields"])
        self.assertNotIn("age",report["groups"]["ALD"])

    def test_source_columns_dont_use_unknown_names(self):
        allowed=canonical_columns({"ageatclassification":3,"participantname":8,
                                    "bclc":10,"rnaqc":11})
        self.assertEqual(allowed,{"age":3,"bclc":10})

    def test_malformed_50_sample_identity_is_rejected(self):
        a,b,c=synthetic()
        a[0]="unknown"
        with self.assertRaises(ValueError):audit(a,b,c)

    def test_multibiopsy_field_disagreement_not_assumed_balanced(self):
        a,b,c=synthetic()
        b[48][5]="82"
        report=audit(a,b,c)
        self.assertTrue(report["PILOT_ONLY"])
        self.assertEqual(report["groups"]["ALD"]["age"]["nonmissing_n"],11)

if __name__=="__main__":unittest.main()
