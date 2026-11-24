"""Stage25 patient-bootstrap & leave-one-out quantitative fail-closed tests."""
import copy
import pathlib
import sys
import unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from validate_ng_stage25_patient_influence import validate,SOURCE_SHA

def demo():
    a={
        "schema":"hcc-ng-stage24-descriptive-spearman-v1",
        "source_RDS_sha256":SOURCE_SHA,
        "mapping_scope":"source_RNA_x_current_UniProt_reviewed_strict_one_to_one",
        "crosswalk_protein_gene_pair_n":5343,
        "patient_mean_48":{"median":0.284086},
    }
    b={
        "schema":"hcc-ng-stage25-patient-influence-v1",
        "source_RDS_sha256":SOURCE_SHA,
        "mapping_scope":a["mapping_scope"],
        "strict_gene_protein_pair_n":5343,
        "patient_n":48,"leave_one_patient_out_iterations":48,
        "primary_median_spearman_rho":0.284086,
        "leave_one_out_median_min":0.280,"leave_one_out_median_max":0.289,
        "leave_one_out_median_q05":0.281,"leave_one_out_median_q95":0.288,
        "leave_one_out_max_absolute_shift":0.004914,
        "leave_one_out_max_shift_gt_0_05_review_flag":False,
        "patient_bootstrap_resamples":80,"patient_bootstrap_seed":20261031,
        "patient_bootstrap_median_rho_q025":0.20,
        "patient_bootstrap_median_rho_q975":0.35,
        "technical_preprocessing_robustness_proven":False,
        "original_ng_correlation_reproduced":False,
        "gene_level_FDR_or_per_gene_CI_calculated":False,
        "clinical_or_causal_validity_proven":False,
        "patient_gene_rho_or_expression_values_exported":False,
    }
    return a,b

class Tests(unittest.TestCase):
    def test_valid_remains_partial(self):
        a,b=demo()
        self.assertEqual(validate(a,b)["V07"],"PARTIAL")
        self.assertEqual(validate(a,b)["V03"],"PARTIAL")
    def test_unpinned_input_rejected(self):
        a,b=demo();b["source_RDS_sha256"]="wrong"
        with self.assertRaises(ValueError):validate(a,b)
    def test_missing_patient_repetition_rejected(self):
        a,b=demo();b["patient_n"]=50
        with self.assertRaises(ValueError):validate(a,b)
    def test_overstated_low_influence_rejected(self):
        a,b=demo();b["leave_one_out_max_absolute_shift"]=0.001
        with self.assertRaises(ValueError):validate(a,b)
    def test_bootstrap_outside_correlation_bounds_rejected(self):
        a,b=demo();b["patient_bootstrap_median_rho_q975"]=1.1
        with self.assertRaises(ValueError):validate(a,b)
    def test_individual_gene_export_or_clinical_pass_rejected(self):
        a,b=demo();b["patient_gene_rho_or_expression_values_exported"]=True
        with self.assertRaises(ValueError):validate(a,b)
        a,b=demo();b["clinical_or_causal_validity_proven"]=True
        with self.assertRaises(ValueError):validate(a,b)
    def test_inconsistent_review_flag_rejected(self):
        a,b=demo();b["leave_one_out_max_shift_gt_0_05_review_flag"]=True
        with self.assertRaises(ValueError):validate(a,b)

if __name__=="__main__":unittest.main()
