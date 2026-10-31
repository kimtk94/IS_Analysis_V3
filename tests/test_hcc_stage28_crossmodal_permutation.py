"""Stage28 strict source-pinned permutation metrics and claim safeguards."""
import copy
import pathlib
import sys
import unittest

sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from verify_stage28_crossmodal_permutation import postflight,COMMIT,FULL,FULLZ

def fixture():
    shape={"M":3519,"CN":20339,"T":19768,"P":5585,"PP":7889}
    modal=[]
    for k,n in shape.items():
        x={"modality":k,"features_common_n":n,"train_feature_n":(n+1)//2,
           "heldout_feature_n":n//2,"sample_n":50,
           "train_map_is_one_to_one":True,
           "reverse_heldout_map_is_one_to_one":True,
           "train_vs_reverse_map_agreement_n":50,
           "unchanged_named_positions_n":2,
           "reordered_data_persisted":False,
           "sample_permutation_or_gene_identifiers_exported":False}
        for p in ("test","reverse_test"):
            x[p+"_finite_feature_n"]=n//2
            x[p+"_rho_before_median"]=0.02
            x[p+"_rho_after_median"]=0.97
            x[p+"_rho_after_q10"]=0.9
            x[p+"_rho_after_q90"]=0.99
            x[p+"_median_within_feature_rho_delta"]=0.95
            x[p+"_frac_rho_gt_095_before"]=0.0
            x[p+"_frac_rho_gt_095_after"]=0.88
            x[p+"_frac_rho_negative_after"]=0.01
        modal.append(x)
    return {"schema":"hcc-stage28-crossmodal-trainheldout-permutation-v1",
       "source_commit":COMMIT,"full_sha256":FULL,"fullz_sha256":FULLZ,
       "synthetic_permutation_and_heldout_recovery_passed":True,
       "each_modality_reported_independently":True,
       "technical_alignment_is_not_a_biological_replication":True,
       "publisher_original_biopsy_identity_reconciled":False,
       "original_suter_preprocessing_generation_reproduced":False,
       "sample_labels_persistently_corrected":False,
       "per_sample_mapping_and_per_gene_values_exported":False,
       "modalities":modal,
       "crossmodal_comparisons":[
           {"modality_a":a,"modality_b":b,
            "matching_sample_map_positions_n":50,
            "both_maps_bijective":True}
           for a,b in (("T","P"),("T","PP"),("P","PP"),("M","CN"))],
    }

class Stage28Tests(unittest.TestCase):
    def test_valid_review_is_still_partial(self):
        x=postflight(fixture())
        self.assertEqual(x["V03"],"PARTIAL")
        self.assertEqual(x["original_biopsy_identity_metadata_verification"],"HOLD")
        self.assertFalse(x["clinical_significance_claim_allowed"])
    def test_source_sha_drift_rejected(self):
        x=fixture();x["fullz_sha256"]="f"*64
        with self.assertRaises(ValueError):postflight(x)
    def test_missing_sample_rejected(self):
        x=fixture();x["modalities"][2]["sample_n"]=49
        with self.assertRaises(ValueError):postflight(x)
    def test_sample_mapping_relabel_forbidden(self):
        x=fixture();x["sample_labels_persistently_corrected"]=True
        with self.assertRaises(ValueError):postflight(x)
    def test_permutations_not_required_to_be_identical_for_valid_audit(self):
        x=fixture()
        x["crossmodal_comparisons"][0]["matching_sample_map_positions_n"]=1
        x["modalities"][2]["train_map_is_one_to_one"]=False
        self.assertEqual(postflight(x)["V03"],"PARTIAL")
    def test_impossible_crossmodal_count_rejected(self):
        x=fixture();x["crossmodal_comparisons"][2]["matching_sample_map_positions_n"]=51
        with self.assertRaises(ValueError):postflight(x)
    def test_no_export_of_raw_ids_or_mappings(self):
        x=fixture();x["per_sample_mapping_and_per_gene_values_exported"]=True
        with self.assertRaises(ValueError):postflight(x)

if __name__=="__main__":unittest.main()
