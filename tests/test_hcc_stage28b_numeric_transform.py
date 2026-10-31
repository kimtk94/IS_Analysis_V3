"""Stage28b quantitative reconstruction source and claim guards."""
import pathlib,sys,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from verify_stage28b_numeric_transform import verify,FULL,FULLZ,COMMIT
def fixture():
    ns={"T":19768,"P":5585,"PP":7889}
    layers=[]
    for key,n in ns.items():
        layers.append({"modality":key,"common_features_n":n,
          "training_features_n":(n+1)//2,"test_features_n":n//2,
          "sample_n":50,"mapping_bijective":True,
          "unchanged_named_positions_n":2,
          "heldout_median_abs_z_residual_before":0.8,
          "heldout_median_abs_z_residual_after":0.0,
          "heldout_p95_abs_z_residual_after":0.0,
          "heldout_max_abs_z_residual_after":0.0000001,
          "heldout_cell_fraction_abs_lt_1e6_after":1.0,
          "heldout_cell_fraction_abs_lt_1e4_after":1.0,
          "heldout_feature_fraction_row_max_abs_lt_1e5_after":1.0,
          "sample_feature_identifiers_exported":False})
    return {
      "schema":"hcc-stage28b-fullz-numeric-reconstruction-v1",
      "source_commit":COMMIT,"full_sha256":FULL,"fullz_sha256":FULLZ,
      "synthetic_numeric_reconstruction_passed":True,"layers":layers,
      "source_original_Ng_tumor_normal_2022_values_reproduced":False,
      "metadata_verified_canonical_sample_ids":False,
      "source_files_overwritten":False,
      "public_artifact_contains_sample_gene_or_permutation_vector":False}
class TestStage28B(unittest.TestCase):
    def test_technical_reconstruction_not_biological_claim(self):
        x=verify(fixture());self.assertEqual(x["V03"],"PARTIAL")
        self.assertEqual(x["clinical_patient_biopsy_correct_label_orientation"],"HOLD")
    def test_original_source_hash_drift_fails(self):
        x=fixture();x["full_sha256"]="f"*64
        with self.assertRaises(ValueError):verify(x)
    def test_claimed_patient_labels_fails(self):
        x=fixture();x["metadata_verified_canonical_sample_ids"]=True
        with self.assertRaises(ValueError):verify(x)
    def test_feature_count_guard(self):
        x=fixture();x["layers"][0]["common_features_n"]=19767
        with self.assertRaises(ValueError):verify(x)
    def test_source_export_guard(self):
        x=fixture();x["layers"][1]["sample_feature_identifiers_exported"]=True
        with self.assertRaises(ValueError):verify(x)
    def test_not_forced_to_show_zscore_exactness(self):
        x=fixture();x["layers"][0]["heldout_median_abs_z_residual_after"]=0.5
        x["layers"][0]["heldout_p95_abs_z_residual_after"]=1.0
        x["layers"][0]["heldout_max_abs_z_residual_after"]=2.0
        x["layers"][0]["heldout_cell_fraction_abs_lt_1e6_after"]=0
        self.assertEqual(verify(x)["V03"],"PARTIAL")
if __name__=="__main__":unittest.main()
