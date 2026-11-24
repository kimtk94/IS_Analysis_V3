import pathlib,sys,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from verify_ng_stage29_hccdata_anchor import verify
def fixture():
    return {"schema":"hcc-stage29-small-HCCdata-rds-anchor-v1",
      "source_commit":"438f52271d0a86a6701baa7a96e80c516ba575c2",
      "source_HCCdata_git_blob_sha":"45538b8452a1d2e3844b596f6bd50f1086f0bf21",
      "source_full_sha256":"461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a",
      "source_fullz_sha256":"9c1e4e05472a0aa477b81aa7b4c4e1b27809df9a2b71c1f57e247bb031263007",
      "modalities":[{
          "modality":m,"small_data_feature_n":200,"small_data_sample_n":50,
          "small_data_source_numeric_matrix":True,
          "small_data_original_sample_axis_is_rows":True,
          "small_data_exact_50_biopsy_comparison_eligible":True,
          "small_full_common_feature_n":100,"small_full_common_sample_labels_n":50,
          "small_full_rho_finite_feature_n":100,
          "small_full_same_label_per_feature_rho_median":.9,
          "small_full_fraction_gene_rho_gt_095":.5,
          "small_fullz_common_feature_n":100,"small_fullz_common_sample_labels_n":50,
          "small_fullz_rho_finite_feature_n":100,
          "small_fullz_same_label_per_feature_rho_median":.1,
          "small_fullz_fraction_gene_rho_gt_095":0.,
          "biopsy_or_gene_identifiers_exported":False}
          for m in ("T","P","PP")],
      "HCCdata_is_independent_patient_cohort":False,
      "HCCdata_authoritative_biopsy_identity_proven":False,
      "source_sample_labels_corrected":False,
      "clinical_or_gene_effect_inference_performed":False}
class TestAnchor(unittest.TestCase):
    def test_technical_partial(self):
        self.assertEqual(verify(fixture())["V03"],"PARTIAL")
    def test_uninformative_names_not_forced(self):
        d=fixture()
        for x in d["modalities"]:
            for q in ("small_full","small_fullz"):
                x[q+"_common_feature_n"]=0
                x[q+"_rho_finite_feature_n"]=0
                x[q+"_same_label_per_feature_rho_median"]=None
                x[q+"_fraction_gene_rho_gt_095"]=None
        self.assertEqual(verify(d)["V03"],"PARTIAL")
    def test_canonical_label_proof_forbidden(self):
        d=fixture();d["HCCdata_authoritative_biopsy_identity_proven"]=True
        with self.assertRaises(ValueError):verify(d)
    def test_source_leak(self):
        d=fixture();d["modalities"][1]["biopsy_or_gene_identifiers_exported"]=True
        with self.assertRaises(ValueError):verify(d)
    def test_digest_guards(self):
        d=fixture();d["source_HCCdata_git_blob_sha"]="2"*40
        with self.assertRaises(ValueError):verify(d)
if __name__=="__main__":unittest.main()
