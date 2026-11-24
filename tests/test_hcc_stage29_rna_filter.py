"""Stage29 RNA 588-feature source filter scientific/privacy constraints."""
import pathlib,sys,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from verify_ng_stage29_rna_filter import verify,FULL,FULLZ,COMMIT
def fixture(zero=True):
    def stats(n,zero_count):
        return {"n":n,"zero_sd_n":zero_count,
          "sd_lt_1e_minus12_n":zero_count,"sd_lt_1e_minus6_n":zero_count,
          "sd_min":0.0,"sd_p10":0.0,"sd_median":0.0,"sd_p90":0.0,"sd_max":1.0}
    return {"schema":"hcc-stage29-rna-fullz-exclusion-qc-v1",
      "source_commit":COMMIT,"full_sha256":FULL,"fullz_sha256":FULLZ,
      "rna_full_feature_n":20356,"rna_fullz_feature_n":19768,
      "strict_matching_feature_n":19768,"rna_full_only_feature_n":588,
      "feature_sample_n":50,
      "excluded_rna_feature_sd":stats(588,588 if zero else 500),
      "retained_rna_feature_sd":stats(19768,0),
      "rna_full_nonzero_sd_feature_n":19768 if zero else 19856,
      "rna_full_gt_1e_minus12_sd_feature_n":19768 if zero else 19856,
      "all_588_excluded_features_exactly_constant":zero,
      "retained_constant_feature_n":0,
      "feature_selection_author_code_verified":False,
      "per_gene_or_patient_values_exported":False,
      "source_rds_files_overwritten":False}
class TestRNAFilter(unittest.TestCase):
    def test_technical_only(self):
        out=verify(fixture())
        self.assertTrue(out["exact_zero_SD_filter_equivalence_supported"])
        self.assertEqual(out["V03"],"PARTIAL")
    def test_not_forcing_constant_filter(self):
        out=verify(fixture(False))
        self.assertFalse(out["exact_zero_SD_filter_equivalence_supported"])
    def test_forged_reproduction(self):
        x=fixture();x["feature_selection_author_code_verified"]=True
        with self.assertRaises(ValueError):verify(x)
    def test_wrong_sha(self):
        x=fixture();x["fullz_sha256"]="1"*64
        with self.assertRaises(ValueError):verify(x)
    def test_wrong_cardinality(self):
        x=fixture();x["rna_full_only_feature_n"]=587
        with self.assertRaises(ValueError):verify(x)
    def test_wrong_consistency(self):
        x=fixture();x["excluded_rna_feature_sd"]["zero_sd_n"]=300
        with self.assertRaises(ValueError):verify(x)
if __name__=="__main__":unittest.main()
