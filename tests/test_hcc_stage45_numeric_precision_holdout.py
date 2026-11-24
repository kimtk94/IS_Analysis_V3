import pathlib,sys,unittest
import numpy as np
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage45_numeric_precision_holdout import heldout_parity,digits_of_decimal,verify

class Stage45Tests(unittest.TestCase):
    def test_rounded_permuted_original_quantitative_source(self):
        rng=np.random.default_rng(45)
        original=rng.normal(size=(80,330))
        perm=rng.permutation(330)
        rounded=np.round(original,4)[:,perm]
        out=heldout_parity(original,rounded)
        self.assertEqual(out["independent_inferred_mapping_agreement_n"],330)
        self.assertEqual(out["raw_numeric_heldout_comparisons_n"],40*330)
        self.assertTrue(out["all_heldout_prefix_cells_agree_at_4dp_rounding_scale"])
        self.assertLessEqual(out["raw_numeric_absolute_difference_max"],.000051)
        self.assertGreater(out["negative_rotated_column_control_absolute_difference_median"],.1)
        self.assertFalse(out["official_PDC_aliquot_ID_to_columns_verified"])
        self.assertFalse(out["any_source_cell_or_column_mapping_exported"])
    def test_independent_random_source_rejected(self):
        rng=np.random.default_rng(100)
        a=rng.normal(size=(60,330))
        b=rng.normal(size=(60,330))
        out=heldout_parity(a,b)
        self.assertFalse(out["heldout_prefix_4dp_scale_consistency_99_9pct_candidate"])
        self.assertFalse(out["all_heldout_prefix_cells_agree_at_4dp_rounding_scale"])
    def test_affine_not_rounded_source_rejected(self):
        rng=np.random.default_rng(83)
        a=rng.normal(size=(64,330))
        b=np.round(a*1.1+1.2,4)
        out=heldout_parity(a,b)
        self.assertFalse(out["all_heldout_prefix_cells_agree_at_4dp_rounding_scale"])
    def test_decimal_precision_parser(self):
        self.assertEqual(digits_of_decimal("-0.1234"),4)
        self.assertEqual(digits_of_decimal("1.20"),2)
        self.assertEqual(digits_of_decimal("1.2e-03"),1)
        self.assertEqual(digits_of_decimal("3"),0)
        self.assertIsNone(digits_of_decimal("privateID"))
    def test_full_matrix_or_patient_claim_forbidden(self):
        x={"schema":"hcc-stage45-PDC-CDAP-heldout-precision-parity-v1",
           "source":{"status":"ACCESS_HOLD","raw_sample_feature_or_quant_values_exported":False,
                     "official_sample_IDs_joined":False,"external_RNA_protein_matched":False}}
        self.assertEqual(verify(x)["G30_D"],"STOP_NO_INDEPENDENT_MATCHED_PATIENT_RNA_PROTEIN")
        x["source"]["official_sample_IDs_joined"]=True
        with self.assertRaises(ValueError):verify(x)

if __name__=="__main__":
    unittest.main()
