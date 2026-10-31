import pathlib,sys,unittest
import numpy as np
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage44_numeric_affine_heldout import evaluate,verify

class Stage44Tests(unittest.TestCase):
    def test_affine_gene_wise_transform_on_original_column_permutation(self):
        rng=np.random.default_rng(41)
        p=rng.normal(size=(80,330))
        perm=rng.permutation(330)
        slope=rng.uniform(.65,1.85,size=(80,1))
        offset=rng.normal(size=(80,1))*2.7
        c=(p*slope+offset)[:,perm]
        metrics=evaluate(p,c)
        self.assertEqual(metrics["source_column_assignment_train_vs_test_agreement_n"],330)
        self.assertGreater(metrics["mapped_feature_spearman_median"],.99999)
        self.assertGreater(metrics["affine_heldout_sample_R2_median"],.99999)
        self.assertTrue(metrics["numerical_source_transform_technical_GO_candidate"])
        self.assertLess(metrics["exact_untransformed_numeric_cell_fraction_1e_minus6"],.01)
        self.assertFalse(metrics["official_PDC_aliquot_or_case_identity_certified"])
        self.assertFalse(metrics["raw_sample_column_map_or_gene_values_exported"])
    def test_independent_dataset_does_not_look_like_mapping(self):
        rng=np.random.default_rng(1234)
        p=rng.normal(size=(70,330))
        c=rng.normal(size=(70,330))
        q=evaluate(p,c)
        self.assertFalse(q["numerical_source_transform_technical_GO_candidate"])
        self.assertLess(q["source_column_assignment_train_vs_test_agreement_n"],100)
        self.assertLess(q["affine_heldout_sample_R2_median"],.2)
    def test_invalid_source_axes_rejected(self):
        with self.assertRaises(ValueError):evaluate(np.ones((7,330)),np.ones((7,330)))
        with self.assertRaises(ValueError):evaluate(np.ones((30,330)),np.ones((30,329)))
    def test_failing_closed_on_fake_biological_claim(self):
        x={"schema":"hcc-stage44-source-featurewise-affine-holdout-v1",
           "study":{"status":"ACCESS_HOLD",
             "private_source_identifiers_exported":False,
             "PDC_aliquot_case_join_certified":False,"NODE_RNA_matched":False}}
        self.assertEqual(verify(x)["G30_D"],"STOP_NO_INDEPENDENT_SAME_BIOPSY_RNA_PROTEIN")
        x["study"]["PDC_aliquot_case_join_certified"]=True
        with self.assertRaises(ValueError):verify(x)

if __name__=="__main__":
    unittest.main()
