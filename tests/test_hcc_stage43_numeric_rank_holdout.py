import pathlib,sys,unittest
import numpy as np
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage43_numeric_rank_holdout import score_numeric_mapping,exact_unique_feature_rows

class Stage43Tests(unittest.TestCase):
    def test_source_gene_scaled_and_permuted_can_be_heldout_verified(self):
        rng=np.random.default_rng(39)
        x=rng.normal(size=(80,12))
        perm=rng.permutation(12)
        y=x[:,perm]*rng.uniform(0.8,3.0,size=(80,1))+rng.normal(size=(80,1))*3
        r=score_numeric_mapping(x,y)
        self.assertEqual(r["numerically_independently_inferred_column_candidate_agreement_n"],12)
        self.assertGreater(r["heldout_rank_profile_cosine_median"],0.999)
        self.assertFalse(r["source_sample_position_joins_used"])
        self.assertFalse(r["official_PDC_aliquot_case_IDs_certified"])
        self.assertFalse(r["provisional_numerical_source_bijection_candidate"])

    def test_independent_expression_has_no_perfect_holdout_agreement(self):
        rng=np.random.default_rng(211)
        x=rng.normal(size=(60,12))
        y=rng.normal(size=(60,12))
        r=score_numeric_mapping(x,y)
        self.assertLess(r["numerically_independently_inferred_column_candidate_agreement_n"],12)
        self.assertFalse(r["NODE_RNA_patient_join_certified"])

    def test_wrong_dimension_fail_closed(self):
        with self.assertRaises(ValueError):
            score_numeric_mapping(np.ones((10,9)),np.ones((10,10)))
        with self.assertRaises(ValueError):
            score_numeric_mapping(np.ones((7,10)),np.ones((7,10)))

    def test_exact_feature_label_sets_with_330_reporter_columns(self):
        ph=[["Gene/Aliquot"]+[f"Q{i}" for i in range(330)]]
        p=[[f"Gene{j}"]+[str(j*2+i) for i in range(330)] for j in range(8)]
        sh=["Gene"]+[f"Bio{i} Log Ratio" for i in range(330)]
        s=[[f"Gene{j}"]+[str(j*2+i) for i in range(330)] for j in range(8)]
        x,y,d=exact_unique_feature_rows(p,s,sh)
        self.assertEqual(x.shape,(8,330))
        self.assertEqual(y.shape,(8,330))
        self.assertEqual(d["retained_complete_feature_pairs_n"],8)
        self.assertNotIn("Gene0",str(d))

    def test_duplicate_feature_annotations_not_falsely_unique(self):
        p=[[f"G{i}"]+[float(j+i) for j in range(330)] for i in range(4)]
        s=[[f"G{i}"]+[float(j+i) for j in range(330)] for i in range(4)]
        s.append(s[0].copy())
        sh=["Gene"]+[f"Bio{k} Log Ratio" for k in range(330)]
        _,_,d=exact_unique_feature_rows(p,s,sh)
        self.assertEqual(d["retained_complete_feature_pairs_n"],3)

if __name__=="__main__":
    unittest.main()
