import pathlib,sys,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage42_cdap_pdc_numeric_prefix import compare,finite

class Stage42Tests(unittest.TestCase):
    def case(self,permutation=False):
        phead=["Gene/Aliquot"]+[f"PDC{j}" for j in range(330)]
        names=["GeneA","GeneB","GeneC","GeneD","GeneE"]
        pdata=[[name]+[str((idx+1)*100+j) for j in range(330)]
               for idx,name in enumerate(names)]
        shead=["Gene"]+[f"S{j} Log Ratio" for j in range(330)]
        order=list(range(330))
        if permutation:order=order[100:]+order[:100]
        sdata=[[name]+[str((idx+1)*100+j) for j in order]
               for idx,name in enumerate(names)]
        return phead,pdata,shead,sdata
    def test_numeric_permutation_source_candidate_not_biology(self):
        x=compare(*self.case(True))
        self.assertEqual(x["unique_exact_feature_row_pairs_n"],5)
        self.assertEqual(x["unique_exact_numeric_fingerprint_candidates_n"],330)
        self.assertFalse(x["numerical_matches_prove_official_aliquot_identity"])
        self.assertNotIn("GeneA",str(x))
        self.assertNotIn("S100",str(x))
    def test_no_exact_feature_names_prevents_matching(self):
        a,b,c,d=self.case()
        for row in d:row[0]="X_"+row[0]
        x=compare(a,b,c,d)
        self.assertEqual(x["unique_exact_feature_row_pairs_n"],0)
        self.assertEqual(x["unique_exact_numeric_fingerprint_candidates_n"],0)
    def test_duplicate_numeric_profiles_fail_unique(self):
        a,b,c,d=self.case()
        for row in b:row[1]=row[2]
        for row in d:row[1]=row[2]
        x=compare(a,b,c,d)
        self.assertLess(x["unique_exact_numeric_fingerprint_candidates_n"],330)
    def test_missing_numeric_values_hold(self):
        a,b,c,d=self.case()
        for row in d:row[1]=""
        x=compare(a,b,c,d)
        self.assertEqual(x["CDAP_columns_with_complete_numeric_fingerprint_n"],329)
    def test_non_finite_is_not_a_numeric_fingerprint(self):
        self.assertIsNone(finite("nan"));self.assertIsNone(finite("Inf"))
        self.assertIsNone(finite(""));self.assertEqual(finite("1.5"),1.5)
    def test_original_330_shape_required(self):
        a,b,c,d=self.case()
        with self.assertRaises(ValueError):compare(a[:-1],b,c,d)

if __name__=="__main__":unittest.main()
