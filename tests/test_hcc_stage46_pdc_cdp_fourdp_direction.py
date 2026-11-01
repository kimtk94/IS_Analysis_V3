import pathlib,sys,unittest
import numpy as np
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage46_pdc_cdp_fourdp_direction import rounding_direction,validate,digit_histogram_from_original_qdata

class Stage46Tests(unittest.TestCase):
    def test_canonical_4dp_pdc_from_highprecision_original(self):
        rng=np.random.default_rng(460)
        high=rng.normal(size=(82,330))
        pdc=np.round(high,4)
        cdap=high[:,rng.permutation(330)]
        out=rounding_direction(pdc,cdap)
        self.assertEqual(out["independently_inferred_column_mappings_same_n"],330)
        self.assertEqual(out["PDC_equals_CDAP_rounded_to_4dp_fraction"],1.0)
        self.assertTrue(out["PDC_equals_CDAP_4dp_original_feature_prefix_confirmed"])
        self.assertLess(out["CDAP_equals_PDC_rounded_to_4dp_fraction"],0.1)
        self.assertFalse(out["PDC_official_aliquot_case_identifiers_verified"])
    def test_affine_source_change_cannot_pass_rounding(self):
        rng=np.random.default_rng(461)
        x=rng.normal(size=(62,330))
        out=rounding_direction(np.round(x*1.03+0.2,4),x)
        self.assertFalse(out["PDC_equals_CDAP_4dp_original_feature_prefix_confirmed"])
        self.assertLess(out["PDC_equals_CDAP_rounded_to_4dp_fraction"],.1)
    def test_unrelated_numeric_source_does_not_pass(self):
        rng=np.random.default_rng(462)
        x=rng.normal(size=(60,330))
        y=rng.normal(size=(60,330))
        out=rounding_direction(x,y)
        self.assertFalse(out["PDC_equals_CDAP_4dp_original_feature_prefix_confirmed"])
    def test_declared_decimals_are_aggregated_only(self):
        x=digit_histogram_from_original_qdata([["Gene","0.12","0.1234","-1.2345"]])
        self.assertEqual(x,{"2":1,"4":2})
    def test_scientific_gate_must_remain_hold(self):
        q={"schema":"hcc-stage46-PDC-CDAP-fourdp-direction-v1",
           "source":{"status":"SOURCE_ACCESS_HOLD",
             "original_source_sample_or_gene_ids_exported":False,
             "full_source_numeric_matrix_verified":False,
             "PDC_official_biospecimen_join_verified":False,
             "independent_RNA_protein_matched":False}}
        self.assertEqual(validate(q)["G30_D"],"STOP_NO_INDEPENDENT_PATIENT_MATCHED_RNA_PROTEIN")
        q["source"]["independent_RNA_protein_matched"]=True
        with self.assertRaises(ValueError):validate(q)

if __name__=="__main__":
    unittest.main()
