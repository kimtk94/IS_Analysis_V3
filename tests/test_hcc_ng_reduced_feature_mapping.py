"""Synthetic bnInfo reduced mapping provenance contract tests."""
import pathlib,sys,unittest
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"/"hcc"))
from verify_ng_reduced_feature_mapping import check,SRC

def doc():
    return {
      "status":"SOURCE_REDUCED_GENE_CROSSWALK_PRESENT",
      "source_commit":SRC,"HCCdata_sha256":"a"*64,"HCCinfo_sha256":"b"*64,
      "source_mapping_object_class":"bnInfo",
      "samples_per_modality":50,"transcript_feature_total":188,"protein_feature_total":116,
      "T_source_feature_keys_matched_n":188,"P_source_feature_keys_matched_n":116,
      "T_unique_gene_symbols_mapped_n":187,"P_unique_gene_symbols_mapped_n":115,
      "T_P_shared_gene_symbols_n":28,
      "T_P_unambiguous_one_to_one_genes_n":23,
      "full_feature_Ensembl_protein_mapping_verified":False,
      "source_gene_mapping_independently_curated":False,
      "clinical_effects_verified":False,
      "patient_ids_feature_names_expression_values_exported":False
    }
class Tests(unittest.TestCase):
    def test_reduced_mapping_does_not_verify_full_quant_gene_mapping(self):
        r=check(doc())
        self.assertEqual(r["shared_genes_in_selected_188T_116P"],28)
        self.assertEqual(r["full_quantitative_gene_wide_mapping"],"HOLD")
    def test_no_anchors_fails(self):
        d=doc();d["T_P_shared_gene_symbols_n"]=0;d["T_P_unambiguous_one_to_one_genes_n"]=0
        with self.assertRaises(ValueError):check(d)
    def test_raw_samples_flag_fails(self):
        d=doc();d["patient_ids_feature_names_expression_values_exported"]=True
        with self.assertRaises(ValueError):check(d)
    def test_wrong_source_rejected(self):
        d=doc();d["source_commit"]="main"
        with self.assertRaises(ValueError):check(d)
if __name__=="__main__":unittest.main()
