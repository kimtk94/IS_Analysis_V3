"""Full Ng/Suter RNA/proteome namespace audit never licenses a biological claim."""
import pathlib,sys,unittest
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"/"hcc"))
from verify_ng_full_feature_namespace import COMMIT,check
def fixture():
    return {
      "schema":"hcc-ng-full-namespace-v1","source_commit":COMMIT,
      "full_sha256":"a"*64,"info_sha256":"b"*64,
      "n_full_RNA":20356,"n_full_protein":5585,
      "RNA_ENSG_like_n":20356,"protein_UniProt_like_n":5585,
      "protein_piped_UniProt_like_n":0,"RNA_unique_feature_labels_n":20356,
      "protein_unique_feature_labels_n":5585,
      "raw_feature_names_direct_overlap_n":0,"RNA_ENSG_version_suffix_n":0,
      "selected_RNA_lookup_matches_full_n":188,
      "selected_protein_lookup_matches_full_n":116,
      "Ng2022_paper_comparable_genes_benchmark":5490,
      "Ng2022_paper_mRNA_protein_correlation_benchmark":5481,
      "full_ensembl_uniprot_mapping_verified":False,
      "same_biopsy_rna_protein_abundance_effect_tested":False,
      "feature_or_patient_ids_exported":False
    }
class Tests(unittest.TestCase):
    def test_paper_benchmark_does_not_promote_source_matrix_mapping(self):
        r=check(fixture())
        self.assertEqual(r["full_gene_protein_join"],"HOLD")
        self.assertEqual(r["original_paper_benchmark_comparable_n"],5490)
    def test_unpinned_source_rejected(self):
        v=fixture();v["source_commit"]="main"
        with self.assertRaises(ValueError):check(v)
    def test_unverified_quantitative_effect_rejected(self):
        v=fixture();v["same_biopsy_rna_protein_abundance_effect_tested"]=True
        with self.assertRaises(ValueError):check(v)
    def test_omitted_ensembl_labels_rejected(self):
        v=fixture();v["RNA_ENSG_like_n"]=0
        with self.assertRaises(ValueError):check(v)
if __name__=="__main__":unittest.main()
