"""Synthetic source controls for public Ng full UniProt crosswalk."""
import pathlib,sys,unittest
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"/"hcc"))
from probe_ng_full_uniprot_ensembl_crosswalk import (
    canonical_protein,parse_uniprot_tsv,summarize,COMMIT
)
class Tests(unittest.TestCase):
    def test_reference_tsv_gives_Ensembl_gene_not_transcript(self):
        tsv=("Entry\tEnsembl\n"+
             "P04637\tENST00000269305; ENSP00000269305; ENSG00000141510\n"+
             "Q14534\tENST00000000001; ENSG00000000002\n").encode()
        refs,n=parse_uniprot_tsv(tsv)
        self.assertEqual(n,2)
        self.assertEqual(refs["P04637"],{"ENSG00000141510"})
        self.assertEqual(refs["Q14534"],{"ENSG00000000002"})
    def test_no_gene_ids_get_no_mapping(self):
        data=b"Entry\tEnsembl\nP04637\tENST00000269305\n"
        refs,n=parse_uniprot_tsv(data)
        r=summarize({"ENSG00000141510"},["P04637"],refs,n,"a"*64)
        self.assertEqual(r["protein_features_with_Ensembl_gene_in_source_RNA_n"],0)
        self.assertFalse(r["clinical_effect_or_pathway_significance_tested"])
    def test_multi_gene_protein_count_separate(self):
        refs={"P04637":{"ENSG00000141510"},
              "Q14534":{"ENSG00000141510","ENSG00000000002"}}
        r=summarize({"ENSG00000141510","ENSG00000000002"},
                    ["P04637","Q14534"],refs,2,"b"*64)
        self.assertEqual(r["distinct_Ensembl_genes_corresponding_to_proteins_n"],2)
        self.assertEqual(r["protein_features_with_unique_Ensembl_gene_n"],1)
        self.assertEqual(r["multi_gene_protein_feature_count"],1)
        self.assertFalse(r["exact_original_paper_2022_mapping_reproduced"])
        self.assertNotIn("P04637",str(r))
    def test_uniprot_accession_not_gene_symbol_guess(self):
        self.assertEqual(canonical_protein("sp|P04637|P53_HUMAN"),"P04637")
        self.assertEqual(canonical_protein("Q14534-2"),"Q14534")
        self.assertIsNone(canonical_protein("NOT_A_UNIPROT_ACCESSION"))
    def test_bad_TSV_schema_fails_closed(self):
        with self.assertRaises(ValueError):parse_uniprot_tsv(b"Protein\tGene\nFOO\tBAR\n")
if __name__=="__main__":unittest.main()
