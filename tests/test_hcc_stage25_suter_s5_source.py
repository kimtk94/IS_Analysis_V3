"""Stage25 Suter S5 provenance source-only guards, no live network."""
import pathlib
import sys
import unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_suter_s5_preprocessing_source import summarize_text,layer_method_keyword_sections,fetch_pdf,MAX_BYTES,DOI

class Tests(unittest.TestCase):
    def test_key_terms_are_mentions_not_method_proof(self):
        summary=summarize_text(
            "S5 Appendix. Data pre-processing of RNA transcriptome and protein "
            "data. Missing data imputation using KNN was performed in a batch. "
            "Normalization to healthy livers described separately."
        )
        self.assertGreater(summary["rna_or_transcriptome"],0)
        self.assertGreater(summary["protein_or_proteome"],0)
        self.assertGreater(summary["batch_or_combat"],0)
        self.assertGreater(summary["imputation"],0)
    def test_layer_extraction_does_not_transfer_pp_batch_to_rna(self):
        passage=(
            "S5 Appendix. Data pre-processing. "
            "For transcriptome data, the preprocessing steps included upper-quartile normalization to 1000 and log2 transformation. "
            "For proteome data, log2 transformation and normalization by median substraction and removing 50% missing proteins. "
            "For clustering only: imputation using impute. For differential expression analysis we used unimputed values. "
            "For phosphoproteome data, batch correction with edgeR was used. "
            "The CNA data was collected from the original paper."
        )
        a=layer_method_keyword_sections(passage)
        self.assertIn("transcriptome",a)
        self.assertIn("proteome",a)
        self.assertIn("phosphoproteome",a)
        self.assertTrue(a["transcriptome"]["upper_quartile_mentioned"])
        self.assertFalse(a["transcriptome"]["batch_correction_mentioned"])
        self.assertTrue(a["proteome"]["median_subtraction_mentioned"])
        self.assertNotIn("cna",a)
        self.assertTrue(a["proteome"]["imputation_mentioned"])
        self.assertTrue(a["proteome"]["clustering_only_qualifier_mentioned"])
        self.assertTrue(a["phosphoproteome"]["batch_correction_mentioned"])
        self.assertFalse(a["proteome"]["batch_correction_mentioned"])
    def test_unidentifiable_text_is_fail_closed(self):
        with self.assertRaises(ValueError):summarize_text("placeholder")
        with self.assertRaises(ValueError):summarize_text("this is not a method "*50)
    def test_only_bounded_pdfs_allowed(self):
        self.assertTrue(DOI.endswith(".s013"))
        self.assertLess(MAX_BYTES,5_000_000)
if __name__=="__main__":unittest.main()
