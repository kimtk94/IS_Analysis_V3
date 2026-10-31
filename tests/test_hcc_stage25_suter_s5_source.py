"""Stage25 Suter S5 provenance source-only guards, no live network."""
import pathlib
import sys
import unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_suter_s5_preprocessing_source import summarize_text,fetch_pdf,MAX_BYTES,DOI

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
    def test_unidentifiable_text_is_fail_closed(self):
        with self.assertRaises(ValueError):summarize_text("placeholder")
        with self.assertRaises(ValueError):summarize_text("this is not a method "*50)
    def test_only_bounded_pdfs_allowed(self):
        self.assertTrue(DOI.endswith(".s013"))
        self.assertLess(MAX_BYTES,5_000_000)
if __name__=="__main__":unittest.main()
