import unittest
from workflow.is_evidence_guardrails import canonical_variant, coloc_class, ld_audit, functional_priority

class ISEvidenceGuardrailsTests(unittest.TestCase):
    def setUp(self):
        self.row = dict(build="GRCh37", chr="12", pos="112241766", ref="G", alt="A",
                        effect_allele="A", other_allele="G", beta="-0.1",
                        se="0.02", p="1e-8", eaf="0.25", ancestry="EAS")
    def test_functional_requires_evidence(self):
        r = dict(gene="FGF5", mechanism_branch="regulatory", gwas_qc="PASS",
                 coloc_class="STRONG_COLOC", ld_match="MATCHED",
                 functional_context="vascular", functional_support="YES")
        self.assertEqual(functional_priority(r), "TIER_1_CONVERGENT")
        self.assertEqual(functional_priority(dict(r, ld_match="BLOCK")), "QC_UNRESOLVED")
        self.assertEqual(functional_priority(dict(r, functional_support="PENDING")), "TIER_2_FOLLOWUP")
    def test_variant_valid(self):
        r = canonical_variant(self.row)
        self.assertEqual(r["variant_key"], "GRCH37:12:112241766:G:A")
        self.assertEqual(r["qc"], "PASS")
        self.assertEqual(r["effect_orientation"], "ALT")
    def test_missing_allele_block(self):
        self.assertEqual(canonical_variant(dict(self.row, effect_allele="C"))["qc"], "BLOCK")
    def test_palindromic_requires_review(self):
        self.assertIn("PALINDROMIC_REVIEW", canonical_variant(dict(self.row, ref="A", alt="T", effect_allele="T", other_allele="A"))["qc_flags"])
    def test_ld_mismatch(self):
        self.assertEqual(ld_audit(dict(gwas_ancestry="EAS", ld_ancestry="EUR", ld_cohort="1000G", ld_file_verified="YES")), "BLOCK")
    def test_coloc_abf_never_strong(self):
        self.assertEqual(coloc_class(dict(method="ABF", ld_match="MATCHED", harmonization_qc="PASS", pp_h4="0.99", pp_h3="0.001")), "ABF_ONLY")
    def test_signal_strong(self):
        r = dict(method="SUSIE_SIGNAL", ld_match="MATCHED", harmonization_qc="PASS", pp_h4="0.9", pp_h3="0.02", signal_pair_valid="YES")
        self.assertEqual(coloc_class(r), "STRONG_COLOC")
    def test_signal_requires_pair(self):
        r = dict(method="SUSIE_SIGNAL", ld_match="MATCHED", harmonization_qc="PASS", pp_h4="0.9", pp_h3="0.02")
        self.assertEqual(coloc_class(r), "MULTISIGNAL_UNRESOLVED")
    def test_ld_requires_verified_reference(self):
        self.assertEqual(ld_audit(dict(gwas_ancestry="EAS", ld_ancestry="EAS", ld_cohort="JPT", ld_file_verified="NO")), "BLOCK")

if __name__ == "__main__": unittest.main()
