"""Primary source public expression-matrix and GIGASTROKE GWAS beta-unit rules."""
import importlib
import json
import unittest
from pathlib import Path

M=importlib.import_module("scripts.is.audit_p3a_36_source_expression_and_gwas_beta_units")


def synthetic_zenodo():
    return json.dumps({"id":4678936,"doi":"10.5281/zenodo.4678936",
        "metadata":{"title":"RNA-seq counts"},
        "files":[{"key":"GTEx_gene_counts.tsv.gz","size":85000,
                  "checksum":"md5:abc"},
                 {"key":"ROSMAP_counts.tsv.gz","size":60000},
                 {"key":"all_studies_manifest.txt","size":100}]}).encode()


def fake_qtlmap():
    return b"""  2. count_matrix - normalised phenotype matrix. The name of the first column should be phenotype_id.
  3. pheno_meta - specifies the location of your phenotype metadata file.
  4. sample_meta - specifies the location of your sample metadata file.
  5. vcf - Genotype file used for QTL analysis.
  Note that qtlmap uses the intersection of samples present in count_matrix and sample_meta for analysis.
"""


def fake_paper():
    return b"""<?xml version="1.0"?><article>
    <front><article-meta><article-id pub-id-type="doi">10.1038/s41586-022-05165-3</article-id></article-meta></front>
    <body><sec><title>Genetics</title>
      <p>We performed logistic regression within each independent cohort followed by
      an inverse-variance weighted meta-analysis of the genome-wide stroke
      summary statistics, retaining effect sizes and standard errors.</p>
    </sec></body></article>"""


class SourceTests(unittest.TestCase):
    def test_original_zenodo_public_counts_not_qtlmap_normalized_expression(self):
        r=M.parse_zenodo(synthetic_zenodo())
        self.assertEqual(r["n_public_files"],3)
        self.assertEqual(r["n_files_mentioning_study_names"],2)
        self.assertFalse(r["direct_Release7_phenotype_expression_sdY_verified"])
        self.assertFalse(r["direct_Release7_sample_meta_verified"])

    def test_unrelated_zenodo_id_fails(self):
        r=json.loads(synthetic_zenodo())
        r["id"]=98765
        r["doi"]="10.5281/zenodo.98765"
        with self.assertRaises(ValueError):
            M.parse_zenodo(json.dumps(r).encode())

    def test_zenodo_empty_file_listing_fails_closed(self):
        r=json.loads(synthetic_zenodo());r["files"]=[]
        with self.assertRaises(ValueError):
            M.parse_zenodo(json.dumps(r).encode())

    def test_qtlmap_studyfile_inputs_and_sample_intersection(self):
        m=M.parse_qtlmap(fake_qtlmap())
        self.assertTrue(m["original_required_fields_verified"]["VCF"])
        self.assertTrue(m["study_expression_and_sample_intersection_required"])
        self.assertFalse(m["does_not_establish_exact_target_source_expression_sdY"] is False)

    def test_qtlmap_without_sample_intersection_is_not_model_N_evidence(self):
        with self.assertRaises(ValueError):
            M.parse_qtlmap(fake_qtlmap().replace(
                b"intersection of samples present in count_matrix and sample_meta",
                b"all samples"))

    def test_paper_doi_exact_2022(self):
        r=M.paper_evidence(fake_paper())
        self.assertEqual(r["verified_publication_DOI"],M.ORIGINAL_PAPER_DOI)
        self.assertFalse(r["per_exact_GWAS_release_GCST90104540_beta_log_OR_verified"])
        self.assertFalse(r["per_exact_GWAS_release_GCST90104545_beta_log_OR_verified"])

    def test_paper_wrong_doi_rejected(self):
        with self.assertRaises(ValueError):
            M.paper_evidence(fake_paper().replace(b"s41586-022-05165-3",
                                                b"s41586-022-00000-0"))

    def test_original_paper_logistic_does_not_auto_prove_release_beta_unit(self):
        r=M.paper_evidence(fake_paper().replace(
            b"summary statistics, retaining effect sizes and standard errors",
            b"summary statistics on the log odds scale"))
        self.assertTrue(r["mentions_explicit_log_odds_around_method"])
        self.assertFalse(r["per_exact_GWAS_release_GCST90104540_beta_log_OR_verified"])

    def test_tissue_identities_and_allele_N_are_study_specific(self):
        self.assertEqual(M.EXPECTED_QTLS["QTD000131"],
                         ("QTS000015","GTEx","artery_aorta",387,10))
        self.assertEqual(M.EXPECTED_QTLS["QTD000434"],
                         ("QTS000025","ROSMAP","brain_naive",560,15))

    def test_avoid_arbitrary_sdY_and_N_eff_from_source_names(self):
        r=M.parse_zenodo(synthetic_zenodo())
        self.assertTrue(r["public_raw_or_processed_counts_cannot_prove_same_qtlmap_matrix"])
        self.assertFalse(r["direct_Release7_phenotype_expression_sdY_verified"])


if __name__=="__main__":
    unittest.main()
