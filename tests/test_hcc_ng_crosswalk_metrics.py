"""Synthetic Stage23 crosswalk evidence gate tests (no private identifiers)."""
import copy
import pathlib
import sys
import unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from validate_ng_crosswalk_metrics import evaluate

def fixture():
    return {
      "schema":"hcc-ng-full-uniprot-ensembl-live-audit-v1",
      "reviewed_only_annotation_policy":True,
      "annotation_source":"UniProtKB reviewed Homo sapiens",
      "annotation_snapshot_utc":"2026-10-31T05:08:28+00:00",
      "annotation_tsv_sha256":"a"*64,
      "source_RNA_ENSG_n":20356,
      "source_proteome_features_n":5585,
      "source_protein_features_with_recognized_UniProt_accession_n":5530,
      "source_protein_accessions_found_in_uniprot_reviewed_n":5515,
      "protein_features_with_Ensembl_gene_in_source_RNA_n":5390,
      "protein_features_with_unique_Ensembl_gene_n":5347,
      "distinct_Ensembl_genes_corresponding_to_proteins_n":5437,
      "one_gene_one_protein_unambiguous_genes_n":5343,
      "multi_gene_protein_feature_count":43,
      "unmapped_or_noncanonical_source_protein_n":195,
      "source_paper_comparable_gene_benchmark_n":5490,
      "source_paper_correlated_gene_benchmark_n":5481,
      "absolute_deviation_from_paper_comparable_n":53,
      "uniprot_reviewed_entries_with_gene_ENSG_in_TSV_field_n":0,
      "uniprot_reviewed_entries_with_gene_ENSG_after_json_n":5485,
      "canonical_accession_isoform_collapse_excess_features_n":0,
      "accession_specific_json_enrichment":{"batches":139,"batches_completed":139},
      "exact_original_paper_2022_mapping_reproduced":False,
      "clinical_effect_or_pathway_significance_tested":False,
      "patient_sample_gene_accession_or_expression_values_exported":False,
    }

class GateTests(unittest.TestCase):
    def test_aggregate_pass_is_not_biological_pass(self):
        r=evaluate(fixture())
        self.assertEqual(r["feature_to_gene_mapping_n"],5390)
        self.assertEqual(r["distinct_ensembl_gene_union_n"],5437)
        self.assertEqual(r["paper_count_delta"],-53)
        self.assertFalse(r["original_paper_5490_gene_set_replicated"])
        self.assertFalse(r["new_biological_significance"])
    def test_multigene_feature_inconsistency_fails(self):
        d=fixture();d["multi_gene_protein_feature_count"]=42
        with self.assertRaises(ValueError):evaluate(d)
    def test_missing_raw_provenance_fails(self):
        d=fixture();d.pop("uniprot_reviewed_entries_with_gene_ENSG_in_TSV_field_n")
        with self.assertRaises(ValueError):evaluate(d)
    def test_json_incomplete_fails(self):
        d=fixture();d["accession_specific_json_enrichment"]["batches_completed"]=138
        with self.assertRaises(ValueError):evaluate(d)
    def test_false_paper_reproduction_claim_fails(self):
        d=fixture();d["exact_original_paper_2022_mapping_reproduced"]=True
        with self.assertRaises(ValueError):evaluate(d)
    def test_paper_count_cannot_be_interpreted_as_same_ids(self):
        d=fixture()
        d["distinct_Ensembl_genes_corresponding_to_proteins_n"]=5490
        d["absolute_deviation_from_paper_comparable_n"]=0
        r=evaluate(d)
        self.assertFalse(r["original_paper_5490_gene_set_replicated"])
    def test_invalid_source_protein_count_fails(self):
        d=fixture();d["unmapped_or_noncanonical_source_protein_n"]=194
        with self.assertRaises(ValueError):evaluate(d)
    def test_bad_policy_fails(self):
        d=fixture();d["reviewed_only_annotation_policy"]=False
        with self.assertRaises(ValueError):evaluate(d)
    def test_no_raw_expression_data_export(self):
        d=fixture();d["patient_sample_gene_accession_or_expression_values_exported"]=True
        with self.assertRaises(ValueError):evaluate(d)

if __name__=="__main__":
    unittest.main()
