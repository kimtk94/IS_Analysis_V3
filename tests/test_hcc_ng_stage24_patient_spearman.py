"""Stage24 privacy, patient grouping and noninferential QC unit tests."""
import copy
import pathlib
import sys
import unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from probe_ng_full_uniprot_ensembl_crosswalk import strict_gene_protein_pairs
from build_ng_stage24_patient_index import assign_patient_indexes
from validate_ng_stage24_spearman import verify

def fixtures():
    pairs=5343
    stats={"finite_feature_n":pairs,"undefined_constant_feature_n":0,
           "median":0.115,"q10":-0.401,"q90":0.580,
           "fraction_positive":0.623,"fraction_abs_ge_0_5":0.203}
    mapping={"audit_status":"PASS_AGGREGATE_COUNTS_ONLY","strict_one_gene_one_protein_genes_n":pairs}
    corr={
      "schema":"hcc-ng-stage24-descriptive-spearman-v1",
      "source_RDS_sha256":"461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a",
      "mapping_scope":"source_RNA_x_current_UniProt_reviewed_strict_one_to_one",
      "crosswalk_protein_gene_pair_n":pairs,"source_biopsy_n":50,
      "independent_patient_n":48,"patient_repeated_biopsy_groups_n":2,
      "singleton_only_patient_n":46,
      "patient_mean_48":stats.copy(),"singleton_only":{**stats,"median":0.100},
      "biopsy_50_nonindependent_sensitivity":stats.copy(),
      "median_rho_patient_mean_minus_singleton_only":0.015,
      "full_gene_universe_original_Ng_reproduced":False,
      "original_Ng_5481_gene_correlation_reproduced":False,
      "multiple_testing_correction_or_gene_level_pvalues_calculated":False,
      "direct_etiology_or_recurrence_tested":False,
      "scientific_significance_or_causality_proven":False,
      "feature_patient_identifiers_or_expression_exported":False,
    }
    return mapping,corr

class Tests(unittest.TestCase):
    def test_valid_still_no_scientific_claim(self):
        mapping,corr=fixtures()
        review=verify(mapping,corr)
        self.assertEqual(review["V07"],"PARTIAL_DESCRIPTIVE_NO_SIGNIFICANCE")
        self.assertEqual(review["original_Ng_correlation_replication"],"NOT_REPRODUCED")
    def test_original_paper_claim_rejected(self):
        mapping,corr=fixtures();corr["original_Ng_5481_gene_correlation_reproduced"]=True
        with self.assertRaises(ValueError):verify(mapping,corr)
    def test_wrong_denominator_rejected(self):
        mapping,corr=fixtures();corr["crosswalk_protein_gene_pair_n"]=5437
        with self.assertRaises(ValueError):verify(mapping,corr)
    def test_false_independence_rejected(self):
        mapping,corr=fixtures();corr["independent_patient_n"]=50
        with self.assertRaises(ValueError):verify(mapping,corr)
    def test_correlations_out_of_range_rejected(self):
        mapping,corr=fixtures();corr["patient_mean_48"]["q90"]=1.02
        with self.assertRaises(ValueError):verify(mapping,corr)
    def test_export_flag_rejected(self):
        mapping,corr=fixtures();corr["feature_patient_identifiers_or_expression_exported"]=True
        with self.assertRaises(ValueError):verify(mapping,corr)
    def test_unique_crosswalk_omits_multigene_and_duplicate_gene(self):
        t={"ENSG00000111111","ENSG00000222222","ENSG00000333333"}
        p=["P04637","Q14534","P12345","P67890"]
        refs={"P04637":{"ENSG00000111111"},
              "Q14534":{"ENSG00000111111"},
              "P12345":{"ENSG00000222222","ENSG00000333333"},
              "P67890":{"ENSG00000333333"}}
        self.assertEqual(strict_gene_protein_pairs(t,p,refs),
                         [("ENSG00000333333","P67890")])
    def test_patient_index_preserves_order_and_two_repeats(self):
        names=["biopsy"+str(i) for i in range(50)]
        columns={"tumorbiopsyid":0,"patientid":1,"rnaseq":2,"proteomics":3}
        rows=[{0:n,1:"pat"+str(min(i,47)),2:"Y",3:"Y"} for i,n in enumerate(names)]
        idx=assign_patient_indexes(names,rows,columns)
        self.assertEqual(len(idx),50)
        self.assertEqual(len(set(idx)),48)
        self.assertEqual(idx[-1],idx[-2])
    def test_mismatched_biopsy_fails(self):
        names=["biopsy"+str(i) for i in range(50)]
        cols={"tumorbiopsyid":0,"patientid":1,"rnaseq":2,"proteomics":3}
        rows=[{0:n,1:"pat"+str(min(i,47)),2:"Y",3:"Y"} for i,n in enumerate(names)]
        rows.pop()
        with self.assertRaises(ValueError):assign_patient_indexes(names,rows,cols)

if __name__=="__main__":unittest.main()
