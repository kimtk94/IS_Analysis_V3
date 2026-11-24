"""Stage29 source metadata, gene correlation, and scientific claims do not conflate."""
import pathlib,sys,unittest,copy
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from verify_ng_stage29_parallel import validate,COMMIT,FULL,FULLZ,XLSX
def published():
    return {"schema":"hcc-stage29-original-fig2b-rho-extraction-v1",
        "source_xlsx_sha256":XLSX,"source_sheet":"Fig2b",
        "physical_source_worksheet_index":8,"published_ENSG_n":5481,
        "published_rho_numeric_n":5481,"source_spearman_median":0.287059,
        "source_gene_and_per_gene_effect_values_exported":False,
        "original_Ng_2022_quantitation_reproduced":False}
def effects():
    cohort=[]
    for k in ("naive_50_biopsies","patient_mean_48","singleton_46"):
        cohort.append({
          "cohort":k,"matched_effect_gene_n":5253,
          "source_median_rho":0.284,"published_same_gene_median_rho":0.287,
          "between_studies_gene_rho_pearson":0.7,"between_studies_gene_rho_spearman":0.6,
          "median_absolute_gene_rho_difference":0.1,
          "median_signed_gene_rho_difference":-0.002,
          "gene_rho_difference_q10":-0.4,"gene_rho_difference_q90":0.4,
          "gene_effect_sign_concordance_fraction":0.7,
          "gene_effect_sign_discordance_fraction":0.3,
          "fraction_abs_gene_rho_difference_gt_0_2":0.2})
    return {"schema":"hcc-stage29-original-fig2b-per-gene-effect-benchmark-v1",
      "source_git_commit":COMMIT,"source_full_RDS_sha256":FULL,
      "official_original_published_gene_n":5481,
      "current_source_strict_RNA_protein_pair_n":5343,
      "matching_original_published_ENSG_n":5253,"source_biopsy_n":50,
      "independent_patient_n":48,"singleton_patient_n":46,
      "cohort_metrics":cohort,
      "original_2022_quantitation_recreated":False,
      "original_2022_2026_gene_annotation_identical":False,
      "source_original_patient_tumor_normal_and_batch_controls_verified":False,
      "independent_cohort_replication_completed":False,
      "per_gene_rho_and_patient_identifiers_exported":False}
def filtered():
    return {"schema":"hcc-stage29-source-filter-588-transcripts-v1",
       "source_commit":COMMIT,"full_sha256":FULL,"fullz_sha256":FULLZ,
       "full_RNA_features_n":20356,"fullz_RNA_features_n":19768,
       "retained_RNA_features_n":19768,"omitted_RNA_features_n":588,
       "fullz_RNA_missing_from_full_n":0,
       "dropped_in_bottom_588_by_SD_n":584,
       "dropped_in_bottom_588_by_rowMean_n":470,
       "SD_strict_threshold_partition_possible":False,
       "mean_strict_threshold_partition_possible":False,
       "original_author_filter_rule_source_verified":False,
       "filtered_gene_identifiers_or_source_expression_exported":False,
       "clinical_significance_evaluated":False}
def labels():
    return {"schema":"hcc-stage29-patient-group-label-orientation-v1",
      "source_full_sha256":FULL,"source_fullz_sha256":FULLZ,
      "source_original_supplement_patient_groups_verified":True,
      "official_source_biopsies_n":50,"official_unique_patients_n":48,
      "original_label_sets_identical_across_all_six_T_P_PP_layers":True,
      "original_full_T_P_PP_same_order":True,
      "original_fullz_T_P_PP_same_order":True,
      "full_T_vs_fullz_T_sample_name_order_equal":False,
      "crossmodal_inferred_permutation_identical":True,
      "same_biopsy_label_position_best_match_n":2,
      "same_official_patient_group_under_candidate_crossfile_mapping_n":2,
      "different_official_patient_group_under_candidate_crossfile_mapping_n":48,
      "within_same_official_patient_group_but_different_biopsy_position_n":0,
      "authoritative_orientation_of_source_biopsy_labels_verified":False,
      "original_preprocessing_creator_reproduced":False,
      "source_sample_labels_or_patient_IDs_exported":False,
      "biological_effects_inferred":False}
class Contract(unittest.TestCase):
    def test_all_synthetic_sources_technical_only(self):
        for name,data in (("fig2",published()),("effects",effects()),
                          ("filter",filtered()),("labels",labels())):
            v=validate(name,data)
            self.assertEqual(v["V03"],"PARTIAL")
            self.assertEqual(v["V07"],"PARTIAL")
            self.assertFalse(v["independent_hcc_biological_replication"])
    def test_paper_column_provenance(self):
        d=published();d["source_sheet"]="SuppFig13g"
        with self.assertRaises(ValueError):validate("fig2",d)
    def test_paper_rho_values_private(self):
        d=published();d["source_gene_and_per_gene_effect_values_exported"]=True
        with self.assertRaises(ValueError):validate("fig2",d)
    def test_same_patient_sensitivity_required(self):
        d=effects();d["cohort_metrics"].pop()
        with self.assertRaises(ValueError):validate("effects",d)
    def test_external_cohort_claim_forbidden(self):
        d=effects();d["independent_cohort_replication_completed"]=True
        with self.assertRaises(ValueError):validate("effects",d)
    def test_588_cardinality_required(self):
        d=filtered();d["omitted_RNA_features_n"]=587
        with self.assertRaises(ValueError):validate("filter",d)
    def test_premature_label_orientation_forbidden(self):
        d=labels();d["authoritative_orientation_of_source_biopsy_labels_verified"]=True
        with self.assertRaises(ValueError):validate("labels",d)
    def test_patient_group_count_inconsistency_rejected(self):
        d=labels();d["different_official_patient_group_under_candidate_crossfile_mapping_n"]=47
        with self.assertRaises(ValueError):validate("labels",d)
if __name__=="__main__":unittest.main()
