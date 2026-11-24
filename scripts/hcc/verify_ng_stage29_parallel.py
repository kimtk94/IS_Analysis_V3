#!/usr/bin/env python3
"""Stage29 fail-closed source and scientific-claim checks; aggregate-only."""
import argparse,json,math,pathlib
FULL="461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a"
FULLZ="9c1e4e05472a0aa477b81aa7b4c4e1b27809df9a2b71c1f57e247bb031263007"
XLSX="16376438aa16a09a13db417a4ae1f4cc16b87ec2c47630aeae651004501be420"
COMMIT="438f52271d0a86a6701baa7a96e80c516ba575c2"
def fail_false(d,keys):
    for k in keys:
        if d.get(k) is not False:raise ValueError("disallowed claim or leak: "+k)
def validate(kind,d):
    if kind=="fig2":
        if (d.get("schema")!="hcc-stage29-original-fig2b-rho-extraction-v1" or
            d.get("source_xlsx_sha256")!=XLSX or d.get("source_sheet")!="Fig2b" or
            d.get("physical_source_worksheet_index")!=8 or
            d.get("published_ENSG_n")!=5481 or d.get("published_rho_numeric_n")!=5481 or
            abs(d.get("source_spearman_median",99)-0.287059)>1e-5):
            raise ValueError("Fig2b XLSX provenance invalid")
        fail_false(d,("source_gene_and_per_gene_effect_values_exported","original_Ng_2022_quantitation_reproduced"))
    elif kind=="effects":
        if (d.get("schema")!="hcc-stage29-original-fig2b-per-gene-effect-benchmark-v1" or
            d.get("source_git_commit")!=COMMIT or d.get("source_full_RDS_sha256")!=FULL or
            d.get("official_original_published_gene_n")!=5481 or
            d.get("current_source_strict_RNA_protein_pair_n")!=5343 or
            d.get("matching_original_published_ENSG_n")!=5253 or
            d.get("source_biopsy_n")!=50 or d.get("independent_patient_n")!=48 or
            d.get("singleton_patient_n")!=46):raise ValueError("source effect cohort mismatch")
        by={x.get("cohort"):x for x in d.get("cohort_metrics",[])}
        if set(by)!={"naive_50_biopsies","patient_mean_48","singleton_46"}:
            raise ValueError("three sensitivity cohorts required")
        for x in by.values():
            if x.get("matched_effect_gene_n")!=5253:raise ValueError("invalid original gene intersection")
            for k in ("source_median_rho","published_same_gene_median_rho","between_studies_gene_rho_pearson",
                "between_studies_gene_rho_spearman","median_absolute_gene_rho_difference",
                "median_signed_gene_rho_difference","gene_rho_difference_q10","gene_rho_difference_q90"):
                y=x.get(k)
                if type(y) not in (int,float) or not math.isfinite(y) or abs(y)>2:
                    raise ValueError("gene effect invalid "+k)
            for k in ("gene_effect_sign_concordance_fraction",
                      "gene_effect_sign_discordance_fraction",
                      "fraction_abs_gene_rho_difference_gt_0_2"):
                y=x.get(k)
                if type(y) not in (int,float) or not math.isfinite(y) or not 0<=y<=1:
                    raise ValueError("gene effect fraction invalid "+k)
            if abs(x["gene_effect_sign_concordance_fraction"]+
                   x["gene_effect_sign_discordance_fraction"]-1)>1e-6:
                raise ValueError("gene effect concordance sum invalid")
        fail_false(d,("original_2022_quantitation_recreated",
                      "original_2022_2026_gene_annotation_identical",
                      "source_original_patient_tumor_normal_and_batch_controls_verified",
                      "independent_cohort_replication_completed",
                      "per_gene_rho_and_patient_identifiers_exported"))
    elif kind=="filter":
        if (d.get("schema")!="hcc-stage29-source-filter-588-transcripts-v1" or
            d.get("source_commit")!=COMMIT or d.get("full_sha256")!=FULL or
            d.get("fullz_sha256")!=FULLZ or d.get("full_RNA_features_n")!=20356 or
            d.get("fullz_RNA_features_n")!=19768 or d.get("retained_RNA_features_n")!=19768 or
            d.get("omitted_RNA_features_n")!=588 or d.get("fullz_RNA_missing_from_full_n")!=0):
            raise ValueError("588 excluded RNA source mismatch")
        for k in ("dropped_in_bottom_588_by_SD_n","dropped_in_bottom_588_by_rowMean_n"):
            n=d.get(k)
            if type(n) is not int or not 0<=n<=588:raise ValueError("feature filter count invalid")
        for k in ("SD_strict_threshold_partition_possible","mean_strict_threshold_partition_possible"):
            if type(d.get(k)) is not bool:raise ValueError("feature separation flag missing")
        fail_false(d,("original_author_filter_rule_source_verified",
                      "filtered_gene_identifiers_or_source_expression_exported",
                      "clinical_significance_evaluated"))
    elif kind=="labels":
        if (d.get("schema")!="hcc-stage29-patient-group-label-orientation-v1" or
            d.get("source_full_sha256")!=FULL or d.get("source_fullz_sha256")!=FULLZ or
            d.get("official_source_biopsies_n")!=50 or d.get("official_unique_patients_n")!=48 or
            d.get("source_original_supplement_patient_groups_verified") is not True or
            d.get("original_label_sets_identical_across_all_six_T_P_PP_layers") is not True or
            d.get("crossmodal_inferred_permutation_identical") is not True):
            raise ValueError("source labels not anchored")
        for key in ("original_full_T_P_PP_same_order",
                    "original_fullz_T_P_PP_same_order",
                    "full_T_vs_fullz_T_sample_name_order_equal"):
            if type(d.get(key)) is not bool:raise ValueError("source array order flag absent: "+key)
        n=d.get("same_biopsy_label_position_best_match_n")
        a=d.get("same_official_patient_group_under_candidate_crossfile_mapping_n")
        b=d.get("different_official_patient_group_under_candidate_crossfile_mapping_n")
        c=d.get("within_same_official_patient_group_but_different_biopsy_position_n")
        if any(type(k) is not int or not 0<=k<=50 for k in (n,a,b,c)) or a+b!=50 or c>a:
            raise ValueError("sample match arithmetic invalid")
        fail_false(d,("authoritative_orientation_of_source_biopsy_labels_verified",
                      "original_preprocessing_creator_reproduced",
                      "source_sample_labels_or_patient_IDs_exported",
                      "biological_effects_inferred"))
    else:raise ValueError("unknown kind")
    return {"schema":"hcc-stage29-postflight-v1","kind":kind,
       "source_based_technical_QC":"PASS","V03":"PARTIAL","V07":"PARTIAL",
       "independent_hcc_biological_replication":False,
       "original_quantitation_reconstructed":False,
       "patient_identifier_or_gene_effect_leakage":False}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--kind",required=True,choices=("fig2","effects","filter","labels"))
    p.add_argument("--input",required=True);p.add_argument("--out",required=True)
    a=p.parse_args();d=json.loads(pathlib.Path(a.input).read_text(encoding="utf-8"))
    o=pathlib.Path(a.out);o.parent.mkdir(parents=True,exist_ok=True)
    x=validate(a.kind,d);o.write_text(json.dumps(x,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(x,indent=2))
if __name__=="__main__":main()
