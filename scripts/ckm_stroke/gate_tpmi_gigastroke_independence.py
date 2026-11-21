#!/usr/bin/env python3
"""Strict study-level provenance gate, not a participant overlap estimator.

Sources:
  Yang et al. 2025, Nature, TPMI DOI 10.1038/s41586-025-09680-x
  Mishra et al. 2022, Nature, GIGASTROKE DOI 10.1038/s41586-022-05165-3
  EBI GWAS Catalog GCST90104545 EAS AIS discovery sample metadata.
  TPMI PheWeb Phecode 433.21.
Literature mentions separate cohorts; that alone cannot prove disjoint
individuals. The outcome discovery and replication must not be pooled.
"""
import argparse
import json
from pathlib import Path

SOURCES={
    "TPMI":"https://www.nature.com/articles/s41586-025-09680-x",
    "GIGASTROKE":"https://www.nature.com/articles/s41586-022-05165-3",
    "TPMI_PHECODE":"https://pheweb.ibms.sinica.edu.tw/pheno/433.21",
    "EAS_AIS_CATALOG":"https://www.ebi.ac.uk/gwas/rest/api/v2/studies/GCST90104545",
}

def gate(evidence):
    """Only accept documented participation nonoverlap, not a guess."""
    if not isinstance(evidence,dict):
        raise TypeError("source evidence must be an object")
    needed=("TPMI_source_cohort_ids","GIGASTROKE_EAS_AIS_discovery_cohort_ids",
            "dataset_level_overlap_report","participant_nonoverlap_provenance")
    observed={k:evidence.get(k) for k in needed}
    exact_cohort_lists=all(isinstance(observed[k],list) and bool(observed[k])
            for k in needed[:2])
    independently_supported=(exact_cohort_lists and
        isinstance(observed["dataset_level_overlap_report"],dict) and
        observed["dataset_level_overlap_report"].get("evaluated") is True and
        observed["dataset_level_overlap_report"].get("overlap_count") == 0 and
        isinstance(observed["participant_nonoverlap_provenance"],dict) and
        observed["participant_nonoverlap_provenance"].get("validated_by_dataset_custodians") is True and
        bool(observed["participant_nonoverlap_provenance"].get("source_id")))
    return {
        "schema":"TPMI-GIGASTROKE-EAS-AIS-sample-independence-evidence-v1",
        "source_references":SOURCES,
        "source_TPMI_Phecode":"433.21",
        "source_GIGASTROKE_discovery":"GCST90104545",
        "TPMI_PheWeb_reported_cases":9249,
        "TPMI_PheWeb_reported_controls":304660,
        "TPMI_Phenotype_N":313909,
        "TPMI_PLINK2_OBS_CT":204099,
        "TPMI_case_control_N_in_regression_verified":False,
        "GIGASTROKE_EAS_AIS_discovery_cases":19032,
        "GIGASTROKE_EAS_AIS_discovery_controls":237242,
        "GIGASTROKE_EAS_AIS_replication_cases":8304,
        "GIGASTROKE_EAS_AIS_replication_controls":174855,
        "discovery_and_replication_not_pooled":True,
        "study_specific_lists_present":exact_cohort_lists,
        "custodian_validated_nonoverlap":independently_supported,
        "original_study_names_are_not_nonoverlap_proof":True,
        "source_effect_snp_eligible_by_numeric_QC":1507261,
        "numeric_QC_is_not_replication":True,
        "population_sample_Neff_approved":False,
        "scientific_status":"REVIEW_NONOVERLAP_PROVEN" if independently_supported else
                             "HOLD_INDEPENDENCE_NOT_PROVEN",
        "downstream_LDSC_MR_coloc_enabled":False,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    result=gate({})
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"cohort_list_documented":result["study_specific_lists_present"],
                     "custodian_proven_no_overlap":result["custodian_validated_nonoverlap"],
                     "status":result["scientific_status"],
                     "LDSC_MR_coloc_enabled":False},ensure_ascii=False))
if __name__=="__main__":main()
