#!/usr/bin/env python3
"""BBJ/Koyanagi/GIGASTROKE/KCPS2 cohort-overlap provenance matrix.

This is *literature-level* cohort-membership evidence, not participant-level ID
matching. PGS BBJ evaluation is excluded from GIGASTROKE training GWAS, but
that does NOT imply all BBJ-associated GWAS summaries are disjoint.
"""
import csv,json
from pathlib import Path
BASE=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
COHORTS=[
 {"study":"Koyanagi_2024_Japanese_alcohol","population":"Japanese","source":"Japanese alcohol GWAS meta",
  "sample_size_max":175672,"BBJ_contributing_subjects":134993,"other_components":"HERPACC,J-MICC,JPHC,TMM,Nagahama",
  "source_ref":"https://pmc.ncbi.nlm.nih.gov/articles/PMC10816704/",
  "bbj_presence":"VERIFIED_PUBLISHED","participant_level_independence":"UNKNOWN"},
 {"study":"SakaueKanai_2021_BBJ_BP","population":"Japanese","source":"BBJ PheWeb",
  "sample_size_max":145515,"BBJ_contributing_subjects":"","other_components":"BBJ",
  "source_ref":"https://pheweb.jp/variant/12-112241766-G-A",
  "bbj_presence":"VERIFIED","participant_level_independence":"UNKNOWN"},
 {"study":"BBJ_ischemic_stroke","population":"Japanese","source":"BBJ GWAS",
  "sample_size_max":"","BBJ_contributing_subjects":"","other_components":"BBJ",
  "source_ref":"https://biobankjp.org/en/",
  "bbj_presence":"VERIFIED","participant_level_independence":"UNKNOWN"},
 {"study":"GIGASTROKE_2022_EAS_AIS","population":"East Asian","source":"GIGASTROKE EAS AIS meta GWAS",
  "sample_size_max":"","BBJ_contributing_subjects":"","other_components":"Multiple EAS stroke sources",
  "source_ref":"https://pmc.ncbi.nlm.nih.gov/articles/PMC9524349/",
  "bbj_presence":"CONTRIBUTION_TO_META_UNRESOLVED_FROM_STUDY_LEVEL_SUMMARY",
  "participant_level_independence":"UNKNOWN"},
 {"study":"GIGASTROKE_2022_BBJ_PGS_evaluation","population":"Japanese",
  "source":"BBJ held-out PGS evaluation",
  "sample_size_max":41929,"BBJ_contributing_subjects":41929,"other_components":"1470 AIS cases / 40459 controls",
  "source_ref":"https://pmc.ncbi.nlm.nih.gov/articles/PMC9524349/",
  "bbj_presence":"VERIFIED_HELD_OUT_FROM_GIGASTROKE_META",
  "participant_level_independence":"BBJ_PGS_CASES_CONTROLS_EXCLUDED_FROM_GIGASTROKE_META"},
 {"study":"KCPS2_2025_Korean_multi_trait","population":"Korean","source":"KCPS2 n153950",
  "sample_size_max":153950,"BBJ_contributing_subjects":0,
  "other_components":"Korean Cancer Prevention Study-II Biobank",
  "source_ref":"https://www.nature.com/articles/s41467-025-59950-5",
  "bbj_presence":"NO_BBJ_IN_SOURCE_STUDY",
  "participant_level_independence":"NO_INDIVIDUAL_OVERLAP_WITH_JAPANESE_BBJ_ASSUMED_NOT_ATTESTED"}
]
PAIRWISE=[
 ("Koyanagi_2024_Japanese_alcohol","SakaueKanai_2021_BBJ_BP",
  "SHARED_BBJ_COHORT","NO","UNKNOWN","NO_INDEPENDENT_TWO_SAMPLE_MR"),
 ("Koyanagi_2024_Japanese_alcohol","BBJ_ischemic_stroke",
  "SHARED_BBJ_COHORT","NO","UNKNOWN","NO_INDEPENDENT_TWO_SAMPLE_MR"),
 ("Koyanagi_2024_Japanese_alcohol","GIGASTROKE_2022_EAS_AIS",
  "EXACT_META_SOURCE_OVERLAP_UNRESOLVED","NO","UNKNOWN","NO_INDEPENDENT_TWO_SAMPLE_MR"),
 ("SakaueKanai_2021_BBJ_BP","BBJ_ischemic_stroke",
  "SHARED_BBJ_COHORT","NO","UNKNOWN","NO_INDEPENDENT_TWO_SAMPLE_MR"),
 ("Koyanagi_2024_Japanese_alcohol","KCPS2_2025_Korean_multi_trait",
  "SEPARATE_POPULATION_STUDIES","YES_POPULATION","UNKNOWN","ALCOHOL_EXPOSURE_CROSS_COHORT_COMPARISON_ONLY"),
 ("KCPS2_2025_Korean_multi_trait","GIGASTROKE_2022_EAS_AIS",
  "EXACT_PARTICIPANT_OVERLAP_UNRESOLVED","NO","UNKNOWN","NO_INDEPENDENT_STROKE_REPLICATION"),
 ("GIGASTROKE_2022_BBJ_PGS_evaluation","GIGASTROKE_2022_EAS_AIS",
  "EXPLICIT_STUDY_HELD_OUT","YES_FOR_PGS_EVALUATION_ONLY","NO_BY_STUDY_METHOD",
  "HELD_OUT_PGS_EVALUATION_NOT_EQUIVALENT_TO_BBJ_GWAS_REPLICATION")
]
def build(out):
    if len(set(x["study"] for x in COHORTS))!=len(COHORTS):
        raise ValueError("Duplicate study")
    study={x["study"] for x in COHORTS}
    if any(a not in study or b not in study for a,b,*_ in PAIRWISE):
        raise ValueError("Unrecognized study")
    out.mkdir(exist_ok=True,parents=True)
    with (out/"IS_RS671_SOURCE_COHORT_PROVENANCE.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(COHORTS[0]),delimiter="\t")
        w.writeheader();w.writerows(COHORTS)
    pairdata=[dict(source_a=a,source_b=b,published_membership_risk=risk,
         source_cohorts_completely_independent_confirmed=independent,
         actual_participant_count_overlap=actual,methodology_gate=gate)
         for a,b,risk,independent,actual,gate in PAIRWISE]
    with (out/"IS_RS671_PAIRWISE_COHORT_OVERLAP_GATE.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(pairdata[0]),delimiter="\t")
        w.writeheader();w.writerows(pairdata)
    total=175672;BBJ=134993
    result={"cohort_records":len(COHORTS),"pairwise_gates":len(pairdata),
      "koyanagi_total_drinking_n":total,"koyanagi_BBJ_n_as_published":BBJ,
      "BBJ_share_of_drinking_sample_count_approx":round(BBJ/total,6),
      "individual_participant_overlap_between_studies_measured":False,
      "gigastroke_BBJ_PGS_evaluation_held_out":True,
      "held_out_PGS_implies_GIGASTROKE_EAS_summary_excludes_all_BBJ":False,
      "two_sample_MR_independence_attested":False,
      "final_status":"COHORT_MEMBERSHIP_RISK_CHARACTERIZED_NOT_ACTUAL_OVERLAP_QUANTIFIED"}
    (out/"IS_RS671_COHORT_OVERLAP_AUDIT_SUMMARY.json").write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    return result
if __name__=="__main__":build(BASE)
