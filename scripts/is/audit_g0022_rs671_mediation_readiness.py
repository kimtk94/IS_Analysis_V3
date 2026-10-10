#!/usr/bin/env python3
"""ALDH2 rs671 causal-pathway prerequisites; blocker-aware machine-readable gates.

Do not estimate MR, mediation, sex-interaction, or QTL coloc from only one SNP.
All files read-only; new TSV/JSON are derived evidence documents only.
"""
import argparse,csv,json
from pathlib import Path
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
PATHWAYS=[
 ("catalytic_loss_of_function","rs671 functional missense ALDH2 catalytic function",
  "EXTERNAL_BIOCHEMICAL_EVIDENCE","Catalytic assay evidence is external literature","No formal stroke direct effect estimate"),
 ("alcohol_intake","rs671 changes drinking tolerance, quantity and pattern",
  "EXTERNAL_HUMAN_GENE_ENVIRONMENT_ASSOCIATION","Sex-stratified alcohol GWAS beta/SE and measurement units",
  "Cannot separately identify alcohol-mediated stroke effect"),
 ("blood_pressure","rs671 may affect BP via alcohol and other pathways",
  "BIOLOGICAL_HYPOTHESIS","EAS sex-stratified SBP/DBP genetic associations and mediator timing",
  "Cannot estimate BP mediation from rs671 stroke summary alone"),
 ("direct_vascular","ALDH2 detoxification or oxidative stress influences vessels",
  "BIOLOGICAL_HYPOTHESIS","Independent tissue/biochemical perturbation and BP/alcohol-adjusted analysis",
  "No direct vascular causal attribution"),
 ("transcription","rs671 and LD variants associate with JCTF blood ALDH2 mRNA",
  "EAS_JCTF_VARIANT_LEVEL_EQTL","Complete unfiltered cis-eQTL summary and matched-ancestry LD",
  "No multi-signal AIS-eQTL colocalization"),
 ("measured_plasma_protein","rs671 associates with JCTF Olink-measured ALDH2 protein",
  "EAS_JCTF_VARIANT_LEVEL_PQTL","Orthogonal LC-MS protein abundance or independent binding epitope test",
  "Measured pQTL may reflect missense epitope rather than real abundance"),
 ("stroke_subtype_specific","AIS/SVS/CES/LAS summary effect estimates may differ",
  "EAS_STUDY_LEVEL_ASSOCIATION","Disjoint subtype/control cohort covariance and adequately powered sample N",
  "Cannot compare overlapping case-control subgroup effect differences formally"),
 ("sex_specific_effect","ALDH2 influence may vary by sex and drinking behavior",
  "EXTERNAL_GENE_ENVIRONMENT_EVIDENCE","Sex-stratified GWAS and QTL and gene-by-drinking interaction N",
  "Sex-stratified stroke effects not available in current canonical sources")
]
def run(root):
    p=json.loads((root/"G0022_RS671_STROKE_SUBTYPE_SUMMARY.json").read_text())
    j=json.loads((root/"G0022_EAS_JCTF_INTEGRATED_EVIDENCE_SUMMARY.json").read_text())
    cond=json.loads((root/"G0022_AIS_RS671_CONDITIONAL_DIAGNOSTIC_SUMMARY.json").read_text())
    if p["rows"]!=24 or not p["all_24_study_variant_effects_protective"]:
        raise RuntimeError("Source phenotype reconciliation incomplete")
    if j["molecular_colocalization_complete"] or j["causal_ALDH2_stroke_mechanism_resolved"]:
        raise RuntimeError("Unexpected molecular causal claim")
    if cond["reference_n"]!=504:
        raise RuntimeError("Unrecognized LD provenance")
    out=[]
    for key,hypothesis,evidence,needed,limit in PATHWAYS:
        out.append({
           "pathway_key":key,"hypothesis":hypothesis,
           "existing_evidence_type":evidence,"required_to_test_causality":needed,
           "specific_unresolved_limitation":limit,
           "independent_instrument_count_verified":0,
           "matched_ALDH2_AIS_gene_coloc":False,
           "sex_interaction_result_available":False,
           "status":"HYPOTHESIS_AND_DATA_READINESS_ONLY_NOT_CAUSAL"})
    path=root/"G0022_RS671_MEDIATION_PATHWAY_READINESS.tsv"
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter="\t")
        w.writeheader();w.writerows(out)
    result={
      "source_locus":"G0022","variant":"12:112241766:G:A",
      "hypotheses_considered":len(out),
      "observed_external_evidence_Japanese_ALDH2_eQTL_variants":j["East_Asian_JCTF_QTL_variants_with_ALDH2_eQTL"],
      "observed_external_evidence_Japanese_ALDH2_pQTL_variants":j["East_Asian_JCTF_QTL_variants_with_ALDH2_pQTL"],
      "eas_AIS_OR_A_allele":p["rs671_AIS_OR"],
      "available_stroke_subtype_GWAS_sources":p["rs671_phenotypes"],
      "available_sex_stratified_stroke_GWAS":False,
      "available_harmonized_alcohol_consumption_GWAS":False,
      "available_harmonized_blood_pressure_GWAS":False,
      "ancestry_independent_instruments_for_MVMR_verified":0,
      "eligible_single_rs671_MR_Egger":False,
      "eligible_single_rs671_MR_pleiotropy_test":False,
      "eligible_coloc_susie_e_p":False,
      "eligible_multivariable_mediation":False,
      "valid_mediation_effects_computed":0,
      "scientific_gate":"BLOCKED_DO_NOT_ESTIMATE_RS671_TO_QTL_TO_STROKE_MEDIATION_WITHOUT_COMPLETE_INPUTS"}
    (root/"G0022_RS671_MEDIATION_READINESS_SUMMARY.json").write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    return result
if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--root",type=Path,default=ROOT)
    args=parser.parse_args()
    run(args.root)
