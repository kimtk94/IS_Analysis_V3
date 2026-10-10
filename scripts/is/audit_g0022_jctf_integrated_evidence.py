#!/usr/bin/env python3
"""G0022 integrated coding / Japanese molecular-QTL / ancestry coverage ledger.

Ranks variant evidence for thesis discussion WITHOUT promoting any causal gene,
mediation model, stroke MR, or valid coloc. All sources must be separately stamped.
"""
import argparse,csv,json
from pathlib import Path
BASE=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
ALT=BASE/"alternate_qtl_sources_v1"
JOB=BASE/"jctf_japan_omics_variants_v1"
GTEX=BASE/"eqtl_catalogue_multitissue_v1"
def load(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def build(base,alt,job,gtex):
    prior=json.loads((base/"G0022_FULL_AIS_PUBLICATION_GATE_SUMMARY.json").read_text())
    if prior.get("causal_gene_publication_ready") or prior.get("causal_signal_publication_ready"):
        raise ValueError("Source publication gate changed")
    ref=json.loads((alt/"G0022_NON_GTEX_QTL_CS_COVERAGE_SUMMARY.json").read_text())
    if ref.get("cs_variant_dataset_tests")!=24 or ref.get("source_test_status",{})!={"POSITION_ABSENT":24}:
        raise ValueError("Non-GTEx QTL source coverage changed; reanalysis necessary")
    gt=json.loads((gtex/"G0022_FULL_AIS_FOUR_CS_QTL_SOURCE_COVERAGE_SUMMARY.json").read_text())
    if gt.get("by_coverage_status",{})!={"NOT_TESTED_VARIANT_POSITION_ABSENT":16}:
        raise ValueError("GTEx assay status changed")
    c=json.loads((base/"G0022_AIS_RS671_CONDITIONAL_DIAGNOSTIC_SUMMARY.json").read_text())
    if c.get("conditional_gws_p_lt_5e8")!=0 or c.get("reference_n")!=504:
        raise ValueError("Single-lead LD evidence changed")
    j=json.loads((job/"G0022_JCTF_4SNP_QTL_EVIDENCE_SUMMARY.json").read_text())
    if j["ALDH2_eQTL_variant_records"]!=4 or j["ALDH2_pQTL_variant_records"]!=4:
        raise ValueError("Incomplete EAS JCTF QTL")
    if j["colocalization_numerically_calculated"] or j["validated_causal_stroke_genes"]:
        raise RuntimeError("Unjustified JCTF causal claim")
    ld={x["variant_id"]:x for x in load(base/"G0022_AIS_FOUR_CS_RS671_LD.tsv")}
    data=load(job/"G0022_JCTF_4SNP_JAPANESE_EQTL_PQTL_EVIDENCE.tsv")
    e={}
    for r in data:
        key=(r["GWAS_variant_GRCh37"],r["gene_symbol"],r["qtl_category"])
        if key in e:raise ValueError("Ambiguous Japanese QTL gene-category")
        e[key]=r
    if len(ld)!=4:raise ValueError("GWAS CS dimension not four")
    output=[]
    for v in sorted(ld):
        eqtl=e.get((v,"ALDH2","eQTL"))
        pqtl=e.get((v,"ALDH2","pQTL"))
        if eqtl is None or pqtl is None:raise RuntimeError("Missing ALDH2 QTL for GWAS CS variant")
        output.append({
          "gwas_variant_grch37":v,"nearest_positional_gene_only":ld[v]["annotation"],
          "gwas_AIS_CSPIP":eqtl["gwas_ais_full_locus_pip"],
          "rs671_1000gEAS_r2":ld[v]["r2_to_rs671"],
          "Japanese_JCTF_ALDH2_eQTL_p":eqtl["jctf_allele_association_p"],
          "Japanese_JCTF_ALDH2_eQTL_SuSiE_PIP":eqtl["jctf_susie_variant_pip"],
          "Japanese_JCTF_ALDH2_pQTL_p":pqtl["jctf_allele_association_p"],
          "Japanese_JCTF_ALDH2_pQTL_SuSiE_PIP":pqtl["jctf_susie_variant_pip"],
          "N_non_JCTF_QTL_source_datasets_without_position":10,
          "Japanese_variant_molecular_association":"PRESENT_EAS_JCTF",
          "causal_mechanism":"NOT_DISENTANGLED_FROM_RS671_CODING_ALCOHOL_BP_PATHWAYS",
          "assay_caveat":"ALDH2_MISSENSE_CAN_AFFECT_PROTEIN_ASSAY_BINDING",
          "independent_stroke_replication":"NOT_ESTABLISHED",
          "full_locus_QTL_coloc":"NOT_TESTED_ONLY_VARIANT_LEVEL_ASSOCIATION",
          "gene_causality":"NOT_ESTABLISHED"})
    with (base/"G0022_EAS_JCTF_INTEGRATED_VARIANT_EVIDENCE.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(output[0]),delimiter="\t")
        w.writeheader();w.writerows(output)
    report={
       "genetic_variant_credible_set":4,
       "variant_rs671_coding":"12:112241766:G:A",
       "conditional_external_LD_GWS_after_rs671":0,
       "other_public_QTL_datasets_queried":10,
       "non_JCTF_QTL_4_SNP_source_checks_without_position":40,
       "East_Asian_JCTF_QTL_variants_with_ALDH2_eQTL":4,
       "East_Asian_JCTF_QTL_variants_with_ALDH2_pQTL":4,
       "Japanese_JCTF_eQTL_sample_size_study_reported":1019,
       "Japanese_JCTF_pQTL_sample_size_study_reported":1384,
       "JCTF_release_filtered_associations_only":True,
       "JCTF_variant_PIPs_do_not_prove_shared_stroke_causal_signal":True,
       "molecular_colocalization_complete":False,
       "causal_ALDH2_stroke_mechanism_resolved":False,
       "research_priority":"ALDH2_RS671_CODING_VS_PROTEIN_ABUNDANCE_VS_ALCOHOL_BP_MEDIATION",
       "result_status":"INTEGRATED_CANDIDATE_EVIDENCE_NOT_CAUSAL_VALIDATION"}
    (base/"G0022_EAS_JCTF_INTEGRATED_EVIDENCE_SUMMARY.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    return report
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--base",type=Path,default=BASE)
    p.add_argument("--alt",type=Path,default=ALT)
    p.add_argument("--job",type=Path,default=JOB)
    p.add_argument("--gtex",type=Path,default=GTEX)
    a=p.parse_args()
    build(a.base,a.alt,a.job,a.gtex)
