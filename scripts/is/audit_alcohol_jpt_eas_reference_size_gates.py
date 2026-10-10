#!/usr/bin/env python3
"""Research-science fail-closed audit of 1KG JPT104 vs EAS504 LD mismatch.

Includes exact source-aligned LD reconstruction, matched SNP set sample-N
sensitivity, fixed reproducibility seeds, and caveat that s is not a
universal pass threshold or a causal test.
"""
import csv,json,math
from pathlib import Path
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/alcohol_conditional_susie_v1")
LOCI=("ADH1B","ALDH2")
def table(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def audit(root=ROOT):
    rows=table(root/"ALCOHOL_ALDH2_ADH1B_EAS504_VS_JPT104_LD_MISMATCH.tsv")
    fixed=table(root/"ALCOHOL_EAS504_JPT104_FIXED_IDENTICAL_SNP_LD_MISMATCH.tsv")
    controls=table(root/"ALCOHOL_EAS104_RANDOM_DOWNSAMPLING_LD_MISMATCH.tsv")
    summary=json.loads((root/"ALCOHOL_JPT104_EAS104_PANEL_SIZE_CONTROL_SUMMARY.json").read_text())
    if len(rows)!=12 or len(fixed)!=4 or len(controls)!=24 or summary["repeats"]!=12:
        raise ValueError("EAS/JPT LD mismatch controls incomplete")
    if any(x["MR_or_causal_finemapping_performed"]!="FALSE" for x in controls):
        raise ValueError("Causal inference unexpectedly claimed")
    result=[]
    for locus in LOCI:
        jqc=json.loads((root/locus/"JPT104_sensitivity"/"input_qc.json").read_text())
        eqc=json.loads((root/locus/"input_qc.json").read_text())
        originally=json.loads((root/locus/"ALCOHOL_SUSIE_DIAGNOSTIC_AND_GATE.json").read_text())
        if (jqc["status"]!="JPT104_LD_SOURCE_RECONSTRUCTED_EAS_BASELINE_CROSSCHECK_PASS" or
            jqc["max_reconstructed_existing_EAS_signed_LD_absolute_difference"]>1e-9 or
            eqc["reference_EAS_n"]!=504 or jqc["JPT_n"]!=104):
            raise ValueError("Source 1000G individual genotypes or allele alignment uncertain")
        cases=[x for x in fixed if x["locus"]==locus]
        if len(cases)!=2 or {x["reference"] for x in cases}!={"EAS504","JPT104"}:
            raise ValueError("Missing true fixed same-SNP set JPT/EAS LD comparison")
        e=next(x for x in cases if x["reference"]=="EAS504")
        j=next(x for x in cases if x["reference"]=="JPT104")
        if not all(x["source_genomic_variants_identical_across_references"]=="TRUE" for x in cases):
            raise ValueError("JPT/EAS source SNPs differ, cannot compare")
        if abs(float(e["LD_mismatch_s"])-float(originally["LD_mismatch_s"]))>1e-6:
            raise ValueError("Original primary mismatch unable to reproduce")
        if any(int(x["used_snp_count"])!=eqc["variants"] for x in cases):
            raise ValueError("Different JPT/EAS SNP sets vs original reference")
        gs=[x for x in controls if x["locus"]==locus]
        if len(gs)!=12 or {int(x["replicate_index"]) for x in gs}!=set(range(1,13)):
            raise ValueError("Reference downsample control replicates incomplete")
        if {int(x["seed"]) for x in gs}!={20261010}:
            raise ValueError("Unproven reproducible 104 sample draw provenance")
        ss=[float(x["diagnostic_s"]) for x in gs]
        if any(not math.isfinite(x) or x<0 or x>1 for x in ss):
            raise ValueError("N104 LD mismatch invalid")
        lo,hi=min(ss),max(ss)
        jpt=float(j["LD_mismatch_s"])
        row={
          "locus":locus,
          "GWAS_Japanese_per_variant_N_median":originally["median_sample_n"],
          "original_SNPs_in_EAS_signed_LD":eqc["variants"],
          "original_EAS504_s":float(e["LD_mismatch_s"]),
          "JPT104_s":jpt,
          "same_SNP_set_for_EAS504_JPT104_and_EAS104":True,
          "EAS104_random_panel_size_repeats":len(ss),
          "EAS104_s_median":sorted(ss)[len(ss)//2-1:len(ss)//2+1],
          "EAS104_s_min":lo,"EAS104_s_max":hi,
          "JPT104_s_within_EAS104_random_control_range":lo<=jpt<=hi,
          "LD_reconstructed_from_existing_EAS_snp_genotypes":True,
          "sample_size_experiment_ancestry_matched_Japanese_LD_proven":False,
          "ALDH2_mismatch_causally_explained":False,
          "validated_credible_sets":0,
          "causal_alcohol_AIS_or_BP_mediation_estimated":False}
        a,b=row["EAS104_s_median"];row["EAS104_s_median"]=(a+b)/2
        result.append(row)
    if any(x["validated_credible_sets"]!=0 for x in result):
        raise ValueError("Unexpected causal fine-map promotion")
    with (root/"ALCOHOL_ALDH2_ADH1B_REF_SIZE_GATE.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(result[0]),delimiter="\t")
        w.writeheader();w.writerows(result)
    verdict={
      "status":"REAL_1KG_JPT104_VS_EAS504_DOWNSAMPLED_SIZE_SENSITIVITY_COMPLETE",
      "genotype_allele_source_AUDIT":"EAS504_SIGNED_LD_RECONSTRUCTION_EXACT_PASS",
      "comparison_basis":"EXACT_IDENTICAL_SNP_LIST_NO_ANCESTRY_SPECIFIC_AF_REMOVAL",
      "n_loci":2,"n_JPT_subjects":104,"n_EAS_subjects":504,
      "n_EAS104_control_repeats":12,
      "ADH1B_results":result[0],"ALDH2_results":result[1],
      "ALDH2_finemap_blocked":True,
      "ADH1B_only_exploratory_source_finemap":True,
      "panel_size_104_may_explain_jpt_s_inflation":True,
      "panel_size_104_proven_only_cause_of_jpt_s":False,
      "sample_replicates_are_independent_of_each_other":False,
      "genomewide_ALDH2_2225_IS_universe_changed":False,
      "cohort_matched_adequate_Japanese_LD_available":False,
      "susie_finite_reference_correction_applied":False,
      "quantitative_two_sample_MR_or_mediation_performed":False,
      "causal_alcohol_to_AIS_pleiotropy_excluded":False}
    (root/"IS_ALCOHOL_JPT_EAS_REFERENCE_SIZE_SCIENTIFIC_GATE.json").write_text(json.dumps(verdict,indent=2))
    print(json.dumps(verdict,indent=2));return verdict
if __name__=="__main__":audit()
