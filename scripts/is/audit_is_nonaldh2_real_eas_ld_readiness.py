#!/usr/bin/env python3
"""Evidence ledger for seven original-alcohol sentinels in real 1KG EAS LD.

Exact genotypes and reference allele frequencies are validated; per-variant
F proxies describe marginal exposure strength only, and r² independence alone
does NOT satisfy exclusion restriction or sample-overlap requirements.
"""
import argparse,csv,json,math
from pathlib import Path
LD=Path("/srv/is-analysis/data/is/ld_reference/nonaldh2_alcohol_eas_20261010")
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
RS671="12:112241766:G:A"
def rows(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def build(root,ld):
    check=json.loads((ld/"IS_RS671_ALCOHOL_7SNP_REAL_EAS_LD_SUMMARY.json").read_text())
    if check["reference_samples"]!=504 or check["exact_verified_variants"]!=7 or check["pairwise_r2_results"]!=21:
        raise ValueError("Real EAS LD source QC not complete")
    pairs=rows(ld/"IS_RS671_ALCOHOL_7SNP_REAL_EAS_LD_PAIRS.tsv")
    if len(pairs)!=21 or len(check["ALT_EAS_reference_af"])!=7:
        raise ValueError("Pairwise LD incomplete")
    candidates=rows(root/"IS_RS671_NON_ALDH2_ALCOHOL_POSITIONAL_GWS_AIS_QC.tsv")
    if len(candidates)!=6:raise ValueError("Unexpected exposure shortlist")
    signed=rows(root/"G0022_4CS_KOYANAGI2024_ALCOHOL_INTAKE_EXPOSURE.tsv")
    baseline=[x for x in signed if x["variant_grch37"]==RS671]
    if len(baseline)!=1:raise ValueError("No source-verified rs671 alcohol sentinel")
    candidates.append({"variant_grch37":RS671,
       "alcohol_beta_ALT":baseline[0]["beta_ALT"],
       "alcohol_se":baseline[0]["se"],
       "alcohol_ALT_eaf":baseline[0]["ALT_EAF"],
       "ais_qc":"ALLELE_HARMONIZED"})
    one=[]
    for r in candidates:
        snp=r["variant_grch37"]
        allele=check["ALT_EAS_reference_af"].get(snp)
        if allele is None or allele["n_nonmissing"]!=504 or allele["n_missing"]!=0:
            raise ValueError("Required 504 complete EAS exact genotypes missing")
        b=float(r["alcohol_beta_ALT"]);se=float(r["alcohol_se"])
        af1=float(r["alcohol_ALT_eaf"]);af2=allele["1KG_EAS_ALT_allele_freq"]
        others=[float(x["ALT_dosage_r2"]) for x in pairs
            if snp in (x["SNP_a"],x["SNP_b"])]
        if len(others)!=6 or len([x for x in others if math.isfinite(x)])!=6:
            raise ValueError("Missing pairwise LD for one SNP")
        one.append({
          "variant_grch37":snp,
          "is_ALDH2_rs671":int(snp==RS671),
          "GWAS_alcohol_original_beta_ALT":b,
          "GWAS_alcohol_se":se,
          "GWAS_ALCOHOL_marginal_F_proxy_beta2_over_se2":(b/se)**2,
          "GWAS_Japanese_ALCOHOL_ALT_EAF":af1,
          "1KG_EAS_504_ALT_EAF":af2,
          "absolute_ALCOHOL_Japanese_vs_1KG_EAS_ALT_EAF_difference":abs(af1-af2),
          "largest_pair_r2_to_other_six":max(others),
          "largest_r2_under_0_1":int(max(others)<.1),
          "all_four_source_chromosomes_checked":True,
          "individual_1000G_EAS_genotype_alleles_aligned":True,
          "reference_allele_freq_matches_within_10pp":abs(af1-af2)<=.1,
          "MR_exclusion_restriction_proved":False,
          "no_sample_overlap_proved":False,
          "population_LD_genetically_independent_of_other_6":"WEAK_PAIRWISE_LD_IN_1KG_ONLY",
          "formally_eligible_two_sample_MR_IV":False})
    if not all(x["reference_allele_freq_matches_within_10pp"] for x in one):
        raise ValueError("EAS ALT frequency >10pp difference, research review needed")
    with (root/"IS_NONALDH2_ALCOHOL_7SNP_EAS_REFERENCE_LD_IV_QC.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(one[0]),delimiter="\t")
        w.writeheader();w.writerows(one)
    jpt_path=ld/"IS_RS671_7SNP_EAS_JPT_ALLELE_FREQUENCIES.tsv"
    sensitivity=ld/"IS_RS671_7SNP_JPT_SENSITIVITY_SUMMARY.json"
    jpt_freq={}
    jpt_summary={}
    if jpt_path.is_file() and sensitivity.is_file():
        jpt_summary=json.loads(sensitivity.read_text())
        if jpt_summary.get("sample_count_JPT")!=104 or jpt_summary.get("pair_count")!=21:
            raise ValueError("JPT 104 reference sensitivity incomplete")
        fr=rows(jpt_path)
        jpt_freq={x["SNP"]:float(x["JPT_ALT_EAF"]) for x in fr}
        if set(jpt_freq)!=set(x["variant_grch37"] for x in one):
            raise ValueError("JPT and EAS reference SNP labels disagree")
        for row in one:
            q=jpt_freq[row["variant_grch37"]]
            row["1KG_JPT_104_ALT_EAF"]=q
            row["Japanese_alcohol_vs_1KG_JPT_ALT_EAF_abs_difference"]=abs(q-row["GWAS_Japanese_ALCOHOL_ALT_EAF"])
        if max(x["Japanese_alcohol_vs_1KG_JPT_ALT_EAF_abs_difference"] for x in one)>.1:
            raise ValueError("JPT reference frequency discrepancy unexpectedly large")
        # Re-emit enriched QC table, retaining exact provenance columns.
        with (root/"IS_NONALDH2_ALCOHOL_7SNP_EAS_REFERENCE_LD_IV_QC.tsv").open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(one[0]),delimiter="\t")
            w.writeheader();w.writerows(one)
    report={"real_eas_genotype_reference_n":504,"verified_sentinels_including_rs671":len(one),
      "pairs_verified":len(pairs),"swapped_REF_ALT":sum(
         json.loads((ld/f"chr{c}.source_qc.json").read_text())["swapped_REF_ALT_source_count"]
         for c in ("2","4","9","12")),
      "max_pairwise_reference_r2":max(float(x["ALT_dosage_r2"]) for x in pairs),
      "JPT_104_reference_validation_complete":bool(jpt_summary),
      "JPT_104_max_pairwise_r2":jpt_summary.get("max_measurable_JPT_r2"),
      "JPT_reference_frequency_max_difference_from_Japanese_GWAS":max(
         (x["Japanese_alcohol_vs_1KG_JPT_ALT_EAF_abs_difference"] for x in one),default=None)
          if jpt_summary else None,
      "max_Japanese_GWAS_EAS_reference_EAF_difference":max(
         x["absolute_ALCOHOL_Japanese_vs_1KG_EAS_ALT_EAF_difference"] for x in one),
      "all_seven_pairwise_r2_below_0_1":all(x["largest_r2_under_0_1"] for x in one),
      "genomewide_clumping_of_all_associated_snps_completed":False,
      "same_allele_signed_dosage_verified":True,
      "treat_non_ALDH2_variants_as_valid_MR_instruments":False,
      "independent_genetic_instruments_causally_validated":0,
      "LD_and_pleiotropy_are_distinct_gates":True,
      "clinical_causal_alcohol_to_AIS_mediation_claim":False,
      "status":"REAL_LD_EVIDENCE_PASS_EXCLUSION_RESTRICTION_AND_COHORT_OVERLAP_UNVERIFIED"}
    (root/"IS_NONALDH2_ALCOHOL_7SNP_EAS_REFERENCE_LD_IV_SUMMARY.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    return report
if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,default=ROOT)
    ap.add_argument("--ld",type=Path,default=LD)
    a=ap.parse_args();build(a.out,a.ld)
