#!/usr/bin/env python3
"""Japan/Korea rs671 cross-cohort *association direction* replication ledger.

Two ancestries, incomparable transformed exposure beta units, no joint z,
no fixed-effect meta analysis, no mediation. Korean BP is included only when
actual verified source TSV + audit exists.
"""
import argparse,csv,json,math
from pathlib import Path
BASE=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
SNP="12:112241766:G:A"
def csvrows(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def build(base):
    ja=csvrows(base/"G0022_4CS_KOYANAGI2024_ALCOHOL_INTAKE_EXPOSURE.tsv")
    ko=csvrows(base/"G0022_RS671_KCPS2_KOREAN_ALCOHOL_AMOUNT.tsv")
    j=[x for x in ja if x["variant_grch37"]==SNP and x["harmonized_effect_allele_ALT"]=="A"]
    k=[x for x in ko if x["effect_allele_harmonized"]=="A"
       and x["variant"]=="ALDH2_rs671_GRCh37_12_112241766_G_A"]
    if len(j)!=1 or len(k)!=1:raise ValueError("Original Japanese/Korean rs671 allele anchors unavailable")
    ko_meta=json.loads((base/"G0022_RS671_KCPS2_KOREAN_SOURCE_AUDIT.json").read_text())
    if not ko_meta["selected_ZIP_member_CRC32_verified"] or not ko_meta["selected_member_SHA256_verified"]:
        raise ValueError("Korean association source ZIP member not verified")
    samples=[{
       "study":"Koyanagi_2024_Japanese","ancestry":"Japanese","phenotype":"alcohol_intake",
       "effect_allele_A_beta":float(j[0]["beta_ALT"]),
       "se":float(j[0]["se"]),"p":float(j[0]["p"]),
       "frequency_A":float(j[0]["ALT_EAF"]),
       "variant_n":int(j[0]["metaGWAS_variant_n"]),
       "source_trait_scale":"log2(grams_per_day+1)",
       "source_credibility":"FULL_ZENODO_SOURCE_MD5_VERIFIED"},
     {
       "study":"Jee_2025_KCPS2_Korean","ancestry":"Korean","phenotype":"alcohol_intake",
       "effect_allele_A_beta":float(k[0]["KCPS2_alcohol_beta_A"]),
       "se":float(k[0]["KCPS2_alcohol_se"]),
       "p":float(k[0]["KCPS2_ALCO_AMOUNT_p"]),
       "frequency_A":float(k[0]["KCPS2_af_A"]),
       "variant_n":int(k[0]["KCPS2_sample_size_at_variant"]),
       "source_trait_scale":"inverse_rank_normal_transformed",
       "source_credibility":"ORIGINAL_ZIP_SELECTED_MEMBER_CRC32_SHA256_VERIFIED"}]
    included=[]
    for phenotype in ("SBP","DBP"):
        a=base/f"G0022_RS671_KCPS2_KOREAN_{phenotype}_GWAS_SUMMARY.json"
        b=base/f"G0022_RS671_KCPS2_KOREAN_{phenotype}_GWAS.tsv"
        if a.exists()!=b.exists():raise RuntimeError("KCPS2 BP summary/source partial state")
        if not a.exists():continue
        m=json.loads(a.read_text());rows=csvrows(b)
        if len(rows)!=1 or not m["selected_member_CRC32_verified"] or not m["selected_member_SHA256_verified"]:
            raise ValueError("KCPS2 BP source audit missing")
        z=rows[0]
        if z["harmonized_effect_allele"]!="A" or z["variant"]!=SNP:
            raise ValueError("KCPS2 BP allele mismatch")
        samples.append({"study":"Jee_2025_KCPS2_Korean","ancestry":"Korean","phenotype":phenotype,
          "effect_allele_A_beta":float(z["KCPS2_beta_A"]),
          "se":float(z["KCPS2_se"]),"p":float(z["KCPS2_p"]),
          "frequency_A":float(z["KCPS2_af_A"]),
          "variant_n":int(z["KCPS2_per_variant_N"]),
          "source_trait_scale":"inverse_rank_normal_transformed_not_mmHg",
          "source_credibility":"ORIGINAL_ZIP_SELECTED_MEMBER_CRC32_SHA256_VERIFIED"})
        included.append(phenotype)
    if any(not math.isfinite(x["effect_allele_A_beta"]) for x in samples):
        raise ValueError("Invalid source effect")
    p=base/"G0022_RS671_KCPS2_JAPAN_CROSSCOHORT_COMPARISON.tsv"
    with p.open("w",newline="") as f:
        w=csv.DictWriter(f,delimiter="\t",fieldnames=list(samples[0]))
        w.writeheader();w.writerows(samples)
    verdict={
       "Japanese_alcohol_beta_ALT_A":samples[0]["effect_allele_A_beta"],
       "Korean_alcohol_beta_ALT_A":samples[1]["effect_allele_A_beta"],
       "alcohol_association_direction_both_negative":samples[0]["effect_allele_A_beta"]<0 and samples[1]["effect_allele_A_beta"]<0,
       "japanese_A_frequency":samples[0]["frequency_A"],
       "korean_A_frequency":samples[1]["frequency_A"],
       "A_frequency_difference_Japanese_minus_Korean":round(samples[0]["frequency_A"]-samples[1]["frequency_A"],6),
       "Japanese_Korean_beta_magnitude_test_valid":False,
       "trait_transformations_harmonized_to_common_units":False,
       "Korean_BP_traits_verified":included,
       "KCPS2_SNP_independent_of_Japanese_subjects_by_cohort_identity":True,
       "independent_EAS_AIS_stroke_replication":False,
       "causal_alcohol_to_BP_to_stroke_mediation_tested":False,
       "scope":"ALCOHOL_EXPOSURE_DIRECTION_REPLICATION_ACROSS_EAS_POPULATIONS_NO_CAUSAL_MEDIATION"}
    (base/"G0022_RS671_KCPS2_JAPAN_CROSSCOHORT_SUMMARY.json").write_text(json.dumps(verdict,indent=2))
    print(json.dumps(verdict,indent=2))
    return verdict
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,default=BASE)
    build(p.parse_args().root)
