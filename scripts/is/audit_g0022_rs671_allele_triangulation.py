#!/usr/bin/env python3
"""Triangulate rs671-A allele associations in real Japanese alcohol/BP + EAS AIS.

This is NOT mediation or drug target validation. Requires *two* full MD5-
verified alcohol GWAS files and original per-variant allele harmonization.
Exposures use different phenotype scales; they cannot be numerically compared
or pooled without an identified causal model and study-overlap covariance.
"""
import argparse,csv,json,math
from pathlib import Path
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
SNP="12:112241766:G:A"
SOURCES={
 "alcohol_intake":ROOT/"G0022_4CS_KOYANAGI2024_ALCOHOL_INTAKE_EXPOSURE.tsv",
 "drinking_status":ROOT/"G0022_4CS_KOYANAGI2024_DRINKING_STATUS_EXPOSURE.tsv",
 "SBP":ROOT/"G0022_RS671_BBJ_BLOOD_PRESSURE_EFFECTS.tsv",
 "AIS":ROOT/"G0022_RS671_SIX_GWAS_SUBTYPE_EFFECTS.tsv",
}
def rows(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def process(root):
    r={}
    for exposure in ("alcohol_intake","drinking_status"):
        source=root/("G0022_4CS_KOYANAGI2024_"+exposure.upper()+"_EXPOSURE.tsv")
        manifest=root/("G0022_KOYANAGI2024_"+exposure.upper()+"_AUDIT.json")
        if not source.exists() or not manifest.exists():
            raise FileNotFoundError(f"FULL_EXPOSURE_INPUT_NOT_VERIFIED:{exposure}")
        meta=json.loads(manifest.read_text())
        if not meta.get("official_MD5_verified") or not meta.get("rs671_present"):
            raise ValueError("Alcohol source MD5 or rs671 allele QC missing")
        hits=[x for x in rows(source) if x["variant_grch37"]==SNP]
        if len(hits)!=1 or hits[0]["harmonized_effect_allele_ALT"]!="A":
            raise ValueError("Alcohol rs671 duplicate or allele mismatch")
        v=hits[0]
        b=float(v["beta_ALT"]);se=float(v["se"]);p=float(v["p"])
        if not all(math.isfinite(x) for x in (b,se,p)) or se<=0 or not 0<=p<=1:
            raise ValueError("Nonfinite exposure association")
        r[exposure]={"trait":exposure,"cohort":"Japanese Koyanagi 2024",
             "effect_allele":"A","beta":b,"se":se,"p":p,
             "per_variant_n":int(v["metaGWAS_variant_n"]),
             "effect_units":"log2(g/day+1)" if exposure=="alcohol_intake" else
                        "binary_drinker_status_ORIENTATION_NOT_VERIFIED",
             "independent_of_BBJ_stroke":"NO_NOT_PROVEN",
             "causal_pathway":"NOT_IDENTIFIED"}
    bp=[x for x in rows(root/"G0022_RS671_BBJ_BLOOD_PRESSURE_EFFECTS.tsv") if x["phenotype"]=="SBP"]
    if len(bp)!=1 or bp[0]["variant_grch37"]!=SNP:raise ValueError("Missing Japanese BP rs671 A")
    v=bp[0]
    r["SBP"]={"trait":"SBP","cohort":"BBJ Sakaue Kanai 2021",
         "effect_allele":"A","beta":float(v["beta_in_reported_trait_scale"]),
         "se":float(v["se"]),"p":float(v["p"]),"per_variant_n":int(v["n_total"]),
         "effect_units":"original_transformed_BP_non_mmHg",
         "independent_of_BBJ_stroke":"NO_NOT_PROVEN","causal_pathway":"NOT_IDENTIFIED"}
    stroke=[x for x in rows(root/"G0022_RS671_SIX_GWAS_SUBTYPE_EFFECTS.tsv") if x["trait"]=="EAS_AIS"]
    if len(stroke)!=1 or stroke[0]["variant"]!=SNP:raise ValueError("Wrong EAS AIS rs671")
    v=stroke[0]
    r["AIS"]={"trait":"AIS","cohort":"GIGASTROKE East Asian 2022",
         "effect_allele":"A","beta":float(v["ALT_beta"]),"se":float(v["se"]),
         "p":float(v["p"]),"per_variant_n":"",
         "effect_units":"log_odds_AIS","independent_of_BBJ_stroke":"NO_NOT_PROVEN",
         "causal_pathway":"NOT_IDENTIFIED"}
    if any(q["effect_allele"]!="A" for q in r.values()) or set(r)!=set(SOURCES):
        raise ValueError("Effect allele / required trait contract")
    output=root/"G0022_RS671_ALCOHOL_BP_AIS_ALLELE_TRIANGULATION.tsv"
    out=list(r.values())
    with output.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter="\t")
        w.writeheader();w.writerows(out)
    summary={"traits_harmonized_ALT_A":list(r),
       "alcohol_intake_A_decreases_log2_grams_per_day":r["alcohol_intake"]["beta"]<0,
       "SBP_A_decreases_reported_BP_scale":r["SBP"]["beta"]<0,
       "EAS_AIS_A_decreases_log_odds":r["AIS"]["beta"]<0,
       "drinking_status_case_direction_verified":False,
       "rs671_A_as_alcohol_only_valid_IV":False,
       "single_SNP_pleiotropy_test_possible":False,
       "EAS_AIS_and_BBJ_exposure_cohort_independent":False,
       "causal_mediation_percent_calculated":False,
       "causal_ALDH2_stroke_gene_validated":False,
       "status":"GENETIC_ASSOCIATION_DIRECTION_CONCORDANCE_ONLY"}
    (root/"G0022_RS671_ALCOHOL_BP_AIS_TRIANGULATION_SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    for x in out:print(x["trait"],"A beta",x["beta"],"se",x["se"],"p",x["p"],x["effect_units"])
    return summary
if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",type=Path,default=ROOT)
    args=ap.parse_args()
    process(args.root)
