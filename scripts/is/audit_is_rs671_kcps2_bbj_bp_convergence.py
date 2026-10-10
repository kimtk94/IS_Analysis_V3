#!/usr/bin/env python3
"""Korean KCPS2 vs Japanese BBJ rs671-A blood-pressure direction audit.

Compatible allele/genome coordinates and direction; NOT common mmol/Hg effect
unit or disjoint outcome samples. Require actual source-verification manifests.
"""
import csv,json
from pathlib import Path
BASE=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
def records(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def run(root):
    src={x["phenotype"]:x for x in records(root/"G0022_RS671_BBJ_BLOOD_PRESSURE_EFFECTS.tsv")}
    x=[];paired=[]
    for trait in ("SBP","DBP"):
        report=root/f"G0022_RS671_KCPS2_KOREAN_{trait}_GWAS_SUMMARY.json"
        path=root/f"G0022_RS671_KCPS2_KOREAN_{trait}_GWAS.tsv"
        if not path.exists() or not report.exists():
            continue
        r=json.loads(report.read_text())
        if (not r.get("selected_member_CRC32_verified")
            or not r.get("selected_member_SHA256_verified")):
            raise ValueError("Unverified original Korean BP member")
        z=records(path)
        if len(z)!=1 or z[0]["variant"]!="12:112241766:G:A" or z[0]["harmonized_effect_allele"]!="A":
            raise ValueError("Korean rs671 allele gate")
        j=src.get(trait)
        if not j or j["variant_grch37"]!="12:112241766:G:A" or j["effect_allele"]!="ALT_A":
            raise ValueError("Japanese rs671 allele gate")
        a=float(z[0]["KCPS2_beta_A"]);b=float(j["beta_in_reported_trait_scale"])
        x.append({"BP_trait":trait,"variant_grch37":"12:112241766:G:A",
         "effect_allele":"A","korean_KCPS2_beta_A":a,"korean_se":float(z[0]["KCPS2_se"]),
         "korean_p":float(z[0]["KCPS2_p"]),"korean_n":int(z[0]["KCPS2_per_variant_N"]),
         "japanese_BBJ_beta_A":b,"japanese_se":float(j["se"]),
         "japanese_p":float(j["p"]),"japanese_n":int(j["n_total"]),
         "sign_consistent":int((a<0)==(b<0)),
         "between_cohort_beta_magnitude_comparable":False,
         "cohort_participant_overlap_attested":False,
         "BP_mediation_causal_claim":False})
        paired.append(trait)
    if not x:raise FileNotFoundError("No completed KCPS2 BP source traits")
    with (root/"G0022_RS671_KCPS2_BBJ_BP_DIRECTION_COMPARISON.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(x[0]),delimiter="\t");w.writeheader();w.writerows(x)
    result={"verified_BP_trait_pairs":paired,
      "korean_japanese_BP_beta_direction_concordant":all(x["sign_consistent"] for x in x),
      "korean_BBJ_BP_beta_units_identical":False,
      "AIS_mediated_alcohol_effect_identified":False,
      "same_stroke_study_independence_verified":False,
      "status":"CROSS_JAPAN_KOREA_BP_ALLELE_ASSOCIATION_DIRECTION_ONLY"}
    (root/"G0022_RS671_KCPS2_BBJ_BP_CONVERGENCE_SUMMARY.json").write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    return result
if __name__=="__main__":run(BASE)
