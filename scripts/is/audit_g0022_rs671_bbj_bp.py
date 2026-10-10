#!/usr/bin/env python3
"""rs671 PheWeb public Japanese BP GWAS: provenance, effect allele, study overlap.

Reads full official BBJ PheWeb variant-page embedded JSON; PheWeb phenotype
beta is ALT A (page ref/alt). No mmHg-scale assumption and no mediation MR.
"""
import argparse,csv,hashlib,json,math,re
from pathlib import Path
SOURCE=Path("/srv/is-analysis/data/is/gwas_exposure/japanese_bp_sakaue2021/rs671_BBJ_PheWeb.html")
OUT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
PHENOS=("SBP","DBP","MAP","PP")
def source_variant(html):
    m=re.search(r"window\.variant\s*=\s*(\{[^\n]+\});",html)
    if not m:raise ValueError("Missing official PheWeb variant payload")
    obj=json.loads(m.group(1))
    if (str(obj.get("chrom")),int(obj.get("pos")),"".join((obj.get("ref"),obj.get("alt"))))!=("12",112241766,"GA"):
        raise ValueError("Wrong variant or REF/ALT on BP source")
    if not isinstance(obj["phenos"],list):raise ValueError("Missing phenotype array")
    return obj
def extract(source,out):
    obj=source_variant(source.read_text())
    found={x["phenocode"]:x for x in obj["phenos"] if x.get("phenocode") in PHENOS}
    if set(found)!=set(PHENOS):raise ValueError("Unexpected missing BBJ BP panel")
    records=[]
    for ph in PHENOS:
        r=found[ph]
        if r["citation"]!="SakaueKanai2021" or r["category"]!="Blood pressure":
            raise ValueError("Wrong trait study provenance")
        b=float(r["beta"]);se=float(r["sebeta"]);p=float(r["pval"]);af=float(r["af"]);n=int(r["num_samples"])
        if not all(math.isfinite(x) for x in (b,se,p,af)) or se<=0 or p<=0 or p>1 or n<100000 or not(0<af<1):
            raise ValueError("Invalid GWAS summary")
        records.append({"phenotype":ph,"full_name":r["phenostring"],"GWAS":r["citation"],
          "study_ancestry":"Japanese BBJ","source_id":r["phenocode"],"variant_grch37":"12:112241766:G:A",
          "effect_allele":"ALT_A","beta_in_reported_trait_scale":b,"se":se,"p":p,"n_total":n,
          "bbj_alt_allele_frequency":af,"original_outcome_units":"SOURCE_TRANSFORMED_TRAIT_VERIFY_BEFORE_MR",
          "raw_mmHg_effect_known":False,"no_direct_causal_BP_to_stroke_effect_claim":True,
          "study_overlap_with_BBJ_stroke":"POSSIBLE_NOT_DISPROVEN"})
    out.mkdir(parents=True,exist_ok=True)
    with (out/"G0022_RS671_BBJ_BLOOD_PRESSURE_EFFECTS.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter="\t")
        w.writeheader();w.writerows(records)
    result={"source_url":"https://pheweb.jp/variant/12-112241766-G-A",
      "source_content_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
      "variant":"12:112241766:G:A","effect_allele":"A",
      "BP_traits":len(records),"min_BP_N":min(x["n_total"] for x in records),
      "max_BP_N":max(x["n_total"] for x in records),
      "all_BP_beta_ALT_A_negative":all(x["beta_in_reported_trait_scale"]<0 for x in records),
      "n_source_GxE_combination_phenotypes":sum(str(x.get("phenocode","")).endswith("_GxE") for x in obj["phenos"]),
      "GxE_combo_interaction_exposure":"UNSPECIFIED_DO_NOT_ASSUME_ALCOHOL",
      "stroke_cohort_independence":"UNVERIFIED",
      "mmHg_trait_units_assumed":False,"mediation_identified":False,
      "evidence_status":"GENETIC_BP_ASSOCIATION_NOT_CAUSAL_MEDIATION"}
    (out/"G0022_RS671_BBJ_BP_PHEWEB_SOURCE_AUDIT.json").write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    for r in records:print(r["phenotype"],r["beta_in_reported_trait_scale"],r["se"],r["p"])
    return result
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,default=SOURCE)
    p.add_argument("--out",type=Path,default=OUT)
    a=p.parse_args();extract(a.source,a.out)
