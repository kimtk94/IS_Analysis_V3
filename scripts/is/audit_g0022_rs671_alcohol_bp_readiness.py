#!/usr/bin/env python3
"""Integrate verified Japanese BP and open alcohol-GWAS source readiness for G0022.

Fail closed until alcohol GWAS source MD5 passes and ALT allele effects validated.
The BBJ BP study is partially overlapping with BBJ ischemic stroke;
causal mediation, MR, or genetic interactions not computed here.
"""
import argparse,csv,hashlib,json
from pathlib import Path
BASE=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
EXPOSURE=Path("/srv/is-analysis/data/is/gwas_exposure/japanese_alcohol_2024")
def md5(path):
    h=hashlib.md5()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):h.update(block)
    return h.hexdigest()
def build(base,exposure):
    bp=json.loads((base/"G0022_RS671_BBJ_BP_PHEWEB_SOURCE_AUDIT.json").read_text())
    gwas=json.loads((base/"G0022_RS671_STROKE_SUBTYPE_SUMMARY.json").read_text())
    qtl=json.loads((base/"G0022_EAS_JCTF_INTEGRATED_EVIDENCE_SUMMARY.json").read_text())
    metadata=json.loads((exposure/"ZENODO_10038152_METADATA.json").read_text())
    source={x["key"]:x for x in metadata["files"]}
    if (bp.get("BP_traits")!=4 or bp.get("variant")!="12:112241766:G:A"
        or bp.get("mmHg_trait_units_assumed") is not False):
        raise ValueError("BP source or scale contract invalid")
    if gwas["rows"]!=24 or qtl["molecular_colocalization_complete"]:
        raise ValueError("GWAS/Japanese QTL evidence provenance changed")
    files=[
       ("daily_alcohol_intake_log2_gplus1","1_Alcohol_intake_Unstratified.tsv.gz",154570),
       ("drinking_status_never_ever","5_Drinking_Unstratified.tsv.gz",175672)]
    items=[]
    for trait,name,n in files:
        m=source[name]
        path=exposure/name;part=exposure/(name+".part")
        if not m["checksum"].startswith("md5:"):
            raise ValueError("Official Zenodo MD5 absent")
        valid=(path.is_file() and path.stat().st_size==m["size"]
             and md5(path)==m["checksum"][4:])
        assay_id=("ALCOHOL_INTAKE" if trait.startswith("daily") else "DRINKING_STATUS")
        extracted_summary=base/f"G0022_KOYANAGI2024_{assay_id}_AUDIT.json"
        extracted_file=base/f"G0022_4CS_KOYANAGI2024_{assay_id}_EXPOSURE.tsv"
        direction_verified=False
        if valid and extracted_summary.is_file() and extracted_file.is_file():
            confirmation=json.loads(extracted_summary.read_text())
            if (confirmation.get("official_MD5_verified") is not True
                or confirmation.get("source")!=name):
                raise ValueError("Source-extraction provenance disagreement")
            if confirmation.get("rs671_present"):
                with extracted_file.open() as h:
                    rr=list(csv.DictReader(h,delimiter="\t"))
                verified=[r for r in rr if r["variant_grch37"]=="12:112241766:G:A"
                     and r["harmonized_effect_allele_ALT"]=="A"]
                if len(verified)!=1:raise ValueError("rs671 allele audit disagreement")
                direction_verified=True
        items.append({
          "exposure":trait,"source_dataset":"zenodo_10038152_Koyanagi_2024_v1",
          "japanese_metaGWAS_samples":n,"genome_build":"GRCh37_hg19_VERIFY_IN_HEADER",
          "phenotype_units":("log2(grams_per_day+1)" if trait.startswith("daily") else "case_control_binary"),
          "source_file":name,"source_expected_bytes":m["size"],"source_MD5":m["checksum"],
          "downloaded_full_bytes":path.stat().st_size if path.is_file() else 0,
          "unverified_partial_bytes":part.stat().st_size if part.exists() else 0,
          "official_source_MD5_pass":int(bool(valid)),
          "rs671_ALT_effect_allele_verified":int(direction_verified),
          "BBJ_cohort_included":"YES_AS_PER_ORIGINAL_README",
          "source_status":("SOURCE_MD5_AND_RS671_ALT_HARMONIZED" if direction_verified else
            "SOURCE_MD5_VALID_BUT_SNP_UNCHECKED" if valid else "BLOCKED_FULL_SOURCE_NOT_VERIFIED")})
    out=base/"G0022_RS671_ALCOHOL_BP_EXPOSURE_READINESS.tsv"
    with out.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(items[0]),delimiter="\t")
        w.writeheader();w.writerows(items)
    state={
      "source_rs671_GWAS_phenotypes":gwas["rs671_phenotypes"],
      "verified_bbj_BP_phenotypes":bp["BP_traits"],
      "verified_BP_rs671_A_direction":"LOWER_REPORTED_TRAIT_SCALE",
      "BP_units_converted_to_mmHg":False,
      "alcohol_metaGWAS_expected_source_count":2,
      "alcohol_metaGWAS_full_source_md5_verified":sum(x["official_source_MD5_pass"] for x in items),
      "alcohol_metaGWAS_rs671_direction_allele_validated":sum(x["rs671_ALT_effect_allele_verified"] for x in items),
      "GIGASTROKE_BBJ_study_independent":False,
      "rs671_multi_trait_pleiotropy_not_ruled_out":True,
      "sex_stratified_alcohol_stroke_BP_effects_in_sources":False,
      "additional_clumped_independent_instruments_validated":0,
      "eligibility_for_ALDH2_to_alcohol_to_BP_to_stroke_MVMR":False,
      "eligible_one_SNP_MR_Egger":False,
      "eligible_coloc_susie_with_complete_cis_QTL":False,
      "causal_mediation_effect_estimated":0,
      "next":"After both source MD5 and rs671 ALT effects, descriptively triangulate alcohol consumption and BP with AIS. Causal mediation stays blocked pending independent instruments, overlap and pleiotropy checks.",
      "gate":"EXPOSURE_ASSOCIATIONS_IDENTIFIED_CAUSAL_MEDIATION_BLOCKED"}
    (base/"G0022_RS671_ALCOHOL_BP_READINESS_SUMMARY.json").write_text(json.dumps(state,indent=2))
    print(json.dumps(state,indent=2))
    return state
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--base",type=Path,default=BASE)
    p.add_argument("--exposure",type=Path,default=EXPOSURE)
    a=p.parse_args();build(a.base,a.exposure)
