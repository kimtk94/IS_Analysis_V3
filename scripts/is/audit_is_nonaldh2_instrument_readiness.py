#!/usr/bin/env python3
"""Strength-vs-validity ledger for positional non-ALDH2 alcohol GWAS candidates."""
import csv,json,math
from pathlib import Path
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
GENE={
 "4:100239319:T:C":("ADH1B","direct_ethanol_oxidation"),
 "2:27730940:T:C":("GCKR","metabolism_pleiotropy"),
 "9:38395928:T:C":("ALDH1B1","aldehyde_metabolism_pleiotropy"),
 "9:75461066:T:C":("ALDH1A1_nearby","aldehyde_retinoid_metabolism"),
 "4:39413780:A:G":("KLB","FGF21_endocrine_metabolism"),
 "12:106750302:A:G":("chr12_unresolved","unknown_molecular_pleiotropy")
}
def audit(root):
    src=root/"IS_RS671_NON_ALDH2_ALCOHOL_POSITIONAL_GWS_AIS_QC.tsv"
    with src.open() as f:rows=list(csv.DictReader(f,delimiter="\t"))
    if set(r["variant_grch37"] for r in rows)!=set(GENE):
        raise ValueError("Source candidate identities changed; manual annotation review needed")
    out=[]
    for r in rows:
        v=r["variant_grch37"]
        f=float(r["alcohol_beta_over_se_abs"])**2
        if f<10:raise ValueError("Weak marginal exposure effect")
        match=r["ais_qc"]=="ALLELE_HARMONIZED"
        eaf_diff=(abs(float(r["alcohol_ALT_eaf"])-float(r["ais_ALT_eaf"]))
                  if match else None)
        if eaf_diff is not None and eaf_diff>.1:
            raise ValueError("Excessive East Asian/Japanese allele freq mismatch")
        sym,pleiotropy=GENE[v]
        out.append({"variant_grch37":v,"associated_locus_annotation":sym,
         "approx_marginal_F_beta2_over_se2":round(f,3),
         "has_source_ais_allele_match":int(match),
         "alcohol_beta_ALT":r["alcohol_beta_ALT"],
         "AIS_beta_ALT":r["ais_beta_ALT"],
         "AIS_original_p":r["ais_p"],
         "EAF_delta_japanese_vs_EAS":eaf_diff if eaf_diff is not None else "",
         "source_pleiotropy_risk_category":pleiotropy,
         "strength_F_gt_10_does_not_certify_valid_IV":True,
         "actual_EAS_LD_clumped":False,
         "residual_horizontal_pleiotropy_excluded":False,
         "individual_sample_overlap_attested":False,
         "mr_eligible_instrument":"NO"})
    with (root/"IS_RS671_NONALDH2_INSTRUMENT_STRENGTH_AND_VALIDITY.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter="\t")
        w.writeheader();w.writerows(out)
    matches=[x for x in out if x["has_source_ais_allele_match"]]
    report={
      "positional_sentinels":len(out),"exact_allele_AIS_matches":len(matches),
      "minimum_marginal_F":min(x["approx_marginal_F_beta2_over_se2"] for x in out),
      "max_japan_EAS_ALT_frequency_difference":max(x["EAF_delta_japanese_vs_EAS"] for x in matches),
      "AIS_nominal_p_less_0_05":sum(float(x["AIS_original_p"])<.05 for x in matches),
      "valid_independent_instruments_verified":0,
      "LD_clumped_independent_instruments_verified":0,
      "causal_MR_outcome_calculated":False,
      "pleiotropic_metabolism_loci_present":True,
      "status":"EXPOSURE_SIGNIFICANCE_ONLY_ALL_GENETIC_IV_VALIDITY_GATES_OPEN"}
    (root/"IS_RS671_NONALDH2_INSTRUMENT_READINESS_SUMMARY.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    return report
if __name__=="__main__":audit(ROOT)
