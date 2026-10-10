#!/usr/bin/env python3
"""Strength and EAF proxy-QC for positional non-ALDH2 Japanese alcohol hits.

Single-SNP F_proxy=(beta/se)^2 is association strength, NOT a conditional
multi-instrument F statistic. Physical distance is NOT population LD. EAF
between Japanese exposure and mixed EAS AIS can differ by ancestry;
EAF concordance is not independence proof. No MR effect estimated.
"""
import argparse,csv,json,math
from pathlib import Path
BASE=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
INPUT=BASE/"IS_RS671_NON_ALDH2_ALCOHOL_POSITIONAL_GWS_AIS_QC.tsv"
def audit(source,out):
    with source.open() as f:rows=list(csv.DictReader(f,delimiter="\t"))
    if len(rows)!=6:raise ValueError("Unexpected candidate positional lead count")
    items=[]
    for r in rows:
        b=float(r["alcohol_beta_ALT"]);se=float(r["alcohol_se"])
        if not math.isfinite(b) or not math.isfinite(se) or se<=0:
            raise ValueError("Non-finite alcohol effects")
        strength=(b/se)**2
        if r["independent_genetic_IV_certified"]!="False":
            raise ValueError("Source falsely declared LD independence")
        quality=r["ais_qc"]
        if quality=="ALLELE_HARMONIZED":
            p=float(r["ais_p"]);af_eas=float(r["ais_ALT_eaf"])
            af_jpn=float(r["alcohol_ALT_eaf"])
            delta=abs(af_eas-af_jpn)
            if delta>.10:status="EAF_DIFF_GT_10PP_NEEDS_ANCESTRY_REVIEW"
            else:status="ALLELE_MATCH_EAF_WITHIN_10PP"
        elif quality=="NOT_IN_EAS_AIS_SOURCE":
            status="OUTCOME_SNP_MISSING_NO_BETA_IMPUTED";p=None;delta=None
        else:raise ValueError("Unanticipated EAS AIS allele status "+quality)
        items.append({
           "variant_GRCh37":r["variant_grch37"],
           "Japanese_exposure_beta_ALT":b,"Japanese_exposure_se":se,
           "Japanese_exposure_F_proxy_beta_over_se_squared":strength,
           "Japanese_alcohol_p":r["alcohol_p"],
           "Japanese_alcohol_ALT_EAF":r["alcohol_ALT_eaf"],
           "EAS_AIS_ALT_EAF":r.get("ais_ALT_eaf",""),
           "ALT_EAF_absolute_difference":"" if delta is None else delta,
           "AIS_p":"" if p is None else p,
           "allele_reference_status":status,
           "JPN_vs_EAS_participant_overlap_ruled_out":False,
           "LD_independent_IV_certified":False,
           "causal_exclusion_restriction_certified":False,
           "MR_eligibility":"NOT_YET_ELIGIBLE_POSITIONAL_SNP_AND_PLEIOTROPY_GATE"})
    out.mkdir(exist_ok=True,parents=True)
    with (out/"IS_RS671_NON_ALDH2_ALCOHOL_IV_READINESS.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,delimiter="\t",fieldnames=list(items[0]))
        w.writeheader();w.writerows(items)
    found=[r for r in items if r["allele_reference_status"]=="ALLELE_MATCH_EAF_WITHIN_10PP"]
    result={"candidate_positional_loci":len(items),
      "AIS_exact_allele_EAF_within_10pp":len(found),
      "AIS_SNP_missing":sum(r["allele_reference_status"]=="OUTCOME_SNP_MISSING_NO_BETA_IMPUTED" for r in items),
      "max_abs_ALT_EAF_difference_matching":max(float(r["ALT_EAF_absolute_difference"]) for r in found),
      "F_proxy_min_across_all":min(x["Japanese_exposure_F_proxy_beta_over_se_squared"] for x in items),
      "F_proxy_max_across_all":max(x["Japanese_exposure_F_proxy_beta_over_se_squared"] for x in items),
      "LD_independence_validated":False,
      "pleiotropy_tested":False,
      "sample_overlap_ruled_out":False,
      "formal_multi_variant_conditional_F_available":False,
      "MR_computed":False,
      "status":"EXPOSURE_STRENGTH_PROXIES_ONLY_NOT_VALIDATED_INDEPENDENT_MR_INSTRUMENTS"}
    (out/"IS_RS671_NON_ALDH2_ALCOHOL_IV_READINESS_SUMMARY.json").write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    return result
if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--source",type=Path,default=INPUT)
    parser.add_argument("--out",type=Path,default=BASE)
    opts=parser.parse_args();audit(opts.source,opts.out)
