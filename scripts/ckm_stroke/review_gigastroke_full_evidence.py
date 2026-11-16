#!/usr/bin/env python3
"""Independent fail-closed aggregate review of complete GIGASTROKE x BMI source QC.

Consumes only the bounded JSON report emitted by audit_gigastroke_bmi_allele_harmonization.
Does NOT validate source bytes itself; the producing job separately checks full MD5.
No individual SNP associations are emitted; no LDSC h2/rg/MR is performed.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

EXPECTED_BMI_CANDIDATE=1_007_899
EXPECTED_BMI_SHA="ed8f713827cd3281bea5f88a9cb42a3eabfecb1ee67c912ca7e86f151e4399bf"
PINNED={
    "AIS":("GCST90104540","1bfe4ae8a5042fb24cb0a562b05b0d2f",283_482_604),
    "LAA":("GCST90104542","71ea9fd7c06791a831861053b6a1cfbe",230_156_887),
    "CE":("GCST90104541","c252262586a6ea84b985f30dd7ab32be",263_057_516),
    "SVD":("GCST90104543","f41bc5f55b3cab238e305e3e86363f69",237_118_943),
}
DIRECTIONS=("same","swapped","complement_same","complement_swapped")
OBSERVED_COUNTS={
    "AIS":(7_482_032,1_006_932,967),
    "LAA":(5_784_788,956_929,50_970),
    "CE":(6_659_793,996_662,11_237),
    "SVD":(5_974_028,967_356,40_543),
}
CASES={"AIS":62_100,"LAA":6_399,"CE":10_804,"SVD":6_811}
CONTROLS=1_234_808


def _int(counts,key):
    value=counts.get(key,0)
    if isinstance(value,bool) or not isinstance(value,int) or value<0:
        raise ValueError("nonnegative integer expected for "+key)
    return value


def review(report:dict)->dict:
    if report.get("schema")!="ckm-stroke-gigastroke-eur-full-md5-bmi-harmonization-v1":
        raise ValueError("unexpected source report schema")
    if report.get("technical_status")!="PASS_ALL_REQUESTED_ORIGINAL_MD5_AND_ALLELE_COUNTS":
        raise ValueError("source technical status is not PASS")
    if report.get("ready_for_ldsc") is not False or report.get("scientific_status")!="HOLD":
        raise ValueError("incorrect scientific gate state")
    if report.get("LDSC_h2_or_rg_executed") is not False:
        raise ValueError("LDSC calculations should not have run")
    if report.get("source_archives_retained") is not False:
        raise ValueError("source archives retained contrary to policy")
    bmi=report.get("BMI_reference_panel",{})
    if (bmi.get("candidate_total")!=EXPECTED_BMI_CANDIDATE
        or bmi.get("candidate_plaintext_sha256")!=EXPECTED_BMI_SHA):
        raise ValueError("BMI source candidate changed")
    outcomes=report.get("outcomes")
    if not isinstance(outcomes,list) or len(outcomes)!=4:
        raise ValueError("expected four outcomes")
    if {x.get("phenotype") for x in outcomes}!=set(PINNED):
        raise ValueError("duplicate or missing outcome phenotype")
    summary=[]
    for a in sorted(outcomes,key=lambda x:["AIS","LAA","CE","SVD"].index(x["phenotype"])):
        ph=a["phenotype"]
        acc,expected_md5,size=PINNED[ph]
        if (a.get("accession")!=acc or a.get("original_compressed_bytes")!=size
            or a.get("published_md5")!=expected_md5
            or a.get("raw_md5_verified") is not True
            or a.get("publisher_yaml_md5_sidecar_verified") is not True):
            raise ValueError("publisher original pin, bytes, or MD5 not verified for "+ph)
        if a.get("technical_status")!="PASS_RAW_MD5_AND_ALLELE_CLASSIFICATION_ONLY":
            raise ValueError("not a successful full origin read: "+ph)
        digest=a.get("source_sha256","")
        if len(digest)!=64 or any(ch not in "0123456789abcdef" for ch in digest):
            raise ValueError("missing compressed SHA256 evidence: "+ph)
        if (a.get("outcome_has_rsid_column") is not False
            or a.get("outcome_per_variant_N_available") is not False):
            raise ValueError("unexpected source schema metadata: "+ph)
        c=a.get("counts")
        if not isinstance(c,dict):
            raise ValueError("missing counts: "+ph)
        total=_int(c,"original_outcome_rows")
        cand_hits=_int(c,"outcome_rows_with_BMI_coordinate")
        out_not_bmi=_int(c,"outcome_rows_absent_from_BMI_candidate")
        concordant=_int(c,"unique_harmonizable_BMI_positions")
        unique_coord=_int(c,"unique_outcome_coordinate_matches")
        not_outcome=_int(c,"candidate_coordinates_not_in_outcome")
        collisions=_int(c,"outcome_collision_coordinates_excluded")
        if total!=cand_hits+out_not_bmi:
            raise ValueError("outcome row partition not balanced: "+ph)
        if (unique_coord+not_outcome+collisions!=EXPECTED_BMI_CANDIDATE
            or cand_hits<unique_coord or concordant>unique_coord):
            raise ValueError("BMI coordinate accounting mismatch: "+ph)
        if concordant!=sum(_int(c,"unique_"+d) for d in DIRECTIONS):
            raise ValueError("directional category sum mismatch: "+ph)
        if (_int(c,"unique_orientation_sign_positive")!=
            _int(c,"unique_same")+_int(c,"unique_complement_same")):
            raise ValueError("positive sign count mismatch: "+ph)
        if (_int(c,"unique_orientation_sign_negative")!=
            _int(c,"unique_swapped")+_int(c,"unique_complement_swapped")):
            raise ValueError("negative sign count mismatch: "+ph)
        expected_total,expected_matched,expected_missing=OBSERVED_COUNTS[ph]
        if (total,concordant,not_outcome)!=(expected_total,expected_matched,expected_missing):
            raise ValueError("expected pinned original source counts changed: "+ph)
        cases=CASES[ph]
        # This is a population-level arithmetic illustration, NOT per-SNP N and NOT an approved LDSC input.
        illustrative_neff=4*cases*CONTROLS/(cases+CONTROLS)
        all_one_orientation=(concordant==_int(c,"unique_same"))
        summary.append({
            "phenotype":ph,"accession":acc,
            "publisher_original_md5":expected_md5,
            "original_sha256":digest,
            "original_compressed_bytes":size,
            "original_rows":total,
            "matched_and_allele_aligned":concordant,
            "candidate_not_in_outcome":not_outcome,
            "match_rate_pct":round(100*concordant/EXPECTED_BMI_CANDIDATE,4),
            "unique_same":_int(c,"unique_same"),
            "unique_swapped":_int(c,"unique_swapped"),
            "unique_complement_same":_int(c,"unique_complement_same"),
            "unique_complement_swapped":_int(c,"unique_complement_swapped"),
            "collision_coordinates_excluded":collisions,
            "all_matched_in_same_orientation":all_one_orientation,
            "reported_cases":cases,"reported_controls":CONTROLS,
            "illustrative_total_N":cases+CONTROLS,
            "illustrative_Neff_4_case_control_harmonic":round(illustrative_neff,2),
        })
    return {
        "schema":"ckm-stroke-gigastroke-four-verified-aggregate-review-v1",
        "technical_status":"PASS_PUBLISHER_MD5_AND_AGGREGATE_INVARIANTS",
        "BMI_candidate_rows":EXPECTED_BMI_CANDIDATE,
        "BMI_candidate_plaintext_sha256":EXPECTED_BMI_SHA,
        "outcomes":summary,
        "uniform_allele_orientation_needs_independent_review":
            all(x["all_matched_in_same_orientation"] for x in summary),
        "illustrative_Neff_is_NOT_per_variant_or_approved_for_LDSC":True,
        "scientific_status":"HOLD",
        "ready_for_ldsc":False,
        "LDSC_h2_rg_MR_performed":False,
        "remaining_gates":[
            "Independent ancestry-stratified effect-allele/orientation assessment",
            "Per-outcome effective sample size justification and case/control heterogeneity",
            "GIANT/UKB vs GIGASTROKE participant overlap risk",
            "Per-SNP outcome missingness and effect-size/SE scale QC",
            "LDSC h2/intercept diagnostics before rg interpretation"
        ],
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--evidence",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    args=ap.parse_args()
    try:
        evidence=json.loads(args.evidence.read_text(encoding="utf-8"))
        output=review(evidence)
    except Exception as exc:
        output={
            "technical_status":"HOLD_FAILED_EVIDENCE_REVIEW",
            "error_type":type(exc).__name__,
            "error_detail":str(exc)[:220],
            "scientific_status":"HOLD",
            "ready_for_ldsc":False,
        }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(output,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":output["technical_status"],
                      "outcomes":[{"phenotype":x["phenotype"],
                                   "aligned":x["matched_and_allele_aligned"]}
                                  for x in output.get("outcomes",[])],
                      "science":output["scientific_status"]}))
    if output["technical_status"]!="PASS_PUBLISHER_MD5_AND_AGGREGATE_INVARIANTS":
        raise SystemExit(1)


if __name__=="__main__":
    main()
