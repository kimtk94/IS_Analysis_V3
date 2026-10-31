#!/usr/bin/env python3
"""Conservative postflight for source-pinned Stage28 crossmodal sample mapping."""
import argparse
import json
import math
import pathlib
from itertools import combinations

FULL="461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a"
FULLZ="9c1e4e05472a0aa477b81aa7b4c4e1b27809df9a2b71c1f57e247bb031263007"
COMMIT="438f52271d0a86a6701baa7a96e80c516ba575c2"
MODS={"M","CN","T","P","PP"}

def postflight(v):
    if v.get("schema")!="hcc-stage28-crossmodal-trainheldout-permutation-v1":
        raise ValueError("source schema mismatch")
    for key,sha in (("full_sha256",FULL),("fullz_sha256",FULLZ),("source_commit",COMMIT)):
        if v.get(key)!=sha:raise ValueError("source provenance drift: "+key)
    for key in ("synthetic_permutation_and_heldout_recovery_passed",
                "each_modality_reported_independently",
                "technical_alignment_is_not_a_biological_replication"):
        if v.get(key) is not True:raise ValueError("source QA key is false: "+key)
    for key in ("publisher_original_biopsy_identity_reconciled",
                "original_suter_preprocessing_generation_reproduced",
                "sample_labels_persistently_corrected",
                "per_sample_mapping_and_per_gene_values_exported"):
        if v.get(key) is not False:raise ValueError("overstated provenance: "+key)
    rows=v.get("modalities")
    if not isinstance(rows,list) or len(rows)!=5 or {x.get("modality") for x in rows}!=MODS:
        raise ValueError("five modality diagnostics required")
    wanted={
        "M":3519,"CN":20339,"T":19768,"P":5585,"PP":7889
    }
    for row in rows:
        m=row["modality"]
        if row.get("features_common_n")!=wanted[m] or row.get("sample_n")!=50:
            raise ValueError("unexpected source modality shape")
        n=row["features_common_n"]
        if row.get("train_feature_n")+row.get("heldout_feature_n")!=n:
            raise ValueError("source split incomplete")
        for f in ("train_vs_reverse_map_agreement_n","unchanged_named_positions_n"):
            k=row.get(f)
            if type(k) is not int or not 0<=k<=50:
                raise ValueError("invalid sample match statistic: "+f)
        for f in ("train_map_is_one_to_one","reverse_heldout_map_is_one_to_one"):
            if row.get(f) not in (True,False):raise ValueError("missing one-to-one status")
        for f in ("reordered_data_persisted","sample_permutation_or_gene_identifiers_exported"):
            if row.get(f) is not False:raise ValueError("source leakage risk")
        for prefix in ("test","reverse_test"):
            count=row.get(prefix+"_finite_feature_n")
            if type(count) is not int or not 30<=count<=n:
                raise ValueError("heldout correlation count invalid")
            for name in ("rho_before_median","rho_after_median","rho_after_q10",
                         "rho_after_q90","median_within_feature_rho_delta"):
                fld=prefix+"_"+name
                value=row.get(fld)
                if type(value) not in (int,float) or not math.isfinite(value) or not -2<=value<=2:
                    raise ValueError("invalid Spearman metric: "+fld)
            if not row[prefix+"_rho_after_q10"]<=row[prefix+"_rho_after_median"]<=row[prefix+"_rho_after_q90"]:
                raise ValueError("unsorted rho quantile")
            for name in ("frac_rho_gt_095_before","frac_rho_gt_095_after","frac_rho_negative_after"):
                fld=prefix+"_"+name
                value=row.get(fld)
                if type(value) not in (int,float) or not math.isfinite(value) or not 0<=value<=1:
                    raise ValueError("invalid positive source proportion: "+fld)
    pairs=v.get("crossmodal_comparisons")
    expected={frozenset(a) for a in (("T","P"),("T","PP"),("P","PP"),("M","CN"))}
    if not isinstance(pairs,list) or len(pairs)!=4 or {frozenset((x.get("modality_a"),x.get("modality_b"))) for x in pairs}!=expected:
        raise ValueError("invalid pairwise map comparisons")
    for p in pairs:
        k=p.get("matching_sample_map_positions_n")
        if type(k) is not int or not 0<=k<=50:
            raise ValueError("invalid paired sample mapping count")
        if p.get("both_maps_bijective") not in (True,False):
            raise ValueError("missing pair permutation status")
    bymod={r["modality"]:r for r in rows}
    bypair={frozenset((x["modality_a"],x["modality_b"])):x for x in pairs}
    return {
        "schema":"hcc-stage28-crossmodal-permutation-reviewed-v1",
        "source_RDS_hashes_verified":True,
        "independent_feature_holdout_validation":True,
        "T_P_PP_train_test_sample_mapping_agreement":{
            m:bymod[m]["train_vs_reverse_map_agreement_n"] for m in ("T","P","PP")},
        "T_P_PP_pairwise_identical_mapping_counts":{
            "_".join(p):bypair[frozenset(p)]["matching_sample_map_positions_n"]
            for p in (("T","P"),("T","PP"),("P","PP"))},
        "T_P_PP_heldout_rank_before_after_medians":{
            m:{
                "before":bymod[m]["test_rho_before_median"],
                "after":bymod[m]["test_rho_after_median"]}
            for m in ("T","P","PP")},
        "candidate_sample_permutation_validated_only_for_numeric_source":True,
        "original_biopsy_identity_metadata_verification":"HOLD",
        "preprocessing_original_source_verification":"HOLD",
        "V03":"PARTIAL","V07":"PARTIAL",
        "clinical_significance_claim_allowed":False,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input",required=True)
    p.add_argument("--out",required=True)
    a=p.parse_args()
    v=json.loads(pathlib.Path(a.input).read_text(encoding="utf-8"))
    reviewed=postflight(v)
    out=pathlib.Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(reviewed,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(reviewed,indent=2))
if __name__=="__main__":main()
