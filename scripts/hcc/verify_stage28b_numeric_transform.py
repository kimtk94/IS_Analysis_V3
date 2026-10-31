#!/usr/bin/env python3
"""Stage28b source-z-score reconstruction technical-only postflight."""
import argparse,json,math,pathlib
FULL="461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a"
FULLZ="9c1e4e05472a0aa477b81aa7b4c4e1b27809df9a2b71c1f57e247bb031263007"
COMMIT="438f52271d0a86a6701baa7a96e80c516ba575c2"
EXPECTED={"T":19768,"P":5585,"PP":7889}
def verify(v):
    if v.get("schema")!="hcc-stage28b-fullz-numeric-reconstruction-v1":
        raise ValueError("invalid numeric reconstruction schema")
    for k,x in (("source_commit",COMMIT),("full_sha256",FULL),("fullz_sha256",FULLZ)):
        if v.get(k)!=x:raise ValueError("source pin drift")
    if v.get("synthetic_numeric_reconstruction_passed") is not True:
        raise ValueError("missing independent numeric synthetic")
    for k in ("source_original_Ng_tumor_normal_2022_values_reproduced",
              "metadata_verified_canonical_sample_ids",
              "source_files_overwritten",
              "public_artifact_contains_sample_gene_or_permutation_vector"):
        if v.get(k) is not False:
            raise ValueError("incorrect source identity or patient exposure claim: "+k)
    rows=v.get("layers")
    if not isinstance(rows,list) or len(rows)!=3 or {x.get("modality") for x in rows}!=set(EXPECTED):
        raise ValueError("missing T/P/PP continuous modalities")
    for r in rows:
        if r.get("common_features_n")!=EXPECTED[r["modality"]] or r.get("sample_n")!=50:
            raise ValueError("unexpected source matrix feature/sample count")
        if r.get("training_features_n")+r.get("test_features_n")!=EXPECTED[r["modality"]]:
            raise ValueError("invalid disjoint feature holdout")
        if r.get("mapping_bijective") is not True:
            raise ValueError("sample permutation not one-to-one")
        k=r.get("unchanged_named_positions_n")
        if type(k) is not int or not 0<=k<=50:
            raise ValueError("invalid source column position matches")
        for key in ("heldout_median_abs_z_residual_before",
                    "heldout_median_abs_z_residual_after",
                    "heldout_p95_abs_z_residual_after",
                    "heldout_max_abs_z_residual_after"):
            value=r.get(key)
            if type(value) not in (int,float) or not math.isfinite(value) or value<0:
                raise ValueError("invalid quantitative source residual: "+key)
        if not (r["heldout_median_abs_z_residual_after"]<=
                r["heldout_p95_abs_z_residual_after"]<=
                r["heldout_max_abs_z_residual_after"]):
            raise ValueError("unsorted absolute residual distribution")
        for key in ("heldout_cell_fraction_abs_lt_1e6_after",
                    "heldout_cell_fraction_abs_lt_1e4_after",
                    "heldout_feature_fraction_row_max_abs_lt_1e5_after"):
            value=r.get(key)
            if type(value) not in (int,float) or not math.isfinite(value) or not 0<=value<=1:
                raise ValueError("invalid source cell threshold fraction")
        if r.get("sample_feature_identifiers_exported") is not False:
            raise ValueError("source IDs leakage")
    return {
        "schema":"hcc-stage28b-source-numeric-transform-postflight-v1",
        "source_SHA_and_50_sample_axes_verified":True,
        "independent_feature_holdout_used":True,
        "per_modality_heldout_median_abs_z_residual_after":{
            r["modality"]:r["heldout_median_abs_z_residual_after"] for r in rows},
        "per_modality_heldout_cell_fraction_lt_1e_minus6_after":{
            r["modality"]:r["heldout_cell_fraction_abs_lt_1e6_after"] for r in rows},
        "2022_original_Ng_assay_preprocessing_reproduction":"HOLD",
        "clinical_patient_biopsy_correct_label_orientation":"HOLD",
        "V03":"PARTIAL","V07":"PARTIAL",
        "biological_significance_claimed":False
    }
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",required=True);ap.add_argument("--out",required=True)
    a=ap.parse_args()
    x=json.loads(pathlib.Path(a.input).read_text(encoding="utf-8"))
    y=verify(x)
    p=pathlib.Path(a.out);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(y,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(y,indent=2))
if __name__=="__main__":main()
