#!/usr/bin/env python3
"""Safety/reproducibility check for small HCCdata RDS label evidence."""
import argparse,json,math,pathlib
def verify(v):
  if v.get("schema")!="hcc-stage29-small-HCCdata-rds-anchor-v1" or v.get("source_commit")!="438f52271d0a86a6701baa7a96e80c516ba575c2":
    raise ValueError("source commit mismatch")
  digests={"source_HCCdata_git_blob_sha":"45538b8452a1d2e3844b596f6bd50f1086f0bf21",
   "source_full_sha256":"461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a",
   "source_fullz_sha256":"9c1e4e05472a0aa477b81aa7b4c4e1b27809df9a2b71c1f57e247bb031263007"}
  for k,d in digests.items():
    if v.get(k)!=d:raise ValueError("source fingerprint drift "+k)
  for k in ("HCCdata_is_independent_patient_cohort",
            "HCCdata_authoritative_biopsy_identity_proven",
            "source_sample_labels_corrected","clinical_or_gene_effect_inference_performed"):
    if v.get(k) is not False:raise ValueError("forbidden source inference "+k)
  r=v.get("modalities")
  if not isinstance(r,list) or len(r)!=3 or {x.get("modality") for x in r}!={"T","P","PP"}:
    raise ValueError("three source layers required")
  summaries={}
  for x in r:
    if (type(x.get("small_data_sample_n")) is not int or
        not 1<=x["small_data_sample_n"]<=500 or
        type(x.get("small_data_feature_n")) is not int or
        not 1<=x["small_data_feature_n"]<=1000):
      raise ValueError("small network source structure invalid")
    for key in ("small_data_source_numeric_matrix",
                "small_data_exact_50_biopsy_comparison_eligible"):
      if type(x.get(key)) is not bool:
        raise ValueError("source network matrix evidence missing: "+key)
    if (x["small_data_exact_50_biopsy_comparison_eligible"] and
        (x["small_data_sample_n"]!=50 or
         x["small_full_common_sample_labels_n"]!=50 or
         x["small_fullz_common_sample_labels_n"]!=50 or
         not x["small_data_source_numeric_matrix"])):
      raise ValueError("invalid comparable small network source")
    if x.get("biopsy_or_gene_identifiers_exported") is not False:raise ValueError("source label leak")
    for f in ("small_full_common_feature_n","small_fullz_common_feature_n",
              "small_full_common_sample_labels_n","small_fullz_common_sample_labels_n",
              "small_full_rho_finite_feature_n","small_fullz_rho_finite_feature_n"):
      n=x.get(f)
      upper=50 if "sample_labels_n" in f else x["small_data_feature_n"]
      if type(n)!=int or not 0<=n<=upper:raise ValueError("source matched count impossible")
    for pre in ("small_full","small_fullz"):
      for fld in ("same_label_per_feature_rho_median","fraction_gene_rho_gt_095"):
        key=pre+"_"+fld; z=x.get(key)
        if z is None:continue
        if type(z) not in (int,float) or not math.isfinite(z) or not -1<=z<=1:
          raise ValueError("invalid source rho "+key)
      if x[pre+"_rho_finite_feature_n"]==0 and x[pre+"_same_label_per_feature_rho_median"] is not None:
        raise ValueError("zero finite genes yet median exists")
    summaries[x["modality"]]={
      "same_gene_compared_n":x["small_full_common_feature_n"],
      "source_numeric_matrix":x["small_data_source_numeric_matrix"],
      "source_50_biopsy_comparison_eligible":x["small_data_exact_50_biopsy_comparison_eligible"],
      "HCCdata_vs_full_same_label_rho_median":x["small_full_same_label_per_feature_rho_median"],
      "HCCdata_vs_fullz_same_label_rho_median":x["small_fullz_same_label_per_feature_rho_median"]}
  return {"schema":"hcc-stage29-HCCdata-author-source-technical-postflight-v1",
    "small_source_reference_verified":True,"source_comparisons":summaries,
    "original_Ng_biopsy_identity_orientation":"HOLD",
    "independent_biological_cohort":False,
    "V03":"PARTIAL","V07":"PARTIAL"}
def main():
  p=argparse.ArgumentParser();p.add_argument("--input",required=True);p.add_argument("--out",required=True)
  a=p.parse_args();v=verify(json.loads(pathlib.Path(a.input).read_text()))
  tgt=pathlib.Path(a.out);tgt.parent.mkdir(parents=True,exist_ok=True)
  tgt.write_text(json.dumps(v,indent=2)+"\n")
  print(json.dumps(v,indent=2))
if __name__=="__main__":main()
