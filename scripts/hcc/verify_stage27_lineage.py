#!/usr/bin/env python3
"""Stage27 conservative source-code and numeric rank transformation review."""
import argparse,json,math,pathlib

COMMIT="438f52271d0a86a6701baa7a96e80c516ba575c2"
FULL="461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a"
FULLZ="9c1e4e05472a0aa477b81aa7b4c4e1b27809df9a2b71c1f57e247bb031263007"
def code_review(v):
    if v.get("schema")!="hcc-stage27-source-code-provenance-v1" or v.get("source_commit")!=COMMIT:
        raise ValueError("pinned HCC source mismatch")
    scripts=v.get("audited_named_scripts",{})
    if set(scripts)!={"HCC_prior.R","HCCclutsering.R","HCCfullclustering.R",
                   "helpresfns.R","Tutorial_bnClustOmics.R","README.txt"}:
        raise ValueError("source manifest missing")
    if not v.get("HCCfullz_clustering_script_load_verified"):
        raise ValueError("source loader not found")
    for f,d in scripts.items():
        if len(d.get("sha256",""))!=64 or type(d.get("size_bytes"))!=int or d["size_bytes"]<100:
            raise ValueError("source script not SHA-pinned")
    if v.get("HCCfullz_source_file_writer_mentions_total")!=0:
        raise ValueError("a potential RDS writer must be reviewed manually")
    for key in ("HCCfullz_generation_code_verified_in_source_repository",
       "HCCfullz_generation_code_absent_elsewhere_proven",
       "original_Ng_2022_Fig2b_gene_universe_recreated",
       "source_patient_gene_or_expression_exported"):
        if v.get(key) is not False:raise ValueError("overstated method or privacy claim")
    return {"schema":"hcc-stage27-code-postflight-v1",
      "immutable_public_suter_code_verified":True,
      "clustering_source_HCCfullz_consumer_verified":True,
      "preprocessing_generation_code_in_repository_verified":False,
      "preprocessing_code_cannot_exist_elsewhere_claimed":False,
      "V03":"PARTIAL"}
def ranks_review(v):
    if v.get("schema")!="hcc-stage27-full-fullz-paired-rank-v1" or v.get("source_commit")!=COMMIT or v.get("source_full_sha256")!=FULL or v.get("source_fullz_sha256")!=FULLZ:
        raise ValueError("public source hash or schema mismatch")
    metrics=v.get("layer_metrics")
    if not isinstance(metrics,list) or len(metrics)!=5 or {x.get("modality") for x in metrics}!={"M","CN","T","P","PP"}:
        raise ValueError("missing five modalities")
    for x in metrics:
        if x.get("biopsies_n")!=50 or type(x.get("source_common_features_n"))!=int or x["source_common_features_n"]<20:
            raise ValueError("wrong source sample/feature count")
        n=x.get("rho_finite_features_n");m=x.get("undefined_rho_feature_n")
        if type(n)!=int or type(m)!=int or n<20 or n+m!=x["source_common_features_n"]:
            raise ValueError("invalid complete rank feature count")
        for key in ("spearman_rank_correlation_median","spearman_rank_correlation_p10",
                    "spearman_rank_correlation_p90","pearson_median"):
            val=x.get(key)
            if type(val) not in (int,float) or not math.isfinite(val) or not -1<=val<=1:
                raise ValueError("invalid rho")
        if not x["spearman_rank_correlation_p10"]<=x["spearman_rank_correlation_median"]<=x["spearman_rank_correlation_p90"]:
            raise ValueError("correlation quantiles unsorted")
        for key in ("spearman_fraction_ge_0_99","spearman_fraction_lt_0_9",
                    "spearman_fraction_negative","pearson_fraction_ge_0_99"):
            val=x.get(key)
            if type(val) not in (int,float) or not math.isfinite(val) or not 0<=val<=1:
                raise ValueError("invalid fraction")
        if x.get("feature_or_sample_identifiers_exported") is not False:
            raise ValueError("source feature/patient disclosure")
    for key in ("exact_source_matrix_generation_reproduced",
                "original_Figure2b_gene_universe_reproduced",
                "original_Ng_2022_tumor_normal_relative_quantitation_reproduced",
                "gene_level_FDR_or_biological_significance_tested",
                "per_feature_rhos_or_patient_identifiers_exported"):
        if v.get(key) is not False:raise ValueError("forbidden science/privacy claim")
    return {"schema":"hcc-stage27-rank-postflight-v1",
      "HCCfull_and_fullz_gene_rank_relations_audited":True,
      "T_P_PP_rank_medians":{
         x["modality"]:x["spearman_rank_correlation_median"]
         for x in metrics if x["modality"] in ("T","P","PP")},
      "source_processing_code_reproduction":"HOLD",
      "original_Fig2b_gene_universe":"HOLD",
      "V03":"PARTIAL","V07":"PARTIAL"}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--kind",required=True,choices=("source-code","ranks"))
    p.add_argument("--input",required=True);p.add_argument("--out",required=True)
    a=p.parse_args()
    x=json.loads(pathlib.Path(a.input).read_text(encoding="utf8"))
    y=code_review(x) if a.kind=="source-code" else ranks_review(x)
    f=pathlib.Path(a.out);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(y,indent=2)+"\n",encoding="utf8")
    print(json.dumps(y,indent=2))
if __name__=="__main__":main()
