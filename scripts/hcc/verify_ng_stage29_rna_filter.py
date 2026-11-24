#!/usr/bin/env python3
"""Stage29 fail-closed structural RNA feature filter analysis."""
import argparse,json,math,pathlib
COMMIT="438f52271d0a86a6701baa7a96e80c516ba575c2"
FULL="461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a"
FULLZ="9c1e4e05472a0aa477b81aa7b4c4e1b27809df9a2b71c1f57e247bb031263007"
def verify(v):
    if v.get("schema")!="hcc-stage29-rna-fullz-exclusion-qc-v1":
        raise ValueError("schema mismatch")
    for k,val in (("source_commit",COMMIT),("full_sha256",FULL),("fullz_sha256",FULLZ)):
        if v.get(k)!=val:raise ValueError("source pin failure: "+k)
    if (v.get("rna_full_feature_n"),v.get("rna_fullz_feature_n"),
        v.get("strict_matching_feature_n"),v.get("rna_full_only_feature_n"),
        v.get("feature_sample_n"))!=(20356,19768,19768,588,50):
        raise ValueError("source RNA axes drift")
    for g,n in (("excluded_rna_feature_sd",588),("retained_rna_feature_sd",19768)):
        d=v.get(g)
        if not isinstance(d,dict) or d.get("n")!=n:raise ValueError("SD group count mismatch")
        for f in ("zero_sd_n","sd_lt_1e_minus12_n","sd_lt_1e_minus6_n"):
            q=d.get(f)
            if type(q)!=int or not 0<=q<=n:raise ValueError("invalid SD threshold count")
        for f in ("sd_min","sd_p10","sd_median","sd_p90","sd_max"):
            q=d.get(f)
            if type(q) not in (float,int) or not math.isfinite(q) or q<0:
                raise ValueError("invalid numeric SD")
        if not d["sd_min"]<=d["sd_p10"]<=d["sd_median"]<=d["sd_p90"]<=d["sd_max"]:
            raise ValueError("SD quantiles out of order")
    for f in ("rna_full_nonzero_sd_feature_n","rna_full_gt_1e_minus12_sd_feature_n"):
        q=v.get(f)
        if type(q)!=int or not 0<=q<=20356:raise ValueError("nonzero feature count invalid")
    if v.get("all_588_excluded_features_exactly_constant") not in (True,False):
        raise ValueError("missing zero-SD filter decision")
    if v.get("retained_constant_feature_n")!=v["retained_rna_feature_sd"]["zero_sd_n"]:
        raise ValueError("retained source constant mismatch")
    if v["all_588_excluded_features_exactly_constant"] is True and v["excluded_rna_feature_sd"]["zero_sd_n"]!=588:
        raise ValueError("source zero-SD contradiction")
    for k in ("feature_selection_author_code_verified","per_gene_or_patient_values_exported","source_rds_files_overwritten"):
        if v.get(k) is not False:raise ValueError("unsupported source lineage claim")
    return {
       "schema":"hcc-stage29-rna-feature-filter-postflight-v1",
       "source_588_rna_exclusion_evaluated":True,
       "excluded_588_constant_feature_n":v["excluded_rna_feature_sd"]["zero_sd_n"],
       "kept_19768_constant_feature_n":v["retained_constant_feature_n"],
       "exact_zero_SD_filter_equivalence_supported":
         v["rna_full_nonzero_sd_feature_n"]==19768 and
         v["excluded_rna_feature_sd"]["zero_sd_n"]==588 and
         v["retained_constant_feature_n"]==0,
       "author_feature_selection_code_reproduced":False,
       "V03":"PARTIAL","biological_pass":False}
def main():
    p=argparse.ArgumentParser();p.add_argument("--input",required=True);p.add_argument("--out",required=True)
    a=p.parse_args()
    obj=verify(json.loads(pathlib.Path(a.input).read_text(encoding="utf-8")))
    tgt=pathlib.Path(a.out);tgt.parent.mkdir(parents=True,exist_ok=True)
    tgt.write_text(json.dumps(obj,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(obj,indent=2))
if __name__=="__main__":main()
