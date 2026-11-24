#!/usr/bin/env python3
"""Stage44: source PDC/CDAP numeric transform verification with disjoint holdouts.

The 330-column rank-based mapping is learned with EVEN source feature rows,
challenged with ODD features, and numerical affine transform fitted with EVEN
mapped columns and independently tested on ODD mapped columns. Original source
row identifiers, column names, fitted per-feature values and mapping never
leave CI runner memory. No official aliquot identity or NODE RNA inference.
"""
import argparse,datetime as dt,json,pathlib
import numpy as np
from scipy.optimize import linear_sum_assignment
from audit_stage42_cdap_pdc_numeric_prefix import pdc_prefix,report_prefix
from audit_stage43_numeric_rank_holdout import exact_unique_feature_rows,rowwise_rank_z,similarity
from audit_stage40_cdap_crossassay_source_labels import original_tmt_reports,retry

def quantile(a,p):
    return round(float(np.quantile(a,p)),8)

def evaluate(p,c):
    if p.shape!=c.shape or p.ndim!=2 or p.shape[1]!=330 or p.shape[0]<24:
        raise ValueError("unverified original source sample/feature axis")
    if not np.isfinite(p).all() or not np.isfinite(c).all():
        raise ValueError("nonfinite source feature")
    nfeatures,ncols=p.shape
    rp=rowwise_rank_z(p);rc=rowwise_rank_z(c)
    train=np.arange(0,nfeatures,2)
    test=np.arange(1,nfeatures,2)
    a=similarity(rp[train],rc[train])
    b=similarity(rp[test],rc[test])
    ia,ja=linear_sum_assignment(-a)
    ib,jb=linear_sum_assignment(-b)
    map_a=np.empty(ncols,dtype=int);map_b=np.empty(ncols,dtype=int)
    map_a[ia]=ja;map_b[ib]=jb
    if len(set(map_a.tolist()))!=ncols or len(set(map_b.tolist()))!=ncols:
        raise ValueError("not independent source permutation")
    agreement=int(np.sum(map_a==map_b))
    paired=c[:,map_a]
    # Per-feature monotone/rank signal over ALL independently mapped columns,
    # distinct from sample profile cosine over the training features.
    per_feature_spearman=np.mean(rp*rowwise_rank_z(paired),axis=1)
    shift=np.mean(rp*rowwise_rank_z(np.roll(paired,1,axis=1)),axis=1)
    # Fit a scalar slope and offset per gene in HALF of matched source columns.
    # Evaluate numerical residuals and R2 in the other 165 columns, NO feature
    # selection / sample-position reordering from heldout numerical residuals.
    sample_train=np.arange(0,ncols,2)
    sample_test=np.arange(1,ncols,2)
    xtr=p[:,sample_train];ytr=paired[:,sample_train]
    xmean=xtr.mean(axis=1);ymean=ytr.mean(axis=1)
    xvar=np.mean((xtr-xmean[:,None])**2,axis=1)
    if np.any(xvar<1e-14):raise ValueError("constant feature in sampled columns")
    slope=np.mean((xtr-xmean[:,None])*(ytr-ymean[:,None]),axis=1)/xvar
    intercept=ymean-slope*xmean
    pred=slope[:,None]*p[:,sample_test]+intercept[:,None]
    target=paired[:,sample_test]
    rmse=np.sqrt(np.mean((pred-target)**2,axis=1))
    sst=np.mean((target-target.mean(axis=1,keepdims=True))**2,axis=1)
    valid=sst>1e-14
    if int(valid.sum())<len(valid)*0.9:raise ValueError("source heldout feature variance collapsed")
    heldout_r2=1-rmse[valid]**2/sst[valid]
    heldout_nrmse=rmse[valid]/np.sqrt(sst[valid])
    strict_equal=np.mean(np.abs(p-paired)<=1e-6)
    abs_diff=np.abs(p-paired)
    b_control=float(np.median(shift))
    rank_median=quantile(per_feature_spearman,.5)
    r2_med=quantile(heldout_r2,.5)
    return {
      "source_exact_feature_pairs_n":nfeatures,
      "PDC_and_CDAP_columns_each_n":ncols,
      "independent_feature_training_n":len(train),
      "independent_feature_test_n":len(test),
      "source_column_assignment_train_vs_test_agreement_n":agreement,
      "exact_untransformed_numeric_cell_fraction_1e_minus6":round(float(strict_equal),8),
      "raw_unchanged_cell_absolute_difference_median":quantile(abs_diff,.5),
      "mapped_feature_spearman_q10":quantile(per_feature_spearman,.1),
      "mapped_feature_spearman_median":rank_median,
      "mapped_feature_spearman_q90":quantile(per_feature_spearman,.9),
      "mapped_feature_spearman_above_0_99_fraction":round(float(np.mean(per_feature_spearman>0.99)),8),
      "one_column_rotated_negative_control_rank_median":round(b_control,8),
      "mapped_vs_rotated_rank_median_delta":round(rank_median-b_control,8),
      "independently_heldout_sample_columns_n":len(sample_test),
      "fitted_affine_source_transform_feature_n":len(slope),
      "affine_source_slope_positive_fraction":round(float(np.mean(slope>0)),8),
      "affine_source_slope_q10":quantile(slope,.1),
      "affine_source_slope_median":quantile(slope,.5),
      "affine_source_slope_q90":quantile(slope,.9),
      "affine_heldout_sample_R2_q10":quantile(heldout_r2,.1),
      "affine_heldout_sample_R2_median":r2_med,
      "affine_heldout_sample_R2_q90":quantile(heldout_r2,.9),
      "affine_heldout_sample_normalized_RMSE_median":quantile(heldout_nrmse,.5),
      "num_featurewise_heldout_R2_above_0_99":int(np.sum(heldout_r2>0.99)),
      "numerical_source_transform_technical_GO_candidate":
           agreement==ncols and rank_median>0.99
           and quantile(per_feature_spearman,.1)>0.97
           and r2_med>0.99 and quantile(heldout_r2,.1)>0.97,
      "official_PDC_aliquot_or_case_identity_certified":False,
      "independent_NODE_RNA_biopsy_match_certified":False,
      "source_quant_biological_replicate_claim":False,
      "raw_sample_column_map_or_gene_values_exported":False
    }

def audit():
    out={"pdc_study_id":"PDC000198","status":"ACCESS_HOLD",
         "private_source_identifiers_exported":False,
         "PDC_aliquot_case_join_certified":False,"NODE_RNA_matched":False}
    try:
        urls,source_sha=retry(lambda:original_tmt_reports("PDC000198"))
        out["source_file_catalog_response_SHA256"]=source_sha
        if len(urls)!=1:
            out["status"]="SOURCE_TMT_REPORT_AMBIGUOUS"
            return out
        ph,pr,pdc_sha=retry(pdc_prefix)
        sh,sr,cdap_sha=retry(lambda:report_prefix(urls[0]))
        p,c,features=exact_unique_feature_rows(pr,sr,sh)
        out.update({"PDC_quant_source_prefix_SHA256":pdc_sha,
                    "original_CDAP_TMT_prefix_SHA256":cdap_sha,
                    "source_feature_exact_match_audit":features})
        if len(p)<24:
            out["status"]="INSUFFICIENT_EXACT_FEATURES_FOR_TWO_HOLDOUTS"
            return out
        out.update({"status":"NUMERIC_TRANSFORM_WITH_DOUBLE_HOLDOUT_VALIDATED",
                    "source_numeric_transform_test":evaluate(p,c)})
    except Exception as e:
        out["error_class"]=type(e).__name__
    return out

def verify(out):
    if out.get("schema")!="hcc-stage44-source-featurewise-affine-holdout-v1":
        raise ValueError("schema")
    x=out["study"]
    if any(x.get(k) is not False for k in ("private_source_identifiers_exported",
            "PDC_aliquot_case_join_certified","NODE_RNA_matched")):
        raise ValueError("unapproved identity certification")
    if x["status"]=="NUMERIC_TRANSFORM_WITH_DOUBLE_HOLDOUT_VALIDATED":
        y=x["source_numeric_transform_test"]
        if any(y.get(k) is not False for k in (
            "official_PDC_aliquot_or_case_identity_certified",
            "independent_NODE_RNA_biopsy_match_certified",
            "source_quant_biological_replicate_claim",
            "raw_sample_column_map_or_gene_values_exported")):
            raise ValueError("biological patient mapping overclaimed")
        if y["source_column_assignment_train_vs_test_agreement_n"]>330 or y["PDC_and_CDAP_columns_each_n"]!=330:
            raise ValueError("nonsensical original source mapping counts")
    elif x["status"] not in ("ACCESS_HOLD","SOURCE_TMT_REPORT_AMBIGUOUS",
                              "INSUFFICIENT_EXACT_FEATURES_FOR_TWO_HOLDOUTS"):
        raise ValueError("unexpected source access status")
    return {"schema":"hcc-stage44-original-source-transform-postflight-v1",
      "original_source_featurewise_affine_technical_candidate":
         x.get("source_numeric_transform_test",{}).get("numerical_source_transform_technical_GO_candidate",False),
      "G30_B":"HOLD_OFFICIAL_ALIQUOT_TO_330_PDC_COLUMN_AUTHORIZATION",
      "G30_C":"HOLD_NODE_OEP000321_RNA_SOURCE",
      "G30_D":"STOP_NO_INDEPENDENT_SAME_BIOPSY_RNA_PROTEIN",
      "biological_significance_claim_allowed":False}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",required=True)
    parser.add_argument("--review-out",required=True)
    args=parser.parse_args()
    doc={"schema":"hcc-stage44-source-featurewise-affine-holdout-v1",
         "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
         "study":audit()}
    reviewed=verify(doc)
    for target,obj in ((args.out,doc),(args.review_out,reviewed)):
        p=pathlib.Path(target);p.parent.mkdir(parents=True,exist_ok=True)
        p.write_text(json.dumps(obj,indent=2)+"\n")
    print("HCC_STAGE44_ORIGINAL_TRANSFORM",json.dumps(doc))
    print("HCC_STAGE44_REVIEW",json.dumps(reviewed))

if __name__=="__main__":
    main()
