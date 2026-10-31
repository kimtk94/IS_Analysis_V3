#!/usr/bin/env python3
"""Stage45: independently heldout raw numeric PDC vs original CDAP parity QC.

The 330-column numerical map is trained ONLY on the odd/even partition of
shared source features. Raw, untransformed values from all 330 mapped columns
of the OTHER heldout feature partition are compared with strict predefined
absolute tolerances. A rounded source number is not a clinical sample ID.
No raw source sample names, feature IDs, signed URLs, per-column mapping or
individual quantitative cell values are exported to public GitHub artifacts.
"""
import argparse,collections,datetime as dt,json,pathlib,re
import numpy as np
from scipy.optimize import linear_sum_assignment
from audit_stage40_cdap_crossassay_source_labels import original_tmt_reports,retry
from audit_stage42_cdap_pdc_numeric_prefix import pdc_prefix,report_prefix
from audit_stage43_numeric_rank_holdout import exact_unique_feature_rows,rowwise_rank_z,similarity

TOLS=(1e-6,5e-6,1e-5,5.1e-5,1e-4,1e-3)
def qtile(a,p):
    return round(float(np.quantile(a,p)),9)

def digits_of_decimal(v):
    text=str(v).strip()
    if not re.fullmatch(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?",text):
        return None
    mantissa=text.split("e")[0].split("E")[0]
    return len(mantissa.split(".",1)[1]) if "." in mantissa else 0

def cdap_precision_histogram(rows,header,max_rows=15):
    cols=[i for i,v in enumerate(header)
          if re.fullmatch(r".+ Log Ratio",v,re.I)
          and not re.search(r"Unshared Log Ratio$",v,re.I)]
    if len(cols)!=330:raise ValueError("CDAP reporter column shape")
    hist=collections.Counter()
    for row in rows[:max_rows]:
        for j in cols:
            d=digits_of_decimal(row[j])
            hist["non_numeric" if d is None else str(d)]+=1
    return dict(sorted(hist.items()))

def heldout_parity(p,c):
    if p.shape!=c.shape or p.ndim!=2 or p.shape[1]!=330 or p.shape[0]<24:
        raise ValueError("insufficient original shared feature matrix")
    if not np.isfinite(p).all() or not np.isfinite(c).all():
        raise ValueError("source contains nonfinite values")
    n=p.shape[0]
    rp=rowwise_rank_z(p);rc=rowwise_rank_z(c)
    tr=np.arange(0,n,2);te=np.arange(1,n,2)
    itr,jtr=linear_sum_assignment(-similarity(rp[tr],rc[tr]))
    ite,jte=linear_sum_assignment(-similarity(rp[te],rc[te]))
    p1=np.empty(330,dtype=int);p2=np.empty(330,dtype=int)
    p1[itr]=jtr;p2[ite]=jte
    if len(set(p1.tolist()))!=330 or len(set(p2.tolist()))!=330:
        raise ValueError("source match ambiguous")
    agreement=int(np.sum(p1==p2))
    # Critical: holdout source feature rows are NOT used in deriving p1.
    heldout_p=p[te]
    heldout_c=c[te][:,p1]
    diff=np.abs(heldout_p-heldout_c)
    control=np.abs(heldout_p-np.roll(heldout_c,1,axis=1))
    thresholds={}
    for tol in TOLS:
        key=f"source_numeric_cells_within_{tol:.6g}_fraction"
        thresholds[key]=round(float(np.mean(diff<=tol)),9)
    feature_max=np.max(diff,axis=1)
    per_col_max=np.max(diff,axis=0)
    source_4dp_within=thresholds["source_numeric_cells_within_5.1e-05_fraction"]
    # Independent holdout + source precision rule. The threshold does not
    # claim identity of all source proteins or a clinical sample crosswalk.
    return {
        "original_source_shared_feature_n":n,
        "source_train_feature_n":len(tr),
        "source_independent_heldout_feature_n":len(te),
        "original_source_sample_columns_each_n":330,
        "independent_inferred_mapping_agreement_n":agreement,
        "raw_numeric_heldout_comparisons_n":int(diff.size),
        "raw_numeric_absolute_difference_median":qtile(diff,.5),
        "raw_numeric_absolute_difference_q90":qtile(diff,.9),
        "raw_numeric_absolute_difference_q99":qtile(diff,.99),
        "raw_numeric_absolute_difference_max":round(float(diff.max()),9),
        "heldout_feature_max_absolute_difference_median":qtile(feature_max,.5),
        "heldout_column_max_absolute_difference_median":qtile(per_col_max,.5),
        "negative_rotated_column_control_absolute_difference_median":qtile(control,.5),
        "negative_control_vs_proper_absolute_difference_median_ratio":
          round(float(np.median(control)/max(np.median(diff),1e-15)),4),
        **thresholds,
        "all_heldout_prefix_cells_agree_at_4dp_rounding_scale":
          agreement==330 and source_4dp_within==1.0,
        "heldout_prefix_4dp_scale_consistency_99_9pct_candidate":
          agreement==330 and source_4dp_within>=0.999,
        "full_quantitation_matrix_rows_exhaustively_verified":False,
        "official_PDC_aliquot_ID_to_columns_verified":False,
        "same_patient_NODE_RNA_link_verified":False,
        "any_source_cell_or_column_mapping_exported":False,
    }

def audit():
    out={"pdc_study_id":"PDC000198","status":"ACCESS_HOLD",
         "raw_sample_feature_or_quant_values_exported":False,
         "official_sample_IDs_joined":False,
         "external_RNA_protein_matched":False}
    try:
        urls,digest=retry(lambda:original_tmt_reports("PDC000198"))
        out["original_assembly_catalog_SHA256"]=digest
        if len(urls)!=1:
            out["status"]="TMT_SOURCE_REPORT_COUNT_REVIEW";return out
        head,qrows,qdigest=retry(pdc_prefix)
        shead,srows,sdigest=retry(lambda:report_prefix(urls[0]))
        p,c,stats=exact_unique_feature_rows(qrows,srows,shead)
        out.update({"PDC_original_quant_prefix_SHA256":qdigest,
           "CDAP_original_TMT_prefix_SHA256":sdigest,
           "source_exact_feature_alignment":stats,
           "CDAP_original_logratio_decimal_digits_first_rows":
               cdap_precision_histogram(srows,shead)})
        if len(p)<24:
            out["status"]="INSUFFICIENT_SOURCE_FEATURES";return out
        out["source_independent_numeric_parity"]=heldout_parity(p,c)
        out["status"]="NUMERIC_PREFIX_HELDOUT_PRECISION_QC_FINISHED"
    except Exception as e:
        out["error_class"]=type(e).__name__
    return out

def verify(document):
    if document.get("schema")!="hcc-stage45-PDC-CDAP-heldout-precision-parity-v1":
        raise ValueError("source schema")
    src=document.get("source",{})
    for k in ("raw_sample_feature_or_quant_values_exported","official_sample_IDs_joined","external_RNA_protein_matched"):
        if src.get(k) is not False:raise ValueError("unverified original source patient claim")
    if src.get("status")=="NUMERIC_PREFIX_HELDOUT_PRECISION_QC_FINISHED":
        stats=src["source_independent_numeric_parity"]
        if (stats["original_source_sample_columns_each_n"]!=330 or
            stats["raw_numeric_heldout_comparisons_n"]!=stats["source_independent_heldout_feature_n"]*330):
            raise ValueError("invalid heldout numeric comparison count")
        for flag in ("full_quantitation_matrix_rows_exhaustively_verified",
                     "official_PDC_aliquot_ID_to_columns_verified",
                     "same_patient_NODE_RNA_link_verified",
                     "any_source_cell_or_column_mapping_exported"):
            if stats.get(flag) is not False:raise ValueError("unsupported biological/source claim")
    elif src.get("status") not in ("ACCESS_HOLD","TMT_SOURCE_REPORT_COUNT_REVIEW",
                                   "INSUFFICIENT_SOURCE_FEATURES"):
        raise ValueError("unexpected source status")
    return {
       "schema":"hcc-stage45-original-PDC-CDAP-precision-postflight-v1",
       "first_original_feature_prefix_4dp_parity_candidate":
           src.get("source_independent_numeric_parity",{}).get(
           "all_heldout_prefix_cells_agree_at_4dp_rounding_scale",False),
       "first_original_feature_prefix_99_9pct_4dp_parity_candidate":
           src.get("source_independent_numeric_parity",{}).get(
           "heldout_prefix_4dp_scale_consistency_99_9pct_candidate",False),
       "G30_B":"HOLD_OFFICIAL_PDC_330_QUANT_COLUMNS_TO_ALIQUOT_ID",
       "G30_C":"HOLD_NODE_OEP000321_RNA_SAMPLE_ORIGINAL",
       "G30_D":"STOP_NO_INDEPENDENT_MATCHED_PATIENT_RNA_PROTEIN",
       "biological_claims_allowed":False}

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True);p.add_argument("--review-out",required=True)
    a=p.parse_args()
    doc={"schema":"hcc-stage45-PDC-CDAP-heldout-precision-parity-v1",
         "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
         "source":audit()}
    r=verify(doc)
    for name,payload in ((a.out,doc),(a.review_out,r)):
        file=pathlib.Path(name);file.parent.mkdir(parents=True,exist_ok=True)
        file.write_text(json.dumps(payload,indent=2)+"\n")
    print("HCC_STAGE45_SOURCE_PRECISION",json.dumps(doc))
    print("HCC_STAGE45_REVIEW",json.dumps(r))
if __name__=="__main__":main()
