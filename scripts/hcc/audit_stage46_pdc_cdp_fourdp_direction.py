#!/usr/bin/env python3
"""Stage46: exact four-decimal source rounding direction on independent evidence.

Use source-verified exact feature labels to build paired numeric matrices;
infer all 330 source columns ONLY from 53 feature rows and challenge the
independent other 53 feature rows. Determine whether original PDC numeric
cells equal round(CDAP TMT original cells, 4). Never export original gene or
sample labels, individual values, source URLs, or inferred sample permutation.
Exact subset data equality is NOT a complete source matrix or patient/aliquot
identity validation.
"""
import argparse,collections,datetime as dt,json,pathlib,re
import numpy as np
from scipy.optimize import linear_sum_assignment
from audit_stage40_cdap_crossassay_source_labels import retry,original_tmt_reports
from audit_stage42_cdap_pdc_numeric_prefix import pdc_prefix,report_prefix
from audit_stage43_numeric_rank_holdout import exact_unique_feature_rows,rowwise_rank_z,similarity
from audit_stage45_numeric_precision_holdout import digits_of_decimal

def digit_histogram_from_original_qdata(rows,limit=15):
    out=collections.Counter()
    for row in rows[:limit]:
        for v in row[1:]:
            n=digits_of_decimal(v)
            out["non_numeric" if n is None else str(n)]+=1
    return dict(sorted(out.items()))

def rounding_direction(p,c):
    if p.shape!=c.shape or p.ndim!=2 or p.shape[1]!=330 or p.shape[0]<24:
        raise ValueError("source matrix dimension")
    if not np.isfinite(p).all() or not np.isfinite(c).all():
        raise ValueError("nonfinite source")
    n=p.shape[0]
    train=np.arange(0,n,2);hold=np.arange(1,n,2)
    pr=rowwise_rank_z(p);cr=rowwise_rank_z(c)
    i,j=linear_sum_assignment(-similarity(pr[train],cr[train]))
    i2,j2=linear_sum_assignment(-similarity(pr[hold],cr[hold]))
    mapping=np.empty(330,dtype=int);other=np.empty(330,dtype=int)
    mapping[i]=j;other[i2]=j2
    if len(np.unique(mapping))!=330 or len(np.unique(other))!=330:
        raise ValueError("source numeric map not a bijection")
    x=p[hold]
    y=c[hold][:,mapping]
    forward=np.abs(x-np.round(y,4))
    reverse=np.abs(y-np.round(x,4))
    original=np.abs(x-y)
    shift=np.abs(x-np.round(np.roll(y,1,axis=1),4))
    abs_tol=1e-10
    forward_n=int(np.count_nonzero(forward<=abs_tol))
    reverse_n=int(np.count_nonzero(reverse<=abs_tol))
    negative_n=int(np.count_nonzero(shift<=abs_tol))
    total=int(x.size)
    return {
        "source_original_exact_feature_pairs_n":n,
        "training_original_feature_n":len(train),
        "independent_heldout_original_feature_n":len(hold),
        "original_source_sample_columns_n":330,
        "independently_inferred_column_mappings_same_n":int(np.sum(mapping==other)),
        "independent_heldout_raw_numeric_cells_n":total,
        "PDC_equals_CDAP_rounded_to_4dp_cells_n":forward_n,
        "PDC_equals_CDAP_rounded_to_4dp_fraction":round(forward_n/total,9),
        "CDAP_equals_PDC_rounded_to_4dp_cells_n":reverse_n,
        "CDAP_equals_PDC_rounded_to_4dp_fraction":round(reverse_n/total,9),
        "rotated_negative_control_exact_4dp_cells_n":negative_n,
        "rotated_negative_control_exact_4dp_fraction":round(negative_n/total,9),
        "original_source_unrounded_absolute_difference_max":round(float(original.max()),9),
        "PDC_equals_CDAP_4dp_original_feature_prefix_confirmed":
            forward_n==total and int(np.sum(mapping==other))==330,
        "source_all_feature_rows_verified":False,
        "PDC_official_aliquot_case_identifiers_verified":False,
        "NODE_RNA_biopsy_identity_verified":False,
        "source_values_or_mapping_exported":False}

def probe():
    out={"pdc_study_id":"PDC000198","status":"SOURCE_ACCESS_HOLD",
         "original_source_sample_or_gene_ids_exported":False,
         "full_source_numeric_matrix_verified":False,
         "PDC_official_biospecimen_join_verified":False,
         "independent_RNA_protein_matched":False}
    try:
        files,catalog_sha=retry(lambda:original_tmt_reports("PDC000198"))
        out["source_catalog_SHA256"]=catalog_sha
        if len(files)!=1:
            out["status"]="PDC_ORIGINAL_TMT_FILE_AMBIGUOUS";return out
        ph,pr,pdc_sha=retry(pdc_prefix)
        sh,sr,cdap_sha=retry(lambda:report_prefix(files[0]))
        p,c,features=exact_unique_feature_rows(pr,sr,sh)
        out.update({"PDC_quant_prefix_SHA256":pdc_sha,
          "original_CDAP_TMT_prefix_SHA256":cdap_sha,
          "exact_source_feature_count_audit":features,
          "PDC_first_source_rows_numeric_decimal_digit_counts":
             digit_histogram_from_original_qdata(pr)})
        if len(p)<24:
            out["status"]="INSUFFICIENT_ORIGINAL_SOURCE_FEATURES";return out
        out["independent_source_rounding_direction_evidence"]=rounding_direction(p,c)
        out["status"]="INDEPENDENT_SOURCE_PREFIX_ROUNDING_DIRECTION_VERIFIED"
    except Exception as e:
        out["error_class"]=type(e).__name__
    return out

def validate(data):
    if data.get("schema")!="hcc-stage46-PDC-CDAP-fourdp-direction-v1":
        raise ValueError("source schema")
    source=data["source"]
    if any(source.get(x) is not False for x in (
          "original_source_sample_or_gene_ids_exported",
          "full_source_numeric_matrix_verified",
          "PDC_official_biospecimen_join_verified",
          "independent_RNA_protein_matched")):
        raise ValueError("unsupported scientific assertion")
    if source["status"]=="INDEPENDENT_SOURCE_PREFIX_ROUNDING_DIRECTION_VERIFIED":
        z=source["independent_source_rounding_direction_evidence"]
        if z["independent_heldout_raw_numeric_cells_n"]!=z["independent_heldout_original_feature_n"]*330:
            raise ValueError("invalid source cell count")
        for x in ("source_all_feature_rows_verified",
                  "PDC_official_aliquot_case_identifiers_verified",
                  "NODE_RNA_biopsy_identity_verified",
                  "source_values_or_mapping_exported"):
            if z.get(x) is not False:raise ValueError("source privacy or biology violation")
    elif source["status"] not in ("SOURCE_ACCESS_HOLD","PDC_ORIGINAL_TMT_FILE_AMBIGUOUS",
                                 "INSUFFICIENT_ORIGINAL_SOURCE_FEATURES"):
        raise ValueError("unexpected source status")
    return {"schema":"hcc-stage46-original-fourdp-postflight-v1",
         "source_exact_fourdp_direction_found":
            source.get("independent_source_rounding_direction_evidence",{}).get(
              "PDC_equals_CDAP_4dp_original_feature_prefix_confirmed",False),
         "G30_B":"HOLD_OFFICIAL_PDC_COLUMN_TO_ALIQUOT_CASE_ID",
         "G30_C":"HOLD_NODE_OEP000321_RNA_ORIGINAL_QUANT",
         "G30_D":"STOP_NO_INDEPENDENT_PATIENT_MATCHED_RNA_PROTEIN",
         "biological_or_clinical_claims_allowed":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True);p.add_argument("--review-out",required=True)
    a=p.parse_args()
    doc={"schema":"hcc-stage46-PDC-CDAP-fourdp-direction-v1",
         "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
         "source":probe()}
    result=validate(doc)
    for path,obj in ((a.out,doc),(a.review_out,result)):
        target=pathlib.Path(path);target.parent.mkdir(parents=True,exist_ok=True)
        target.write_text(json.dumps(obj,indent=2)+"\n")
    print("HCC_STAGE46_ROUNDING_DIRECTION",json.dumps(doc))
    print("HCC_STAGE46_REVIEW",json.dumps(result))
if __name__=="__main__":
    main()
