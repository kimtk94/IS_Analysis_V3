#!/usr/bin/env python3
"""Stage43: independent held-out numerical pattern test of PDC vs CDAP.

No source features, clinical IDs, aliases, TMT channels, source sample mapping,
or molecular values exported. Hungarian sample assignments are exploratory
numerical hypotheses; NEVER authoritatively join PDC aliquots or NODE RNA.
"""
import argparse,collections,datetime as dt,hashlib,json,math,pathlib,re
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.stats import rankdata
from audit_stage40_cdap_crossassay_source_labels import retry,original_tmt_reports
from audit_stage42_cdap_pdc_numeric_prefix import pdc_prefix,report_prefix,finite

def exact_unique_feature_rows(quant,report,headers):
    """Only unique identical feature identifiers; do NOT guess feature order."""
    if len(quant)<1 or len(headers)<331:raise ValueError("source numeric axes unavailable")
    rr=[i for i,x in enumerate(headers)
        if re.fullmatch(r".+ Log Ratio",x,re.I)
        and not re.search(r"Unshared Log Ratio$",x,re.I)]
    if len(rr)!=330:raise ValueError("CDAP source 330 Log Ratio fields missing")
    nonratio=[i for i,x in enumerate(headers)
        if i not in rr and not re.search(r"Unshared Log Ratio$",x,re.I)]
    counts=collections.Counter(r[0] for r in quant if isinstance(r[0],str))
    qby={r[0]:r for r in quant if isinstance(r[0],str) and counts[r[0]]==1}
    sby=collections.defaultdict(list)
    multi=0
    for row in report:
        hits={row[i].strip() for i in nonratio[:15]
              if isinstance(row[i],str) and row[i].strip() in qby}
        if len(hits)>1:
            multi+=1;continue
        if len(hits)==1:sby[next(iter(hits))].append(row)
    common=sorted(x for x in qby if len(sby.get(x,[]))==1)
    pa=[];sb=[];dropped=0
    for feature in common:
        v=[finite(x) for x in qby[feature][1:]]
        w=[finite(sby[feature][0][i]) for i in rr]
        if len(v)!=330 or len(w)!=330:raise ValueError("source feature width changed")
        va=np.array([np.nan if z is None else z for z in v],dtype=np.float64)
        wa=np.array([np.nan if z is None else z for z in w],dtype=np.float64)
        if np.mean(np.isfinite(va))<0.95 or np.mean(np.isfinite(wa))<0.95:
            dropped+=1;continue
        for a in (va,wa):
            if not np.all(np.isfinite(a)):
                a[~np.isfinite(a)]=np.nanmedian(a)
        if np.std(va)<1e-8 or np.std(wa)<1e-8:
            dropped+=1;continue
        pa.append(va);sb.append(wa)
    return (np.asarray(pa,dtype=np.float64),np.asarray(sb,dtype=np.float64),
            {"unique_exact_source_feature_name_pairs_n":len(common),
             "source_feature_rows_with_multiple_possible_labels_n":multi,
             "source_feature_rows_dropped_for_missingness_or_constancy_n":dropped,
             "retained_complete_feature_pairs_n":len(pa),
             "individual_gene_ids_or_values_exported":False})

def rowwise_rank_z(x):
    a=rankdata(x,axis=1,method="average").astype(np.float64)
    a-=a.mean(axis=1,keepdims=True)
    sd=a.std(axis=1,keepdims=True)
    if np.any(sd<1e-8):raise ValueError("constant source ranked feature")
    return a/sd

def similarity(a,b):
    if a.shape!=b.shape:raise ValueError("source feature matrix axes differ")
    # Each 330-column sample is a vector over disjoint training features.
    # Center/scale over FEATURES to test ranked cross-assay patterns.
    x=a.T.copy();y=b.T.copy()
    x-=x.mean(axis=1,keepdims=True)
    y-=y.mean(axis=1,keepdims=True)
    xn=np.linalg.norm(x,axis=1,keepdims=True)
    yn=np.linalg.norm(y,axis=1,keepdims=True)
    xn[xn<1e-10]=1.0;yn[yn<1e-10]=1.0
    return (x/xn)@(y/yn).T

def score_numeric_mapping(pdc,cdap):
    if pdc.shape!=cdap.shape or pdc.ndim!=2:
        raise ValueError("source pair matrices incompatible")
    genes,n=pdc.shape
    if n<3 or genes<8:raise ValueError("independent feature holdout too small")
    q=rowwise_rank_z(pdc);s=rowwise_rank_z(cdap)
    tr=np.arange(0,genes,2)
    te=np.arange(1,genes,2)
    left=similarity(q[tr],s[tr])
    right=similarity(q[te],s[te])
    i,j=linear_sum_assignment(-left)
    i2,j2=linear_sum_assignment(-right)
    map_train=dict(zip(i.tolist(),j.tolist()))
    map_test=dict(zip(i2.tolist(),j2.tolist()))
    if len(map_train)!=n or len(set(map_train.values()))!=n:
        raise ValueError("source numeric candidate not bijective")
    hold=np.array([right[k,map_train[k]] for k in range(n)])
    baseline=np.array([right[k,(map_train[k]+1)%n] for k in range(n)])
    agreement=sum(map_train[k]==map_test[k] for k in range(n))
    median=float(np.median(hold))
    base=float(np.median(baseline))
    return {
        "source_shared_gene_pairs_n":genes,"source_assay_columns_each_n":n,
        "train_feature_n":len(tr),"strict_disjoint_holdout_feature_n":len(te),
        "numerically_independently_inferred_column_candidate_agreement_n":agreement,
        "heldout_rank_profile_cosine_median":round(median,6),
        "heldout_rank_profile_cosine_q10":round(float(np.quantile(hold,.1)),6),
        "heldout_rank_profile_cosine_fraction_gt_0_8":round(float(np.mean(hold>0.8)),6),
        "heldout_rotated_control_cosine_median":round(base,6),
        "heldout_delta_vs_rotated_control_median":round(median-base,6),
        "provisional_numerical_source_bijection_candidate":
            genes>=24 and n==330 and agreement>=300 and median>0.8 and median-base>0.3,
        "source_sample_position_joins_used":False,
        "inferred_permutation_or_source_names_exported":False,
        "official_PDC_aliquot_case_IDs_certified":False,
        "NODE_RNA_patient_join_certified":False}

def audit():
    out={"pdc_study_id":"PDC000198","status":"ACCESS_HOLD",
         "original_sample_or_gene_ids_exported":False,
         "official_PDC_aliquot_assignment_certified":False,
         "RNA_protein_same_patient_match_certified":False}
    try:
        urls,digest=retry(lambda:original_tmt_reports("PDC000198"))
        out["official_source_catalog_SHA256"]=digest
        if len(urls)!=1:
            out["status"]="SOURCE_TMT_REPORT_AMBIGUOUS";return out
        ph,pr,pdigest=retry(pdc_prefix)
        sh,sr,sdigest=retry(lambda:report_prefix(urls[0]))
        p,c,details=exact_unique_feature_rows(pr,sr,sh)
        out.update({"PDC_quant_prefix_SHA256":pdigest,
                    "CDAP_TMT_prefix_SHA256":sdigest,
                    "exact_feature_source_audit":details})
        if len(p)<24:
            out["status"]="INSUFFICIENT_CROSS_SOURCE_NUMERIC_FEATURES"
            return out
        out["cross_assay_numerical_source_probe"]=score_numeric_mapping(p,c)
        out["status"]="ORIGINAL_SOURCE_HELDOUT_NUMERIC_PATTERN_CHECKED"
    except Exception as e:
        out["error_class"]=type(e).__name__
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--out",required=True)
    ap.add_argument("--review-out",required=True)
    a=ap.parse_args()
    source=audit()
    if (source.get("original_sample_or_gene_ids_exported") is not False or
        source.get("official_PDC_aliquot_assignment_certified") is not False or
        source.get("RNA_protein_same_patient_match_certified") is not False):
        raise ValueError("unsupported inferred biological patient identity")
    data={"schema":"hcc-stage43-source-PDC-CDAP-heldout-numeric-v1",
          "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),"source":source}
    review={"schema":"hcc-stage43-source-heldout-postflight-v1",
        "technical_numeric_holdout_completed":source["status"]=="ORIGINAL_SOURCE_HELDOUT_NUMERIC_PATTERN_CHECKED",
        "G30_B":"HOLD_NO_AUTHORIZED_330_SOURCE_TO_OFFICIAL_ALIQUOT_CASE_ID",
        "G30_C":"HOLD_NODE_RNA_MATRIX",
        "G30_D":"STOP_NO_INDEPENDENT_PAIRED_RNA_PROTEIN_PATIENT",
        "new_biological_or_clinical_claims_allowed":False}
    for path,obj in ((a.out,data),(a.review_out,review)):
        f=pathlib.Path(path);f.parent.mkdir(parents=True,exist_ok=True)
        f.write_text(json.dumps(obj,indent=2)+"\n")
    print("HCC_STAGE43_HELDOUT_SOURCE",json.dumps(data))
    print("HCC_STAGE43_REVIEW",json.dumps(review))
if __name__=="__main__":main()
