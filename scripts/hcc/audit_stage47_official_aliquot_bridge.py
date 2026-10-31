#!/usr/bin/env python3
"""Stage47: held-out PDC/CDAP quantitative permutation -> official PDC aliquot IDs.

The numeric mapping links *quantitative source columns* to CDAP ratio headers.
The ID bridge is considered supported only when every one of those original
headers exactly names one distinct, official NON-reference aliquot. No ordinal,
fuzzy or similarity-based ID join is accepted. Never print source labels,
individual IDs, numeric values, permutations, or signed download URLs.
"""
import argparse
import collections
import datetime as dt
import json
import pathlib
import re

import numpy as np
from scipy.optimize import linear_sum_assignment

from audit_stage40_cdap_crossassay_source_labels import original_tmt_reports, retry
from audit_stage42_cdap_pdc_numeric_prefix import pdc_prefix, report_prefix
from audit_stage43_numeric_rank_holdout import exact_unique_feature_rows, rowwise_rank_z, similarity
from audit_stage33_original_cdap_sample_bridge import discover, read_bytes, parse_sample_table
from audit_stage39_original_cdap_report_header import expected
from audit_stage31_pdc_tmt_header_templates import run_design
from probe_stage30_gao_catalog import pdc_query

PID="PDC000198"
QUERY='{biospecimenPerStudy(pdc_study_id:"PDC000198"){case_id sample_id aliquot_id aliquot_submitter_id sample_submitter_id}}'
RATIO=re.compile(r"^(.+?) Log Ratio$",re.I)

def numerical_bridge(p,c,source_header):
    if p.shape!=c.shape or p.shape[1]!=330 or p.shape[0]<24:
        raise ValueError("unexpected numeric source axes")
    cols=[m.group(1).strip() for h in source_header if (m:=RATIO.fullmatch(h))]
    if len(cols)!=330 or len(set(cols))!=330:
        raise ValueError("CDAP biological ratio labels not unique")
    a=rowwise_rank_z(p);b=rowwise_rank_z(c)
    train=np.arange(0,p.shape[0],2);hold=np.arange(1,p.shape[0],2)
    i,j=linear_sum_assignment(-similarity(a[train],b[train]))
    k,l=linear_sum_assignment(-similarity(a[hold],b[hold]))
    mapping=np.empty(330,dtype=int);other=np.empty(330,dtype=int)
    mapping[i]=j;other[k]=l
    same=int(np.count_nonzero(mapping==other))
    delta=np.abs(p[hold]-np.round(c[hold][:,mapping],4))
    exact=int(np.count_nonzero(delta<=1e-10))
    if same!=330 or exact!=len(hold)*330:
        raise ValueError("independent numeric source bridge not reproduced")
    return [cols[j] for j in mapping],{
        "pdc_quant_column_n":330,
        "source_numeric_features_n":int(p.shape[0]),
        "independent_heldout_numeric_features_n":int(len(hold)),
        "independent_permutation_agreement_n":same,
        "independent_four_decimal_equal_cells_n":exact,
        "independent_four_decimal_cells_total_n":int(delta.size),
        "numeric_source_bijection_verified":True,
        "numeric_permutation_exported":False,
    }

def alias_evidence(labels, records, reference_ids, designed_aliquots):
    """Only exact official aliquot namespaces can certify source ID labels."""
    if len(labels)!=330 or len(set(labels))!=330:
        raise ValueError("invalid quantified CDAP label universe")
    by_kind={k:collections.defaultdict(set) for k in
             ("aliquot_id","aliquot_submitter_id","sample_id","sample_submitter_id")}
    known=set()
    for r in records:
        if not isinstance(r,dict):continue
        aid=str(r.get("aliquot_id") or "").strip()
        if not aid:continue
        known.add(aid)
        if aid in reference_ids:continue
        for kind,index in by_kind.items():
            value=str(r.get(kind) or "").strip()
            if value:index[value].add(aid)
    report={}
    for kind,index in by_kind.items():
        candidates=[index.get(v,set()) for v in labels]
        unique=[next(iter(s)) for s in candidates if len(s)==1]
        report[kind]={
            "exact_unique_header_matches_n":len(unique),
            "distinct_aliquots_for_exact_matches_n":len(set(unique)),
            "ambiguous_header_matches_n":sum(len(s)>1 for s in candidates),
            "missing_header_matches_n":sum(not s for s in candidates),
            "official_alias_collisions_n":sum(len(s)>1 for s in index.values()),
        }
    strict=report["aliquot_id"]
    submitter=report["aliquot_submitter_id"]
    certified=any(
        m["exact_unique_header_matches_n"]==330
        and m["distinct_aliquots_for_exact_matches_n"]==330
        and m["ambiguous_header_matches_n"]==0
        and len(designed_aliquots-reference_ids)==330
        and set(next(iter(by_kind[k][v])) for v in labels)==designed_aliquots-reference_ids
        for k,m in (("aliquot_id",strict),("aliquot_submitter_id",submitter))
        if all(len(by_kind[k].get(v,set()))==1 for v in labels)
    )
    return {
        "source_official_biospecimen_unique_aliquot_ids_n":len(known),
        "source_official_design_nonreference_aliquots_n":len(designed_aliquots-reference_ids),
        "source_official_reference_aliquots_n":len(reference_ids),
        "alias_namespace_results":report,
        "strict_330_official_aliquot_id_bridge_verified":bool(certified),
        "raw_aliquot_case_sample_ids_or_CDAP_names_exported":False,
    }

def probe():
    result={"pdc_study_id":PID,"status":"SOURCE_ACCESS_HOLD",
            "original_labels_or_values_exported":False,
            "NODE_RNA_same_biopsy_certified":False}
    try:
        files,catalog_digest=retry(lambda:original_tmt_reports(PID))
        if len(files)!=1:
            raise ValueError("original CDAP TMT source not unique")
        _,qrows,pdc_digest=retry(pdc_prefix)
        chead,crows,cdap_digest=retry(lambda:report_prefix(files[0]))
        p,c,feature_meta=exact_unique_feature_rows(qrows,crows,chead)
        labels,bridge=numerical_bridge(p,c,chead)
        companions,_=retry(lambda:discover(PID))
        if len(companions)!=1:raise ValueError("original CDAP sample design not unique")
        original,reporters=parse_sample_table(retry(lambda:read_bytes(companions[0])))
        expected_labels=expected(original,reporters)
        if len(expected_labels)!=330 or set(labels)!=set(expected_labels):
            raise ValueError("original CDAP sample design labels drift")
        design,ref,official_design_digest=retry(lambda:run_design(PID))
        d,biospecimen_digest=retry(lambda:pdc_query(QUERY))
        records=d.get("biospecimenPerStudy")
        if not isinstance(records,list) or not records:
            raise ValueError("official PDC biospecimen list unavailable")
        ali={v["aliquot"] for v in design}
        identity=alias_evidence(labels,records,ref,ali)
        result.update({
            "status":"PDC_NUMERIC_TO_CDAP_LABEL_TO_OFFICIAL_ID_EXACT_AUDITED",
            "PDC_original_quant_prefix_sha256":pdc_digest,
            "CDAP_original_TMT_prefix_sha256":cdap_digest,
            "PDC_source_catalog_sha256":catalog_digest,
            "PDC_official_design_response_sha256":official_design_digest,
            "PDC_biospecimen_response_sha256":biospecimen_digest,
            "original_feature_evidence":feature_meta,
            "independent_numeric_bridge":bridge,
            "official_identity_bridge":identity,
        })
    except Exception as exc:
        result["error_class"]=type(exc).__name__
    return result

def review(obj):
    if obj.get("schema")!="hcc-stage47-official-PDC-aliquot-bridge-v1":
        raise ValueError("unexpected evidence schema")
    source=obj["source"]
    if source.get("original_labels_or_values_exported") is not False or source.get("NODE_RNA_same_biopsy_certified") is not False:
        raise ValueError("unapproved identity disclosure or RNA inference")
    checked=source["status"]=="PDC_NUMERIC_TO_CDAP_LABEL_TO_OFFICIAL_ID_EXACT_AUDITED"
    matched=checked and source["official_identity_bridge"]["strict_330_official_aliquot_id_bridge_verified"]
    return {
        "schema":"hcc-stage47-official-aliquot-postflight-v1",
        "pdc_numeric_to_original_cdap_label_checked":checked,
        "strict_source_330_aliquot_identity_bridge_verified":bool(matched),
        "G30_B":"GO_OFFICIAL_ALIQUOT_BRIDGE" if matched else "HOLD_OFFICIAL_PDC_QUANT_COLUMN_ALIQUOT_ID",
        "G30_C":"HOLD_NODE_OEP000321_RNA_MATRIX_AND_SAMPLE_ID",
        "G30_D":"STOP_NO_SAME_BIOPSY_INDEPENDENT_RNA_PROTEIN",
        "biological_or_clinical_claims_allowed":False,
        "patient_identifiers_exported":False,
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True)
    ap.add_argument("--review-out",required=True)
    args=ap.parse_args()
    data={"schema":"hcc-stage47-official-PDC-aliquot-bridge-v1",
          "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),"source":probe()}
    verdict=review(data)
    for name,item in ((args.out,data),(args.review_out,verdict)):
        path=pathlib.Path(name);path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(item,indent=2)+"\n",encoding="utf-8")
    print("HCC_STAGE47_OFFICIAL_BRIDGE",json.dumps(data))
    print("HCC_STAGE47_REVIEW",json.dumps(verdict))

if __name__=="__main__":
    main()
