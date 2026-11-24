#!/usr/bin/env python3
"""Stage30 direct PDC log2_ratio cross-study header internal alias QC.

Reports ONLY exact overlap counts between public PDC198 proteome and PDC199
phosphoproteome source columns. No raw labels or per-label hashes. Same text
across assay quant headers is NOT biological patient ID authentication.
"""
import argparse,datetime as dt,json,pathlib
from probe_stage30_gao_catalog import PDC
from probe_stage30_gao_quant_api_header_to_biospecimen import get_quant_header

def compare(a,b):
    s1=set(a);s2=set(b)
    if len(s1)!=len(a) or len(s2)!=len(b):raise ValueError("duplicate source header")
    intersection=s1&s2
    return {
        "proteome_column_name_n":len(a),
        "phosphoproteome_column_name_n":len(b),
        "cross_study_common_source_header_token_n":len(intersection),
        "proteome_only_header_token_n":len(s1-s2),
        "phosphoproteome_only_header_token_n":len(s2-s1),
        "identical_source_header_order_position_n":sum(x==y for x,y in zip(a,b)),
        "common_header_token_jaccard":round(len(intersection)/len(s1|s2),9) if s1|s2 else None,
        "this_is_not_official_aliquot_ID_match":True,
        "patient_header_alias_identifiers_exported":False,
    }

def audit():
    source={}
    alias={}
    for pid in PDC:
        try:
            info,head=get_quant_header(pid)
            source[pid]={
               "status":"VERIFIED_BOUNDED_DIRECT_PDC_QUANT_HEADER",
               "candidate_column_n":len(head),
               "first_header_gene_label_present":info["first_header_gene_label_present"],
               "anonymous_quant_numeric_cells_from_first_four_rows_n":
                  info["bounded_postheader_numeric_cells_n"],
               "read_byte_cap":info["read_capped_at_bytes"],
               "aliquot_header_labels_exported":False}
            alias[pid]=head
        except Exception as exc:
            source[pid]={"status":"SOURCE_HEADER_ACCESS_HOLD",
                         "failure_class":type(exc).__name__,
                         "aliquot_header_labels_exported":False}
    both={"status":"CROSS_HEADER_NOT_VERIFIED",
          "patient_header_alias_identifiers_exported":False}
    if len(alias)==2:
        both={"status":"CROSS_HEADER_TOKEN_COMPARISON_VERIFIED",
              **compare(alias[PDC[0]],alias[PDC[1]])}
    return {"schema":"hcc-stage30-gao-direct-quant-crossheader-v1",
        "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
        "studies":source,"comparison":both,
        "official_patient_level_ID_pairing_proven":False,
        "matched_RNA_protein_matrix_available":False,
        "independent_biological_replication_completed":False}

def review(data):
    if data.get("schema")!="hcc-stage30-gao-direct-quant-crossheader-v1" or set(data.get("studies",{}))!=set(PDC):
        raise ValueError("source schema drift")
    for study in data["studies"].values():
        if study.get("status") not in ("VERIFIED_BOUNDED_DIRECT_PDC_QUANT_HEADER","SOURCE_HEADER_ACCESS_HOLD"):
            raise ValueError("invalid source status")
        if study.get("aliquot_header_labels_exported") is not False:raise ValueError("identifier leak")
    cmp=data["comparison"]
    if cmp.get("patient_header_alias_identifiers_exported") is not False:raise ValueError("identifier leak")
    if cmp["status"]=="CROSS_HEADER_TOKEN_COMPARISON_VERIFIED":
        n1=cmp["proteome_column_name_n"];n2=cmp["phosphoproteome_column_name_n"]
        if not all(type(q) is int and q>=0 for q in (n1,n2,cmp["cross_study_common_source_header_token_n"])):
            raise ValueError("invalid quantitative header counts")
        if cmp["cross_study_common_source_header_token_n"]>min(n1,n2):raise ValueError("overlap impossible")
        if cmp["identical_source_header_order_position_n"]>min(n1,n2):raise ValueError("identical position impossible")
    elif cmp["status"]!="CROSS_HEADER_NOT_VERIFIED":raise ValueError("unexpected comparison status")
    for key in ("official_patient_level_ID_pairing_proven","matched_RNA_protein_matrix_available",
                "independent_biological_replication_completed"):
        if data.get(key) is not False:raise ValueError("premature biological claim")
    return {"schema":"hcc-stage30-direct-PDC-cross-header-review-v1",
        "same_internal_quant_header_name_set_tested":cmp["status"]=="CROSS_HEADER_TOKEN_COMPARISON_VERIFIED",
        "shared_internal_source_header_alias_n":cmp.get("cross_study_common_source_header_token_n"),
        "G30_B":"PARTIAL_SOURCE_QUANT_HEADER_SHARED_ALIAS_ONLY",
        "G30_C":"HOLD_ACTUAL_NODE_RNA_MATRIX",
        "G30_D":"STOP_NO_VALIDATED_CASE_RNA_PROTEOME_JOIN",
        "source_aliquot_mapping_complete":False,
        "independent_HCC_biology_replication_pass":False}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True);ap.add_argument("--review-out",required=True)
    opts=ap.parse_args();obj=audit()
    f=pathlib.Path(opts.out);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(obj,indent=2)+"\n")
    result=review(obj)
    target=pathlib.Path(opts.review_out);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
    print("PDC_DIRECT_CROSS_HEADER_AGGREGATES",json.dumps(obj["comparison"]))
if __name__=="__main__":main()
