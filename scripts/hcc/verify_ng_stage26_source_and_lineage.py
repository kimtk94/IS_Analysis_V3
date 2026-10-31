#!/usr/bin/env python3
"""Stage26 aggregate-only source evidence QA; biological/source paper claims HOLD."""
import argparse,json,math,pathlib
SHA="461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a"
COMMIT="438f52271d0a86a6701baa7a96e80c516ba575c2"
def validate_lineage(r):
    if r.get("schema")!="hcc-stage26-full-vs-fullz-lineage-v1" or r.get("source_commit")!=COMMIT or r.get("full_sha256")!=SHA:
        raise ValueError("wrong Suter immutable source")
    layers=r.get("layer_statistics")
    if not isinstance(layers,list) or {a.get("modality") for a in layers}!={"M","CN","T","P","PP"}:
        raise ValueError("layer provenance wrong")
    if len(r.get("fullz_sha256",""))!=64:
        raise ValueError("source fullz not sha-pinned")
    for q in layers:
        if q.get("samples")!=50 or q.get("row_and_column_names_identical") is not True:
            raise ValueError("no matching 50 source biopsies")
        if type(q.get("features"))!=int or q["features"]<20:
            raise ValueError("missing features")
        if q.get("both_finite_n",0)<1000:
            raise ValueError("no source quantitative coverage")
        for k in ("row_z_frac_abs_lt_1e_5","col_z_frac_abs_lt_1e_5"):
            value=q.get(k)
            if type(value) not in (int,float) or not math.isfinite(value) or not 0<=value<=1:
                raise ValueError("bad z residual fraction")
        for k in ("row_z_median_abs_residual","col_z_median_abs_residual"):
            value=q.get(k)
            if type(value) not in (int,float) or not math.isfinite(value) or value<0:
                raise ValueError("bad residual")
    if any(r.get(k) is not False for k in (
        "original_ng_RSEM_TMM_edgeR_reproduced",
        "original_ng_5481_gene_universe_reproduced",
        "source_patient_and_gene_identifiers_exported")):
        raise ValueError("forbidden source or biological claim")
    return {"schema":"hcc-stage26-lineage-postflight-v1",
       "source_immutable_SHA_verified":True,
       "matching_T_P_PP_feature_axes_verified":True,
       "source_z_matrix_numeric_relation_screened":True,
       "layer_diagnostic":[{"modality":x["modality"],
           "source_row_z_standardization_matches_abs_1e_5_95pct":x["row_z_frac_abs_lt_1e_5"]>=0.95,
           "row_z_median_abs_residual":x["row_z_median_abs_residual"],
           "col_z_median_abs_residual":x["col_z_median_abs_residual"]}
           for x in layers],
       "original_RDS_generation_code_and_normal_controls_proven":False,
       "V03":"PARTIAL"}
def validate_paper(r):
    if r.get("schema")!="hcc-stage26-ng-paper-gene-source-probe-v1" or r.get("paper_doi")!="10.1038/s41467-022-29960-8":
        raise ValueError("original source paper identity wrong")
    src=r.get("sources")
    if not isinstance(src,dict) or set(src)!={"publisher_Source_Data","publisher_Supp_Data_4","publisher_Supp_Data_2"}:
        raise ValueError("missing public supplemental source probes")
    for key,rec in src.items():
        if rec.get("status")=="PUBLIC_XLSX_CONTENT_AUDITED":
            if not rec.get("source_url","").startswith(("https://media.springernature.com/","https://pmc.ncbi.nlm.nih.gov/")):
                raise ValueError("not official public source")
            if len(rec.get("sha256",""))!=64:
                raise ValueError("no file digest")
            if type(rec.get("sheet_count_inspected"))!=int or rec["sheet_count_inspected"]<1:
                raise ValueError("no XLSX source worksheets")
            for q in rec["sheets"]:
                if q["status"]=="SCAN_COMPLETED":
                    if type(q["unique_ENSG_like_identifiers_n"])!=int:
                        raise ValueError("missing counted IDs")
                    if q["candidate_exact_5481_count_only"]!=(q["unique_ENSG_like_identifiers_n"]==5481):
                        raise ValueError("unsupported candidate cardinality")
        elif rec.get("status")!="UNVERIFIED":
            raise ValueError("unrecognized source probe result")
    if r.get("exact_published_5481_gene_membership_verified") is not False or r.get("patient_gene_and_protein_identifiers_exported") is not False or r.get("raw_original_XLSX_in_artifact") is not False:
        raise ValueError("prohibited exact original gene membership or public identifiers")
    return {"schema":"hcc-stage26-ng-original-genes-postflight-v1",
       "publisher_xlsx_verified_count":sum(x.get("status")=="PUBLIC_XLSX_CONTENT_AUDITED" for x in src.values()),
       "per_source_status":{k:v["status"] for k,v in src.items()},
       "source_gene_id_count_candidates_are_not_identity_replication":True,
       "original_published_5481_genes_exact_membership":"HOLD",
       "source_2022_original_annotation_identity":"HOLD","V03":"PARTIAL","V07":"PARTIAL"}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--kind",choices=("lineage","paper"),required=True)
    p.add_argument("--input",required=True)
    p.add_argument("--out",required=True)
    a=p.parse_args()
    x=json.loads(pathlib.Path(a.input).read_text(encoding="utf-8"))
    y=validate_lineage(x) if a.kind=="lineage" else validate_paper(x)
    dest=pathlib.Path(a.out);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(y,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(y,indent=2))
if __name__=="__main__":main()
