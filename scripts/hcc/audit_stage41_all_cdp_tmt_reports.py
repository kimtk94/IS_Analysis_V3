#!/usr/bin/env python3
"""Stage41: independently evaluate EVERY official CDAP TMT report candidate.

Protect source labels, file names, signed URLs, and numerical measurements.
Source biological-name overlap across assays is not an aliquot/case/RNA join.
"""
import argparse,datetime as dt,hashlib,json,pathlib,re
from probe_stage30_gao_catalog import PDC
from audit_stage33_original_cdap_sample_bridge import discover,read_bytes,parse_sample_table
from audit_stage39_original_cdap_report_header import expected,compare,header
from audit_stage40_cdap_crossassay_source_labels import retry,graphql,original_tmt_reports,cross

def feature_kind(cols):
    leading=[x.lower() for x in cols if not re.search(r"\s+(?:unshared )?log ratio$",x,re.I)]
    phospho=any(re.search(r"phospho|site|localization|modification",x) for x in leading)
    prot=any(re.search(r"gene|protein|accession|uniprot",x) for x in leading)
    if phospho and prot:return "MIXED_FEATURE_SCHEMA_REVIEW"
    if phospho:return "PHOSPHO_FEATURE_SCHEMA_CANDIDATE"
    if prot:return "PROTEIN_FEATURE_SCHEMA_CANDIDATE"
    return "FEATURE_SCHEMA_UNSPECIFIED"

def source(pid,query=graphql,read_sample=read_bytes,read_report=header):
    out={"pdc_study_id":pid,"status":"ACCESS_HOLD",
         "source_labels_or_file_names_exported":False,
         "official_aliquot_ids_certified":False,"original_PDC_quant_numeric_parity_certified":False,
         "RNA_to_protein_same_patient_verified":False}
    private=None
    try:
        items,s_sha=discover(pid,graphql=query)
        out["original_sample_companion_candidates_n"]=len(items)
        out["CDAP_file_manifest_SHA256"]=s_sha
        if len(items)!=1:
            out["status"]="AMBIGUOUS_OR_MISSING_SAMPLE_SOURCE";return out,None
        raw=retry(lambda:read_sample(items[0]))
        rows,reporters=parse_sample_table(raw)
        labels=expected(rows,reporters)
        if len(rows)!=33 or len(labels)!=len(set(labels)) or len(labels)!=330:
            raise ValueError("unverified original 330 distinct biological CDAP ratios")
        out["original_sample_SHA256"]=hashlib.sha256(raw).hexdigest()
        out["analytical_plex_n"]=len(rows)
        out["source_biological_ratio_labels_n"]=len(labels)
        # Compare source-name universes even when the number of TMT report
        # candidates differs between assays. This is still not patient pairing.
        private={"biological_labels":set(labels),
                 "analytical_samples":{str(x.get("AnalyticalSample") or "").strip() for x in rows}}
        files,catalog_sha=retry(lambda:original_tmt_reports(pid,fetch=query))
        out["original_protein_assembly_catalog_SHA256"]=catalog_sha
        out["original_TMT_report_candidate_files_n"]=len(files)
        report=[]
        if not 1<=len(files)<=4:
            out["status"]="SOURCE_REPORT_CANDIDATE_COUNT_REVIEW";return out,private
        for url in files:
            try:
                cols,prefix=retry(lambda:read_report(url))
                stat=compare(cols,labels)
                report.append({"status":"ORIGINAL_TMT_HEADER_VERIFIED",
                               "source_report_prefix_SHA256":prefix,
                               "feature_column_classification":feature_kind(cols),
                               **stat})
            except Exception as e:
                report.append({"status":"SOURCE_HEADER_ACCESS_HOLD","error_class":type(e).__name__,
                               "official_biospecimen_join_certified":False,
                               "quantDataMatrix_cells_reconciled":False,
                               "raw_column_names_exported":False})
        out["independent_TMT_candidate_report_checks"]=report
        out["report_candidates_with_330_CDAP_label_source_match_n"]=sum(
            r.get("source_330_CDAP_report_header_candidate") is True for r in report)
        out["status"]="ALL_ORIGINAL_TMT_SOURCE_REPORT_CANDIDATES_AUDITED"
    except Exception as e:
        out["error_class"]=type(e).__name__
    return out,private

def review(d):
    if d.get("schema")!="hcc-stage41-CDAP-all-TMT-original-report-candidates-v1":
        raise ValueError("incorrect source schema")
    arr=d.get("studies")
    if not isinstance(arr,list) or {q.get("pdc_study_id") for q in arr}!=set(PDC):
        raise ValueError("study identity changed")
    for r in arr:
        if any(r.get(k) is not False for k in ("source_labels_or_file_names_exported","official_aliquot_ids_certified",
            "original_PDC_quant_numeric_parity_certified","RNA_to_protein_same_patient_verified")):
            raise ValueError("unverified biological data")
        if r["status"] not in ("ACCESS_HOLD","AMBIGUOUS_OR_MISSING_SAMPLE_SOURCE",
                               "SOURCE_REPORT_CANDIDATE_COUNT_REVIEW",
                               "ALL_ORIGINAL_TMT_SOURCE_REPORT_CANDIDATES_AUDITED"):
            raise ValueError("invalid source state")
        if r["status"]=="ALL_ORIGINAL_TMT_SOURCE_REPORT_CANDIDATES_AUDITED":
            candidates=r["independent_TMT_candidate_report_checks"]
            if len(candidates)!=r["original_TMT_report_candidate_files_n"]:
                raise ValueError("not all original TMT candidates checked")
            for c in candidates:
                if c["status"] not in ("ORIGINAL_TMT_HEADER_VERIFIED","SOURCE_HEADER_ACCESS_HOLD"):
                    raise ValueError("invalid report state")
                if c["official_biospecimen_join_certified"] is not False or c["quantDataMatrix_cells_reconciled"] is not False or c["raw_column_names_exported"] is not False:
                    raise ValueError("unjustified biospecimen alias or numeric parity")
    c=d.get("cross_original_CDAP_sample_universe",{})
    if c.get("RNA_matching_proven") is not False or c.get("cross_assay_original_CDAP_names_are_official_patient_IDs") is not False:
        raise ValueError("case-ID inference not justified")
    return {"schema":"hcc-stage41-all-TMT-source-postflight-v1",
            "original_TMT_reports_with_330_CDAP_source_ratio_label_match_n":sum(
                s.get("report_candidates_with_330_CDAP_label_source_match_n",0) for s in arr),
            "cross_assay_original_name_intersection_source_checked":
                c.get("status")=="CROSS_ASSAY_ORIGINAL_CDAP_NAMES_COMPARED",
            "G30_B":"HOLD_PDC_QUANT_INTERNAL_COLUMN_TO_OFFICIAL_ALIQUOT_MAP",
            "G30_C":"HOLD_ORIGINAL_NODE_RNA_SAMPLE_MANIFEST",
            "G30_D":"STOP_NO_SOURCE_AUTHORIZED_MATCHED_RNA_PROTEIN",
            "novel_biology_claims_allowed":False}

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True);p.add_argument("--review-out",required=True)
    args=p.parse_args()
    results=[source(pid) for pid in PDC]
    d={"schema":"hcc-stage41-CDAP-all-TMT-original-report-candidates-v1",
       "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
       "studies":[a for a,b in results],
       "cross_original_CDAP_sample_universe":cross(results[0][1],results[1][1])}
    verdict=review(d)
    for loc,obj in ((args.out,d),(args.review_out,verdict)):
        f=pathlib.Path(loc);f.parent.mkdir(parents=True,exist_ok=True)
        f.write_text(json.dumps(obj,indent=2)+"\n")
    print("HCC_STAGE41_ORIGINAL_SOURCE",json.dumps(d))
    print("HCC_STAGE41_REVIEW",json.dumps(verdict))
if __name__=="__main__":
    main()
