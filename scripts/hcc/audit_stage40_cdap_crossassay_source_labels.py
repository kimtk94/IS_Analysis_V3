#!/usr/bin/env python3
"""Stage40: original CDAP 330 sample-name cross-assay comparison.

Original source CDAP sample.txt plus protein/phosphoprotein TMT report header.
No positional joins, no probabilistic match, no raw file/aliquot identifiers.
Names shared between assays are candidates, NOT an authoritative patient join.
"""
import argparse,datetime as dt,hashlib,json,pathlib,time,urllib.error
from probe_stage30_gao_catalog import PDC,pdc_query
from probe_stage30_gao_quant_header import safe_public_https
from audit_stage33_original_cdap_sample_bridge import discover,read_bytes,parse_sample_table
from audit_stage39_original_cdap_report_header import expected,compare,header
import re

RETRY_CODES=(429,500,502,503,504)
def retry(fn,attempts=3,sleep=time.sleep):
    for n in range(attempts):
        try:return fn()
        except urllib.error.HTTPError as e:
            if e.code not in RETRY_CODES or n+1==attempts:raise
        except (TimeoutError,urllib.error.URLError,ConnectionError):
            if n+1==attempts:raise
        sleep(n+1)

def graphql(q):
    return retry(lambda:pdc_query(q))

def original_tmt_reports(pid,fetch=graphql):
    q=('{filesPerStudy(pdc_study_id:'+json.dumps(pid)+',data_category:"Protein Assembly",offset:0,limit:30)'
        '{file_name file_format data_category signedUrl{url}}}')
    d,digest=fetch(q)
    files=d.get("filesPerStudy")
    if not isinstance(files,list) or len(files)>30 or any(
            not isinstance(v,dict) or v.get("data_category")!="Protein Assembly" for v in files):
        raise ValueError("unverified PDC Protein Assembly catalog")
    urls=[]
    for f in files:
        if not re.search(r"\.(?:tmt\d*|itraq)\.tsv(?:\.gz)?$",str(f.get("file_name") or ""),re.I):
            continue
        signed=f.get("signedUrl") or []
        if isinstance(signed,dict):signed=[signed]
        url=next((x.get("url") for x in signed if isinstance(x,dict) and safe_public_https(x.get("url"))),None)
        if url:urls.append(url)
    return urls,digest

def study(pid,fetch=graphql,download_sample=read_bytes,read_header=header):
    x={"pdc_study_id":pid,"status":"ACCESS_HOLD",
       "source_names_or_urls_exported":False,
       "cross_assay_patient_identity_certified":False,
       "PDC_quant_matrix_to_aliquot_certified":False,
       "external_same_patient_RNA_protein_proven":False}
    local=None
    try:
        found,sha=discover(pid,graphql=fetch)
        x["source_CDAP_catalog_SHA256"]=sha
        x["original_CDAP_sample_companion_candidates_n"]=len(found)
        if len(found)!=1:
            x["status"]="AMBIGUOUS_CDAP_SAMPLE_COMPANION"
            return x,None
        raw=retry(lambda:download_sample(found[0]))
        x["source_CDAP_sample_file_SHA256"]=hashlib.sha256(raw).hexdigest()
        rows,reporters=parse_sample_table(raw)
        labels=expected(rows,reporters)
        if len(rows)!=33 or len(labels)!=330 or len(set(labels))!=330:
            raise ValueError("original CDAP expected 33 plexes and 330 unique biological labels")
        urls,src_sha=original_tmt_reports(pid,fetch)
        x.update({"source_assembly_catalog_SHA256":src_sha,
                  "CDAP_original_plex_rows_n":len(rows),
                  "CDAP_source_ratio_biological_labels_n":len(labels),
                  "original_TMT_quant_report_candidates_n":len(urls)})
        if len(urls)!=1:
            x["status"]="AMBIGUOUS_OR_MISSING_TMT_REPORT"
            return x,None
        cols,digest=retry(lambda:read_header(urls[0]))
        stats=compare(cols,labels)
        x.update({"original_TMT_report_header_prefix_SHA256":digest,
                  "report_header_source_exact_match":stats,
                  "status":"ORIGINAL_CDAP_REPORT_AND_SAMPLE_HEADERS_CHECKED"})
        if stats["source_330_CDAP_report_header_candidate"]:
            local={"biological_labels":set(labels),
                   "analytical_samples":{str(r.get("AnalyticalSample") or "").strip() for r in rows}}
    except Exception as e:
        x["error_class"]=type(e).__name__
    return x,local

def cross(a,b):
    out={"status":"SOURCE_CROSS_ASSAY_HOLD",
         "source_label_sets_shared_no_individual_names_exported":True,
         "cross_assay_original_CDAP_names_are_official_patient_IDs":False,
         "distinct_cases_proven":False,"RNA_matching_proven":False}
    if a is None or b is None:return out
    aset,bset=a["biological_labels"],b["biological_labels"]
    same=aset&bset
    af={x.casefold() for x in aset};bf={x.casefold() for x in bset}
    ana=a["analytical_samples"]&b["analytical_samples"]
    out.update({"status":"CROSS_ASSAY_ORIGINAL_CDAP_NAMES_COMPARED",
        "proteome_biological_labels_n":len(aset),
        "phosphoproteome_biological_labels_n":len(bset),
        "source_exact_matching_biological_names_n":len(same),
        "source_casefold_matching_biological_names_n":len(af&bf),
        "source_CDAP_analytical_sample_exact_name_overlap_n":len(ana),
        "source_same_330_label_universe_candidate":len(aset)==len(bset)==len(same)==330,
        "official_PDC_aliquot_or_case_join_certified":False})
    return out

def review(d):
    if d.get("schema")!="hcc-stage40-Gao-original-CDAP-crossassay-names-v1":
        raise ValueError("invalid evidence schema")
    s=d.get("studies")
    if not isinstance(s,list) or {v.get("pdc_study_id") for v in s}!=set(PDC):
        raise ValueError("source study drift")
    for v in s:
        if any(v.get(k) is not False for k in ("source_names_or_urls_exported","cross_assay_patient_identity_certified",
              "PDC_quant_matrix_to_aliquot_certified","external_same_patient_RNA_protein_proven")):
            raise ValueError("unverified source patient join")
        if v["status"]=="ORIGINAL_CDAP_REPORT_AND_SAMPLE_HEADERS_CHECKED":
            q=v["report_header_source_exact_match"]
            if q["source_expected_biological_ratio_labels_n"]!=330:
                raise ValueError("wrong source CDAP count")
            if q["official_biospecimen_join_certified"] or q["quantDataMatrix_cells_reconciled"]:
                raise ValueError("CDAP source ≠ official PDC aliquot mapping")
        elif v["status"] not in ("ACCESS_HOLD","AMBIGUOUS_CDAP_SAMPLE_COMPANION",
                                  "AMBIGUOUS_OR_MISSING_TMT_REPORT"):
            raise ValueError("unknown source status")
    c=d.get("cross_assay",{})
    if c.get("source_label_sets_shared_no_individual_names_exported") is not True or c.get("RNA_matching_proven") is not False:
        raise ValueError("invalid cross-assay evidence")
    return {"schema":"hcc-stage40-Gao-crossassay-source-postflight-v1",
        "studies_with_330_original_CDAP_TMT_report_labels_n":sum(
            v["status"]=="ORIGINAL_CDAP_REPORT_AND_SAMPLE_HEADERS_CHECKED"
            and v["report_header_source_exact_match"]["source_330_CDAP_report_header_candidate"]
            for v in s),
        "original_CDAP_both_assays_exact_source_name_comparison_completed":
            c.get("status")=="CROSS_ASSAY_ORIGINAL_CDAP_NAMES_COMPARED",
        "G30_B":"HOLD_UNTIL_PDC_QUANT_MATRIX_COLUMN_TO_OFFICIAL_ALIQUOT_JOIN",
        "G30_C":"HOLD_NODE_OEP000321_RNA_ORIGINAL_MATRIX",
        "G30_D":"STOP_NO_SOURCE_AUTHORIZED_RNA_PROTEIN_SAME_PATIENT_JOIN",
        "biological_or_clinical_claims_allowed":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True);p.add_argument("--review-out",required=True)
    a=p.parse_args()
    results=[study(pid) for pid in PDC]
    obj={"schema":"hcc-stage40-Gao-original-CDAP-crossassay-names-v1",
         "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
         "studies":[v[0] for v in results],
         "cross_assay":cross(results[0][1],results[1][1])}
    d=review(obj)
    for target,content in ((a.out,obj),(a.review_out,d)):
        file=pathlib.Path(target);file.parent.mkdir(parents=True,exist_ok=True)
        file.write_text(json.dumps(content,indent=2)+"\n")
    print("HCC_STAGE40_SOURCE",json.dumps(obj))
    print("HCC_STAGE40_REVIEW",json.dumps(d))
if __name__=="__main__":
    main()
