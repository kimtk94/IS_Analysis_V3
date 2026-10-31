#!/usr/bin/env python3
"""Stage39: match source CDAP TMT report header to source sample.txt ratios.
Only aggregate counts/digests exported. No per-sample IDs, URLs or data values.
"""
import argparse,datetime as dt,hashlib,json,pathlib,re,urllib.request
from probe_stage30_gao_catalog import PDC,pdc_query
from probe_stage30_gao_quant_header import safe_public_https
from audit_stage33_original_cdap_sample_bridge import discover,read_bytes,parse_sample_table,RATIO_PATTERN

LIMIT=262144
REPORT_PATTERN=re.compile(r"^(.*?)\\s+(Unshared )?Log Ratio$",re.I)

def expected(rows,reporters):
    names=[]
    for row in rows:
        for expr in (row.get("Ratios") or "").split(","):
            expr=expr.strip()
            if not expr:continue
            match=RATIO_PATTERN.fullmatch(expr)
            if not match:raise ValueError("unparsed CDAP ratio")
            channels=[c for c in reporters if c.upper().replace("TMT_","").replace("TMT-","")==match.group(1).upper()]
            if len(channels)!=1:raise ValueError("ambiguous TMT reporter")
            name=(row.get(channels[0]) or "").strip()
            if not name or name.upper()=="POOL":raise ValueError("unexpected reference numerator")
            names.append(name)
    return names

def compare(header,labels):
    found=set()
    unshared=set()
    for v in header:
        m=REPORT_PATTERN.fullmatch(v)
        if m:
            (unshared if m.group(2) else found).add(m.group(1).strip())
    distinct=set(labels)
    return {
        "source_expected_biological_ratio_labels_n":len(labels),
        "source_distinct_biological_labels_n":len(distinct),
        "quant_report_header_columns_n":len(header),
        "quant_report_log_ratio_labels_n":len(found),
        "quant_report_unshared_log_ratio_labels_n":len(unshared),
        "exact_source_CDAP_label_in_logratio_header_n":len(found&distinct),
        "missing_expected_labels_n":len(distinct-found),
        "unexpected_report_labels_n":len(found-distinct),
        "source_330_CDAP_report_header_candidate":len(labels)==len(distinct)==len(found)==330 and found==distinct,
        "official_biospecimen_join_certified":False,
        "quantDataMatrix_cells_reconciled":False,
        "raw_column_names_exported":False,
    }

def source_file(pid):
    q=('{filesPerStudy(pdc_study_id:'+json.dumps(pid)+',data_category:"Protein Assembly",offset:0,limit:30)'
       '{file_name data_category file_size signedUrl{url}}}')
    data,sha=pdc_query(q)
    arr=data.get("filesPerStudy")
    if not isinstance(arr,list) or any(x.get("data_category")!="Protein Assembly" for x in arr):
        raise ValueError("unverified Protein Assembly")
    files=[]
    for r in arr:
        name=str(r.get("file_name") or "").lower()
        if not re.search(r"\\.(?:tmt\\d*|itraq)\\.tsv$",name):continue
        links=r.get("signedUrl") or []
        if isinstance(links,dict):links=[links]
        u=next((x.get("url") for x in links if isinstance(x,dict) and safe_public_https(x.get("url"))),None)
        if u:files.append(u)
    return files,sha

def header(url):
    if not safe_public_https(url):raise ValueError("unsafe source URL")
    req=urllib.request.Request(url,headers={"Range":"bytes=0-262143","User-Agent":"HCC-stage39-CDAP-TMT-header"})
    with urllib.request.urlopen(req,timeout=35) as resp:
        data=resp.read(LIMIT+1)
    if len(data)>LIMIT or b"\\n" not in data:raise ValueError("bounded header incomplete")
    first=data.split(b"\\n",1)[0].rstrip(b"\\r")
    if len(first)<20 or b"\\x00" in first:raise ValueError("not source TSV header")
    cols=[x.strip() for x in first.decode("utf-8-sig").split("\\t")]
    if not 10<len(cols)<=1500 or len(set(cols))!=len(cols):
        raise ValueError("invalid source header")
    return cols,hashlib.sha256(data).hexdigest()

def one(pid):
    out={"pdc_study_id":pid,"status":"ACCESS_HOLD",
         "raw_sample_names_or_source_urls_exported":False,
         "source_columns_to_official_aliquots_proven":False,
         "same_biopsy_RNA_protein_proven":False}
    try:
        files,sha=source_file(pid)
        out.update({"source_assembly_catalog_sha256":sha,"source_TMT_report_candidates_n":len(files)})
        companions,_=discover(pid)
        if len(companions)!=1:raise ValueError("CDAP sample design not unique")
        raw=read_bytes(companions[0])
        rows,reporters=parse_sample_table(raw)
        labels=expected(rows,reporters)
        out.update({"CDAP_sample_sha256":hashlib.sha256(raw).hexdigest(),
                    "CDAP_ratio_labels_n":len(labels)})
        checks=[]
        for url in files[:4]:
            try:
                cols,digest=header(url)
                checks.append({"status":"SOURCE_TMT_HEADER_COMPARED",
                               "source_report_prefix_sha256":digest,**compare(cols,labels)})
            except Exception as e:
                checks.append({"status":"SOURCE_HEADER_ACCESS_HOLD","error_class":type(e).__name__,
                               "raw_column_names_exported":False,
                               "official_biospecimen_join_certified":False,
                               "quantDataMatrix_cells_reconciled":False})
        out.update({"status":"SOURCE_TMT_REPORTS_CHECKED","reports":checks})
    except Exception as e:
        out["error_class"]=type(e).__name__
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True);p.add_argument("--review-out",required=True)
    a=p.parse_args()
    obj={"schema":"hcc-stage39-original-CDAP-report-header-v1",
         "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
         "studies":[one(pid) for pid in PDC]}
    for x in obj["studies"]:
        if x["raw_sample_names_or_source_urls_exported"] or x["source_columns_to_official_aliquots_proven"] or x["same_biopsy_RNA_protein_proven"]:
            raise ValueError("unsafe identity assertion")
    verdict={"schema":"hcc-stage39-CDAP-source-header-postflight-v1",
             "same_CDAP_report_header_candidates_n":sum(y.get("source_330_CDAP_report_header_candidate",False) for x in obj["studies"] for y in x.get("reports",[])),
             "G30_B":"HOLD_UNVERIFIED_QUANT_COLUMNS_TO_OFFICIAL_ALIQUOTS",
             "G30_C":"HOLD_NODE_RNA_MATRIX","G30_D":"STOP_NO_MATCHED_RNA_PROTEIN",
             "biological_claims_allowed":False}
    for f,data in ((a.out,obj),(a.review_out,verdict)):
        path=pathlib.Path(f);path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(data,indent=2)+"\\n")
    print("HCC_STAGE39_ORIGINAL_SOURCE",json.dumps(obj))
    print("HCC_STAGE39_REVIEW",json.dumps(verdict))
if __name__=="__main__":main()
