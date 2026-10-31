#!/usr/bin/env python3
"""Stage42: bounded ORIGINAL PDC vs CDAP TMT numeric source parity pilot.

Only exact matching feature names and independently matching 4+ feature ratios
can support a provisional column fingerprint. Numerical coincidence NEVER
certifies official aliquot/case IDs or matching NODE OEP000321 RNA biopsies.
No source gene labels, sample aliases, cells or signed URLs are exported.
"""
import argparse,collections,csv,datetime as dt,hashlib,io,json,math,pathlib,re,urllib.request
from audit_stage40_cdap_crossassay_source_labels import original_tmt_reports,retry
from probe_stage30_gao_direct_quant_endpoint import ENDPOINT
from probe_stage30_gao_quant_header import safe_public_https

PDC_LIMIT=2*1024*1024
CDAP_LIMIT=6*1024*1024
FEATURE_LIMIT=120
SOURCE_ROW_LIMIT=3500

def pdc_prefix(pid="PDC000198",opener=urllib.request.urlopen):
    if pid!="PDC000198":raise ValueError("only documented Proteome source pilot")
    query='{quantDataMatrix(pdc_study_id:"PDC000198",data_type:"log2_ratio")}'
    req=urllib.request.Request(ENDPOINT,data=json.dumps({"query":query}).encode(),
          headers={"Content-Type":"application/json","User-Agent":"HCC-Stage42-Numeric-Prefix-QC/1"},method="POST")
    with opener(req,timeout=80) as rsp:
        raw=rsp.read(PDC_LIMIT+1)
    text=raw[:PDC_LIMIT].decode("utf-8","replace")
    m=re.search(r'"quantDataMatrix"\s*:\s*\[',text[:100000])
    if not m:raise ValueError("PDC source quant matrix marker absent")
    pos=m.end();decoder=json.JSONDecoder();rows=[]
    for _ in range(FEATURE_LIMIT+1):
        while pos<len(text) and text[pos] in " \n\r\t,":pos+=1
        try: v,end=decoder.raw_decode(text,pos)
        except (ValueError,IndexError):break
        if not isinstance(v,list) or len(v)!=331:raise ValueError("PDC row shape drift")
        rows.append(v);pos=end
    if len(rows)<5 or len(rows[0])!=331 or len(set(rows[0][1:]))!=330:
        raise ValueError("PDC source first header + 4 rows not available")
    return rows[0],rows[1:],hashlib.sha256(raw[:PDC_LIMIT]).hexdigest()

def report_prefix(url,opener=urllib.request.urlopen):
    if not safe_public_https(url):raise ValueError("source URL not permitted")
    req=urllib.request.Request(url,headers={"Range":f"bytes=0-{CDAP_LIMIT-1}",
        "User-Agent":"HCC-Stage42-original-CDAP-numeric-prefix/1"})
    with opener(req,timeout=70) as rsp:raw=rsp.read(CDAP_LIMIT+1)
    if len(raw)>CDAP_LIMIT:raise ValueError("CDAP report byte cap")
    decoded=raw.decode("utf-8-sig","replace")
    if "\n" not in decoded:raise ValueError("no complete CDAP report header")
    # Never parse final potentially incomplete source record.
    complete=decoded[:decoded.rfind("\n")+1]
    data=csv.reader(io.StringIO(complete),delimiter="\t")
    cols=next(data)
    if not 10<len(cols)<=1500 or len(cols)!=len(set(cols)):
        raise ValueError("CDAP source header drift")
    rows=[]
    for r in data:
        if len(r)!=len(cols):raise ValueError("nonrectangular original report prefix")
        rows.append(r)
        if len(rows)>=SOURCE_ROW_LIMIT:break
    return cols,rows,hashlib.sha256(raw).hexdigest()

def finite(x):
    if x is None or x=="":return None
    try:
        z=float(x)
        return z if math.isfinite(z) else None
    except (ValueError,TypeError,OverflowError):return None

def compare(phead,pdata,shead,sdata):
    if len(phead)!=331:raise ValueError("PDC original sample axis unexpected")
    ratio=[i for i,x in enumerate(shead)
           if re.fullmatch(r".+ Log Ratio",x,re.I)
           and not re.search(r"Unshared Log Ratio$",x,re.I)]
    if len(ratio)!=330:raise ValueError("source report has not exactly 330 Log Ratio measurements")
    nonratio=[i for i in range(len(shead)) if i not in ratio
              and not re.search(r"Unshared Log Ratio$",shead[i],re.I)]
    keyset={v[0] for v in pdata if isinstance(v[0],str) and v[0]}
    candidates=collections.defaultdict(list)
    for row in sdata:
        keys=set(v.strip() for k in nonratio[:15] for v in [row[k]]
                 if isinstance(v,str) and v.strip() in keyset)
        for key in keys:candidates[key].append(row)
    punique={}
    for row in pdata:
        if isinstance(row[0],str) and row[0]:
            punique.setdefault(row[0],[]).append(row)
    paired=[(pp[0],candidates[key][0]) for key,pp in punique.items()
            if len(pp)==1 and len(candidates[key])==1]
    selected=paired[:12]
    # Numeric matching is only a source-parity candidate and cannot establish
    # a biological sample identity without the actual original ID manifest.
    fingerprints_a=collections.defaultdict(list)
    fingerprints_b=collections.defaultdict(list)
    feature_n=len(selected)
    if feature_n>=4:
        for j in range(330):
            q=[finite(row[1+j]) for row,source in selected]
            s=[finite(source[ratio[jj]]) for row,source in selected for jj in []] # deliberately no ordinal alignment
            if all(z is not None for z in q):
                fingerprints_a[tuple(round(z,6) for z in q)].append(j)
            # A separate source report sample column has its OWN name/position;
            # comparing fingerprints across all 330 positions does not assume
            # corresponding indices.
        for jj,col in enumerate(ratio):
            v=[finite(source[col]) for row,source in selected]
            if all(z is not None for z in v):
                fingerprints_b[tuple(round(z,6) for z in v)].append(jj)
    unique=[k for k in fingerprints_a if len(fingerprints_a[k])==1
            and len(fingerprints_b.get(k,[]))==1]
    return {
       "PDC_quant_prefix_complete_feature_rows_n":len(pdata),
       "CDAP_TMT_prefix_complete_feature_rows_n":len(sdata),
       "original_CDAP_report_Log_Ratio_columns_n":len(ratio),
       "CDAP_feature_values_matching_PDC_exact_names_n":len(candidates),
       "unique_exact_feature_row_pairs_n":len(paired),
       "numeric_fingerprint_feature_pairs_used_n":feature_n if feature_n>=4 else 0,
       "PDC_columns_with_complete_numeric_fingerprint_n":sum(len(v) for v in fingerprints_a.values()),
       "CDAP_columns_with_complete_numeric_fingerprint_n":sum(len(v) for v in fingerprints_b.values()),
       "unique_exact_numeric_fingerprint_candidates_n":len(unique),
       "numerical_matches_prove_official_aliquot_identity":False,
       "source_entire_numeric_assay_parity_certified":False,
       "raw_feature_or_sample_names_and_cells_exported":False
    }

def audit():
    out={"pdc_study_id":"PDC000198","status":"ACCESS_HOLD",
         "source_identifiers_or_numeric_values_exported":False,
         "original_330_aliquot_link_proven":False,
         "patient_rna_protein_match_proven":False}
    try:
        urls,digest=retry(lambda:original_tmt_reports("PDC000198"))
        out.update({"official_source_catalog_SHA256":digest,"TMT_report_candidates_n":len(urls)})
        if len(urls)!=1:
            out["status"]="SOURCE_TMT_REPORT_AMBIGUOUS";return out
        ph,pr,pdigest=retry(lambda:pdc_prefix())
        sh,sr,sdigest=retry(lambda:report_prefix(urls[0]))
        out.update({"PDC_quant_prefix_SHA256":pdigest,"CDAP_report_prefix_SHA256":sdigest,
                    "source_prefix_parity":compare(ph,pr,sh,sr),
                    "status":"BOUNDED_SOURCE_NUMERIC_PREFIX_COMPARED"})
    except Exception as e:
        out["error_class"]=type(e).__name__
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True);p.add_argument("--review-out",required=True)
    a=p.parse_args();result={"schema":"hcc-stage42-original-CDAP-PDC-numeric-parity-pilot-v1",
             "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
             "source":audit()}
    x=result["source"]
    if any(x.get(k) is not False for k in ("source_identifiers_or_numeric_values_exported",
                                           "original_330_aliquot_link_proven","patient_rna_protein_match_proven")):
        raise ValueError("unverified scientific inference")
    reviewed={"schema":"hcc-stage42-numeric-parity-pilot-postflight-v1",
         "technical_source_prefix_compared":x["status"]=="BOUNDED_SOURCE_NUMERIC_PREFIX_COMPARED",
         "G30_B":"HOLD_NO_OFFICIAL_QUANT_TO_ALIQUOT_SOURCE_ID",
         "G30_C":"HOLD_NODE_RNA_MATRIX","G30_D":"STOP_UNVERIFIED_MATCHED_RNA_PROTEIN",
         "scientific_claims_allowed":False}
    for dst,obj in ((a.out,result),(a.review_out,reviewed)):
        f=pathlib.Path(dst);f.parent.mkdir(parents=True,exist_ok=True)
        f.write_text(json.dumps(obj,indent=2)+"\n")
    print("HCC_STAGE42_NUMERIC_SOURCE",json.dumps(result))
    print("HCC_STAGE42_REVIEW",json.dumps(reviewed))
if __name__=="__main__":main()
