#!/usr/bin/env python3
"""Source PDC CDAP *.sample.txt internal assay header↔TMT/aliquot bridge probe.

Public metadata + bounded 96 KiB first bytes of maximum three specifically
named sample/experimental-design text files per PDC study. Compare private
tokens to PDC quantDataMatrix header and official aliquot ID namespaces.
Save only aggregate counts/source hashes, not original row records, aliquot
IDs, file names, temporary signed URLs, gene or quantitative values.

Missing companion file in Protein Assembly LIST is ACCESS_HOLD, not evidence
that no file exists in the wider repository/author supplements.
"""
import argparse,datetime as dt,hashlib,json,pathlib,re,urllib.request,collections
from probe_stage30_gao_catalog import PDC,pdc_query,STUDY_FIELDS
from probe_stage30_gao_quant_header import safe_public_https
from probe_stage30_gao_quant_api_header_to_biospecimen import get_quant_header,official_biospecimen

LIMIT=98304

def file_records(pid):
    d,_=pdc_query('{study(pdc_study_id:"'+pid+'"){'+STUDY_FIELDS+'}}')
    s=d.get("study")
    if isinstance(s,list):
        if len(s)!=1:raise ValueError("ambiguous source study")
        s=s[0]
    if not isinstance(s,dict) or s.get("pdc_study_id")!=pid:raise ValueError("source PDC study absent")
    uid=s.get("study_id")
    if not isinstance(uid,str) or not re.fullmatch(r"[0-9a-fA-F-]{32,36}",uid):
        raise ValueError("study version not pinned")
    q='{filesPerStudy(study_id:"'+uid+'",data_category:"Protein Assembly",offset:0,limit:30){file_name file_format data_category signedUrl{url}}}'
    d,digest=pdc_query(q)
    rows=d.get("filesPerStudy")
    if not isinstance(rows,list) or len(rows)>30:raise ValueError("assembly category source shape drift")
    eligible=[]
    for r in rows:
        if r.get("data_category")!="Protein Assembly":raise ValueError("wrong data_category")
        name=str(r.get("file_name") or "")
        if not (re.search(r"\.sample\.txt(?:\.gz)?$",name,re.I) or
                re.search(r"(experimental.?design|sample.?map|sample.?metadata).*\.t(?:sv|xt)$",name,re.I)):
            continue
        urls=r.get("signedUrl") or []
        if isinstance(urls,dict):urls=[urls]
        url=next((x.get("url") for x in urls if isinstance(x,dict) and safe_public_https(x.get("url"))),None)
        if url:eligible.append(url)
    return len(rows),eligible[:3],digest

def sample_table_preview(url):
    if not safe_public_https(url):raise ValueError("untrusted file host")
    req=urllib.request.Request(url,headers={
      "Range":"bytes=0-98303","User-Agent":"HCC-Stage30-CDAP-Design-Bridge/1"})
    with urllib.request.urlopen(req,timeout=25) as rsp:
        b=rsp.read(LIMIT+1)
    b=b[:LIMIT]
    if b[:2]==b"\x1f\x8b":raise ValueError("compressed companion cannot be parsed as TSV")
    txt=b.decode("utf-8","replace")
    lines=[x for x in txt.splitlines() if x.strip() and not x.strip().startswith("#")]
    if not lines:raise ValueError("empty source companion")
    sep="\t" if lines[0].count("\t")>=lines[0].count(",") else ","
    table=[x.split(sep) for x in lines[:1500]]
    return {"row_count_previewed_n":len(table)-1,
        "first_line_column_count_n":len(table[0]),
        "source_prefix_sha256":hashlib.sha256(b).hexdigest(),
        "read_only_bytes_n":len(b),
        "possible_companion_partial_tail":len(b)>=LIMIT,
        "source_names_or_rows_exported":False},table

def compare_rows(table,quant_header,aliases):
    # Any one raw source table row that has both a quant column header and
    # a known aliquot label is a candidate bridge. Not guaranteed unique or
    # directional without exact column schema + author run/channel evidence.
    head=set(quant_header)
    aliq=set(aliases["aliquot_id"])|set(aliases["aliquot_submitter_id"])
    quant_tokens_seen=set()
    aliq_tokens_seen=set()
    both=0
    for row in table[1:]:
        vals=set(x.strip().strip('"') for x in row if x.strip())
        q=vals&head
        a=vals&aliq
        quant_tokens_seen|=q
        aliq_tokens_seen|=a
        if q and a:both+=1
    return {"companion_quant_header_exact_common_tokens_n":len(quant_tokens_seen),
        "companion_official_aliquot_UUID_or_submitter_exact_common_tokens_n":len(aliq_tokens_seen),
        "companion_rows_containing_both_quant_alias_and_aliquot_alias_n":both,
        "mapping_complete_and_unambiguous_verified":False,
        "source_biospecimen_or_quant_alias_exported":False}

def one(pid):
    out={"pdc_study_id":pid,"status":"ACCESS_HOLD",
         "source_file_urls_rows_and_aliquot_names_exported":False,
         "validated_full_assay_col_to_patient_join":False,
         "RNA_patient_join_verified":False}
    try:
        n,files,digest=file_records(pid)
        out.update({"assembly_category_file_n":n,
           "assembly_metadata_response_sha256":digest,
           "sample_design_companion_candidate_n":len(files)})
        if not files:
            out["status"]="NO_SAMPLE_DESIGN_COMPANION_IN_ASSEMBLY_CATEGORY"
            return out
        _,qh=get_quant_header(pid)
        aliases=official_biospecimen(pid)
        candidates=[]
        for u in files:
            try:
                info,rows=sample_table_preview(u)
                info.update(compare_rows(rows,qh,aliases))
                info["status"]="COMPANION_BOUNDED_TABLE_PROBED"
            except Exception as e:
                info={"status":"COMPANION_ACCESS_HOLD",
                      "failure_class":type(e).__name__,
                      "source_names_or_rows_exported":False}
            candidates.append(info)
        out["status"]="SAMPLE_DESIGN_COMPANION_SOURCE_PROBED"
        out["bounded_companion_preview"]=candidates
    except Exception as e:
        out["error_class"]=type(e).__name__
    return out

def verify(data):
    if data.get("schema")!="hcc-stage30-PDC-sample-design-companion-v1":
        raise ValueError("invalid source schema")
    rows=data.get("studies")
    if not isinstance(rows,list) or {x.get("pdc_study_id") for x in rows}!=set(PDC):
        raise ValueError("source identity drift")
    for x in rows:
        if x.get("status") not in ("ACCESS_HOLD","NO_SAMPLE_DESIGN_COMPANION_IN_ASSEMBLY_CATEGORY",
                                 "SAMPLE_DESIGN_COMPANION_SOURCE_PROBED"):
            raise ValueError("invalid companion status")
        for k in ("source_file_urls_rows_and_aliquot_names_exported",
                  "validated_full_assay_col_to_patient_join","RNA_patient_join_verified"):
            if x.get(k) is not False:raise ValueError("unsafe claim or ID exposure")
    return {"schema":"hcc-stage30-sample-design-companion-postflight-v1",
        "pdc_studies_with_bounded_sample_design_companion_n":
          sum(x["status"]=="SAMPLE_DESIGN_COMPANION_SOURCE_PROBED" for x in rows),
        "G30_B":"HOLD_NO_VALIDATED_ASSAY_COLUMN_TO_ALIQUOT_JOIN",
        "G30_C":"HOLD_NODE_RNA_MATRIX",
        "G30_D":"STOP_NO_MATCHED_RNA_PROTEOME_PATIENTS",
        "independent_biological_validation_performed":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True);p.add_argument("--review-out",required=True)
    a=p.parse_args()
    d={"schema":"hcc-stage30-PDC-sample-design-companion-v1",
       "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
       "studies":[one(pid) for pid in PDC]}
    target=pathlib.Path(a.out);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(d,indent=2)+"\n")
    r=verify(d);target2=pathlib.Path(a.review_out);target2.parent.mkdir(parents=True,exist_ok=True)
    target2.write_text(json.dumps(r,indent=2)+"\n")
    print(json.dumps(r,indent=2))
    print("PDC_CDAP_SAMPLE_DESIGN_COMPANION",json.dumps(d["studies"]))
if __name__=="__main__":main()
