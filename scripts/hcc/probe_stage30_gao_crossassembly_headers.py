#!/usr/bin/env python3
"""Bounded comparison of quantitative header columns across two Gao PDC studies.

Only up to four official Protein Assembly files per PDC study; each HTTPS Range
reads <=16KiB. No patient/header gene/aliquot identifier, signed URL,
filename, molecular value or raw response exported. Header hashes and aggregate
set intersections are proxies for cohort comparability, NOT verified joins.
"""
import argparse,datetime as dt,hashlib,json,pathlib,re,urllib.request
from probe_stage30_gao_catalog import PDC,STUDY_FIELDS,pdc_query
from probe_stage30_gao_quant_header import safe_public_https
BUDGET=16384

def sample(u):
    if not safe_public_https(u):return {"status":"URL_HOST_HOLD"},None
    try:
        req=urllib.request.Request(u,headers={"Range":"bytes=0-16383",
            "User-Agent":"HCC-Stage30-Source-Header-Intersection/1"})
        with urllib.request.urlopen(req,timeout=23) as resp:
            raw=resp.read(BUDGET)
            code=resp.status
        if raw[:2]==b"\x1f\x8b":return {"status":"COMPRESSED_HEADER_HOLD","http_status":code},None
        content=raw.decode("utf-8","replace")
        lines=content.splitlines()
        if len(lines)<2:return {"status":"INSUFFICIENT_FIRST_LINES","http_status":code},None
        header=lines[0]
        # Deliberately exclude the non-sample identifier field (first column).
        sep="\t" if header.count("\t")>=header.count(",") else ","
        fields=header.split(sep)
        if len(fields)<3 or len(fields)>1500:
            return {"status":"HEADER_NOT_TABULAR","http_status":code},None
        keys=[x.strip().strip('"') for x in fields[1:] if x.strip()]
        if len(keys)!=len(set(keys)):
            return {"status":"NONUNIQUE_COLUMN_NAMES","http_status":code},None
        nums=0;total=0
        complete=0
        for line in lines[1:6]:
            vals=line.split(sep)
            if len(vals)!=len(fields):continue
            complete+=1
            for cell in vals[1:]:
                total+=1
                try:
                    z=float(cell.strip().strip('"'))
                    if z==z and abs(z)<1e100:nums+=1
                except (ValueError,OverflowError):pass
        return {"status":"BOUNDED_TABULAR_HEADER",
                "header_sha256":hashlib.sha256(header.encode()).hexdigest(),
                "column_n":len(fields),
                "unique_nonfirst_header_token_n":len(keys),
                "complete_numeric_data_rows_previewed_n":complete,
                "numeric_cells_count":nums,
                "numeric_cell_ratio":round(nums/total,6) if total else 0,
                "http_status":code,
                "source_sample_column_ids_exported":False,
                "sample_to_case_mapping_verified":False},set(keys)
    except Exception as exc:
        return {"status":"SOURCE_FILE_ACCESS_HOLD","error_class":type(exc).__name__},None

def study(pid):
    out={"pdc_study_id":pid,"status":"ACCESS_HOLD",
         "source_file_or_patient_identifiers_exported":False}
    internal=[]
    try:
        d,_=pdc_query('{study(pdc_study_id:"'+pid+'"){'+STUDY_FIELDS+'}}')
        a=d.get("study")
        if isinstance(a,list):
            if len(a)!=1:raise ValueError("ambiguous study")
            a=a[0]
        if not isinstance(a,dict) or a.get("pdc_study_id")!=pid:
            raise ValueError("study unavailable")
        uid=a.get("study_id")
        if not isinstance(uid,str) or not re.fullmatch(r"[0-9a-fA-F-]{32,36}",uid):
            raise ValueError("no study UUID")
        data,_=pdc_query('{filesPerStudy(study_id:"'+uid+'",data_category:"Protein Assembly",offset:0,limit:30){file_name file_format data_category signedUrl {url}}}')
        rows=data.get("filesPerStudy")
        if not isinstance(rows,list) or len(rows)>30 or any(x.get("data_category")!="Protein Assembly" for x in rows):
            raise ValueError("wrong assembly source")
        options=[]
        for r in rows:
            name=str(r.get("file_name") or "").lower()
            fmt=str(r.get("file_format") or "").lower()
            urls=r.get("signedUrl") or []
            if isinstance(urls,dict):urls=[urls]
            u=next((x.get("url") for x in urls
                if isinstance(x,dict) and safe_public_https(x.get("url"))),None)
            if not u:continue
            # Prefer tabular protein quantitation reports, not peptide PSMs.
            rank=(5 if fmt=="tsv" else 0)+(
                4 if re.search(r"protein|proteome|abundan|quant",name) else 0)-(
                3 if re.search(r"peptid|psm|fasta|index",name) else 0)
            options.append((rank,u))
        options.sort(key=lambda x:x[0],reverse=True)
        for _,url in options[:4]:
            metric,labels=sample(url)
            internal.append(labels)
            out.setdefault("header_previews",[]).append(metric)
        out.update({"status":"CANDIDATE_HEADER_AUDIT_EXECUTED",
            "assembly_file_metadata_n":len(rows),
            "allowlisted_signed_url_candidates_n":len(options),
            "actual_header_preview_count_n":len(internal)})
    except Exception as exc:
        out["failure_class"]=type(exc).__name__
    return out,internal

def compare(one,two):
    pairs=[];max_match=0;max_jaccard=0.0
    for i,a in enumerate(one):
        if a is None:continue
        for j,b in enumerate(two):
            if b is None:continue
            n=len(a&b);den=len(a|b)
            ratio=n/den if den else 0
            pairs.append({"proteome_local_index":i,"phosphoproteome_local_index":j,
                "common_header_token_n":n,
                "union_header_token_n":den,
                "header_token_jaccard":round(ratio,6)})
            max_match=max(max_match,n);max_jaccard=max(max_jaccard,ratio)
    return {"cross_study_header_sets_compared_n":len(pairs),
            "max_common_header_token_n":max_match,
            "max_header_token_jaccard":round(max_jaccard,6),
            "comparisons":pairs,
            "this_is_not_a_verified_same_patient_join":True,
            "raw_column_label_or_patient_ids_exported":False}

def review(obj):
    if obj.get("schema")!="hcc-stage30-crossassembly-header-overlap-v1":
        raise ValueError("schema drift")
    if {s.get("pdc_study_id") for s in obj.get("studies",[])}!=set(PDC):
        raise ValueError("PDC source drift")
    for s in obj["studies"]:
        if s.get("source_file_or_patient_identifiers_exported") is not False:
            raise ValueError("data leak")
        if s.get("status") not in ("ACCESS_HOLD","CANDIDATE_HEADER_AUDIT_EXECUTED"):
            raise ValueError("unexpected response")
        if s.get("status")=="CANDIDATE_HEADER_AUDIT_EXECUTED" and not 0<=s.get("actual_header_preview_count_n",-1)<=4:
            raise ValueError("preview cap exceeded")
        for p in s.get("header_previews",[]):
            if p.get("source_sample_column_ids_exported") is True or p.get("sample_to_case_mapping_verified") is True:
                raise ValueError("patient ID leak/inference")
    x=obj.get("header_set_overlap")
    if (not isinstance(x,dict) or x.get("this_is_not_a_verified_same_patient_join") is not True or
        x.get("raw_column_label_or_patient_ids_exported") is not False):
        raise ValueError("patient join overstated")
    return {"schema":"hcc-stage30-crossassembly-postflight-v1",
            "candidate_header_overlap_computed":x["cross_study_header_sets_compared_n"]>0,
            "max_header_token_jaccard":x["max_header_token_jaccard"],
            "G30_B":"HOLD_FULL_QUANTITATIVE_GENE_ALIQUOT_MATRIX",
            "G30_C":"HOLD_RNA_NODE_FILE",
            "G30_D":"STOP_NO_SAME_PATIENT_TUMOR_TISSUE_JOIN",
            "G30_E":"NOT_STARTED",
            "biological_replication_verified":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True);p.add_argument("--review-out",required=True)
    a=p.parse_args()
    result=[study(x) for x in PDC]
    cross=compare(result[0][1],result[1][1])
    report={"schema":"hcc-stage30-crossassembly-header-overlap-v1",
        "checked_at_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
        "source_bytes_max_per_candidate":BUDGET,
        "max_candidate_files_per_study":4,
        "studies":[x[0] for x in result],
        "header_set_overlap":cross}
    out=pathlib.Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2)+"\n")
    v=review(report)
    path=pathlib.Path(a.review_out);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(v,indent=2)+"\n")
    print(json.dumps(v,indent=2))
    print("HEADER_CANDIDATES",json.dumps({s["pdc_study_id"]:{
      "status":s["status"],"headers":s.get("header_previews",[]),
      "failure":s.get("failure_class")} for s in report["studies"]}))
    print("HEADER_INTERSECTION_SUMMARY",json.dumps({
      "n_compared":cross["cross_study_header_sets_compared_n"],
      "max_common_tokens":cross["max_common_header_token_n"],
      "max_jaccard":cross["max_header_token_jaccard"]}))
if __name__=="__main__":main()
