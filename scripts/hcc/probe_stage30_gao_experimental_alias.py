#!/usr/bin/env python3
"""PDC experimentalMetadata documented aliquot-run aliases ↔ quant table headers.

Source PDC same-study alias matching only. All raw alias, patient, case,
aliquot, plex and signed file URL strings live exclusively in runner memory.
Report only bounded header geometry and aggregate exact alias mapping counts.
"""
import argparse,collections,datetime as dt,json,pathlib,re
from probe_stage30_gao_catalog import PDC,pdc_query,STUDY_FIELDS
from probe_stage30_gao_crossassembly_headers import sample as read_header
from probe_stage30_gao_quant_header import safe_public_https

def metadata(pid):
    q='{experimentalMetadata(pdc_study_id:"'+pid+'"){study_run_metadata{study_run_metadata_submitter_id aliquot_run_metadata{aliquot_id aliquot_submitter_id label alias}}}}'
    data,sha=pdc_query(q)
    entries=data.get("experimentalMetadata")
    if not isinstance(entries,list):raise ValueError("source metadata array missing")
    alias_to_aliquots=collections.defaultdict(set)
    assign=collections.Counter();group_n=0;alias_n=0
    for entry in entries:
        groups=entry.get("study_run_metadata") or []
        if not isinstance(groups,list):raise ValueError("nested run metadata changed")
        for g in groups:
            group_n+=1
            aliquot_runs=g.get("aliquot_run_metadata") or []
            if not isinstance(aliquot_runs,list):raise ValueError("aliquot run source changed")
            for r in aliquot_runs:
                al=str(r.get("aliquot_id") or "").strip()
                alias=str(r.get("alias") or "").strip()
                lab=str(r.get("label") or "").strip()
                if al and alias:
                    alias_to_aliquots[alias].add(al)
                    alias_n+=1
                assign[lab or "UNKNOWN_LABEL"]+=1
    stats={"source_response_sha256":sha,"experimental_metadata_top_level_records_n":len(entries),
        "study_run_metadata_records_n":group_n,
        "aliquot_run_rows_n":sum(assign.values()),
        "run_alias_assignments_with_aliquot_id_n":alias_n,
        "distinct_aliquot_run_alias_n":len(alias_to_aliquots),
        "aliquot_runs_with_missing_alias_n":sum(assign.values())-alias_n,
        "label_counts":dict(sorted(assign.items())),
        "source_alias_or_aliquot_IDs_exported":False}
    return alias_to_aliquots,stats

def files(pid):
    data,_=pdc_query('{study(pdc_study_id:"'+pid+'"){'+STUDY_FIELDS+'}}')
    st=data.get("study")
    if isinstance(st,list):
        if len(st)!=1:raise ValueError("ambiguous study")
        st=st[0]
    if not isinstance(st,dict) or st.get("pdc_study_id")!=pid:raise ValueError("study missing")
    uid=st.get("study_id")
    if not isinstance(uid,str) or not re.fullmatch(r"[0-9a-fA-F-]{32,36}",uid):
        raise ValueError("no versioned study ID")
    data,_=pdc_query('{filesPerStudy(study_id:"'+uid+'",data_category:"Protein Assembly",offset:0,limit:30){file_name file_format data_category signedUrl{url}}}')
    rr=data.get("filesPerStudy")
    if not isinstance(rr,list) or len(rr)>30 or any(x.get("data_category")!="Protein Assembly" for x in rr):
        raise ValueError("source file list absent")
    scored=[]
    for f in rr:
        name=str(f.get("file_name") or "").lower()
        fmt=str(f.get("file_format") or "").lower()
        signed=f.get("signedUrl") or []
        if isinstance(signed,dict):signed=[signed]
        url=next((x.get("url") for x in signed if isinstance(x,dict) and safe_public_https(x.get("url"))),None)
        if not url:continue
        score= (4 if fmt=="tsv" else 0) +(3 if re.search(r"protein|abundance|quant",name) else 0) - (3 if re.search(r"peptide|psm",name) else 0)
        scored.append((score,url))
    scored.sort(key=lambda z:z[0],reverse=True)
    return scored[:6]

def source(pid):
    out={"pdc_study_id":pid,"status":"ACCESS_HOLD",
         "source_alias_header_aliquot_identifiers_exported":False,
         "matched_RNA_protein_patients_verified":False,
         "complete_quant_numeric_matrix_verified":False}
    try:
        aliases,meta=metadata(pid)
        filelist=files(pid)
        results=[]
        for rank,url in filelist:
            data,keys=read_header(url)
            stat={"header_status":data.get("status"),
                  "header_token_n":len(keys) if keys is not None else 0}
            if keys is not None:
                exact=keys & set(aliases)
                unique_aliquots=set().union(*(aliases[x] for x in exact)) if exact else set()
                ambiguous_aliases=sum(len(aliases[a])!=1 for a in exact)
                stat.update({"exact_experimental_metadata_alias_token_match_n":len(exact),
                    "aliquot_UUIDs_resolved_by_exact_alias_n":len(unique_aliquots),
                    "ambiguous_alias_to_multiple_aliquots_n":ambiguous_aliases,
                    "read_only_header_sample_or_alias_exported":False})
            results.append(stat)
        out.update({"status":"OFFICIAL_EXPERIMENTAL_METADATA_ALIAS_AND_HEADER_CHECKED",
                    "experimental_metadata":meta,
                    "bounded_file_candidates_n":len(filelist),
                    "bounded_header_alias_match":results,
                    "source_quant_header_is_definitively_gene_aliquot_matrix":False})
    except Exception as e:
        out["error_type"]=type(e).__name__
    return out

def review(obj):
    if (obj.get("schema")!="hcc-stage30-experimental-alias-quant-header-v1" or
        {r.get("pdc_study_id") for r in obj.get("studies",[])}!=set(PDC)):
        raise ValueError("source schema")
    for r in obj["studies"]:
        if r.get("status") not in ("ACCESS_HOLD","OFFICIAL_EXPERIMENTAL_METADATA_ALIAS_AND_HEADER_CHECKED"):
            raise ValueError("source alias status unsupported")
        for key in ("source_alias_header_aliquot_identifiers_exported",
                    "matched_RNA_protein_patients_verified",
                    "complete_quant_numeric_matrix_verified",
                    "source_quant_header_is_definitively_gene_aliquot_matrix"):
            if r.get(key) is True:raise ValueError("source identifier leak or unsupported matrix")
        if r["status"]=="OFFICIAL_EXPERIMENTAL_METADATA_ALIAS_AND_HEADER_CHECKED":
            if not 0<=r["bounded_file_candidates_n"]<=6:raise ValueError("source header candidates > limit")
            m=r.get("experimental_metadata")
            if not isinstance(m,dict) or m.get("source_alias_or_aliquot_IDs_exported") is not False:
                raise ValueError("source alias leak")
            for x in r.get("bounded_header_alias_match",[]):
                if x.get("read_only_header_sample_or_alias_exported") is True:
                    raise ValueError("sample header leak")
                for key in ("exact_experimental_metadata_alias_token_match_n",
                            "aliquot_UUIDs_resolved_by_exact_alias_n"):
                    if key in x and (type(x[key]) is not int or not 0<=x[key]<=5000):
                        raise ValueError("implausible exact source alias count")
    return {"schema":"hcc-stage30-alias-header-postflight-v1",
            "studies_with_experimental_alias_source_records_n":
                sum(r["status"]=="OFFICIAL_EXPERIMENTAL_METADATA_ALIAS_AND_HEADER_CHECKED" for r in obj["studies"]),
            "max_exact_alias_match_tokens_by_study":{
                r["pdc_study_id"]:max((x.get("exact_experimental_metadata_alias_token_match_n",0)
                   for x in r.get("bounded_header_alias_match",[])),default=0) for r in obj["studies"]},
            "G30_B":"PARTIAL_SOURCE_METADATA_PROBE_ONLY_HOLD_FULL_NUMERIC_QUANTITY",
            "G30_C":"HOLD_NODE_RNA_MATRIX",
            "G30_D":"STOP_NO_SAME_PATIENT_RNA_PROTEIN_JOIN",
            "G30_E":"NOT_STARTED",
            "biological_replication_verified":False}

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True);p.add_argument("--review-out",required=True)
    args=p.parse_args()
    v={"schema":"hcc-stage30-experimental-alias-quant-header-v1",
       "checked_at_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
       "studies":[source(pid) for pid in PDC]}
    path=pathlib.Path(args.out);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(v,indent=2)+"\n")
    r=review(v)
    path=pathlib.Path(args.review_out);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(r,indent=2)+"\n")
    print(json.dumps(r,indent=2))
    print("EXPERIMENTAL_ALIAS_SOURCE_COUNTS",json.dumps(v["studies"]))
if __name__=="__main__":main()
