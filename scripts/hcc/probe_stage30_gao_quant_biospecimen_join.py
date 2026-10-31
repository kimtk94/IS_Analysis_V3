#!/usr/bin/env python3
"""PDC 198/199 official quant table prefix header ↔ biopsies per-study ID check.

Bounded 16K HTTP Range per at most 6 Protein Assembly files/study; original
PDC case/sample/aliquot identifiers AND signed URLs remain in runner memory.
Aggregate exact-match counts only; this is NOT a full numeric matrix or RNA join.
"""
import argparse,datetime as dt,json,pathlib,re,urllib.request
from probe_stage30_gao_catalog import PDC,pdc_query,STUDY_FIELDS
from probe_stage30_gao_quant_header import safe_public_https

def sample_header(url):
    if not safe_public_https(url):return None,"HOST_HOLD"
    try:
        req=urllib.request.Request(url,headers={"Range":"bytes=0-16383",
            "User-Agent":"HCC-Stage30-Sample-ID-QC-ReadOnly/1"})
        with urllib.request.urlopen(req,timeout=25) as rsp:
            raw=rsp.read(16384)
        if raw[:2]==b"\x1f\x8b":return None,"GZIP"
        txt=raw.decode("utf-8","replace").splitlines()
        if not txt:return None,"EMPTY"
        first=txt[0]
        delimiter="\t" if first.count("\t")>=first.count(",") else ","
        toks=[s.strip().strip('"') for s in first.split(delimiter)]
        if not 2<len(toks)<1500:return None,"NONTABULAR"
        return set(toks[1:]),"HEADER_READ"
    except Exception:
        return None,"ACCESS_HOLD"

def get_biospecimen(pid):
    d,_=pdc_query('{biospecimenPerStudy(pdc_study_id:"'+pid+'"){case_id case_submitter_id sample_id sample_submitter_id aliquot_id aliquot_submitter_id sample_type}}')
    rows=d.get("biospecimenPerStudy")
    if not isinstance(rows,list) or not 1<=len(rows)<=10000:
        raise ValueError("source study biospecimen unavailable")
    out={}
    for r in rows:
        aliquot=(r.get("aliquot_submitter_id") or r.get("aliquot_id") or "").strip()
        samp=(r.get("sample_submitter_id") or r.get("sample_id") or "").strip()
        case=(r.get("case_submitter_id") or r.get("case_id") or "").strip()
        if not aliquot or not samp or not case:
            raise ValueError("missing source biospecimen ID")
        typ=str(r.get("sample_type") or "").lower()
        typ="Normal" if "normal" in typ or "non-tumor" in typ else (
            "Tumor" if "tumor" in typ or "malignant" in typ or "primary" in typ else "Other")
        out[aliquot]={"sample":samp,"case":case,"kind":typ,
            "aliquot_submitter_id":str(r.get("aliquot_submitter_id") or "").strip(),
            "aliquot_id":str(r.get("aliquot_id") or "").strip(),
            "sample_submitter_id":str(r.get("sample_submitter_id") or "").strip(),
            "sample_id":str(r.get("sample_id") or "").strip(),
            "case_submitter_id":str(r.get("case_submitter_id") or "").strip(),
            "case_id":str(r.get("case_id") or "").strip()}
    return out

def study(pid):
    result={"pdc_study_id":pid,"status":"SOURCE_ACCESS_HOLD",
        "gene_aliquot_full_matrix_verified":False,
        "signed_urls_and_raw_identifiers_exported":False,
        "RNA_case_identity_verified":False}
    try:
        index=get_biospecimen(pid)
        d,_=pdc_query('{study(pdc_study_id:"'+pid+'"){'+STUDY_FIELDS+'}}')
        st=d.get("study")
        if isinstance(st,list):
            if len(st)!=1:raise ValueError("ambiguous study")
            st=st[0]
        if not isinstance(st,dict) or st.get("pdc_study_id")!=pid:
            raise ValueError("source study error")
        uid=st.get("study_id")
        if not isinstance(uid,str) or not re.fullmatch(r"[0-9a-fA-F-]{32,36}",uid):
            raise ValueError("uuid unavailable")
        d,_=pdc_query('{filesPerStudy(study_id:"'+uid+'",data_category:"Protein Assembly",offset:0,limit:30){file_name data_category file_format signedUrl {url}}}')
        files=d.get("filesPerStudy")
        if not isinstance(files,list) or any(x.get("data_category")!="Protein Assembly" for x in files):
            raise ValueError("no assembly files")
        scored=[]
        for file in files:
            fmt=str(file.get("file_format") or "").lower()
            name=str(file.get("file_name") or "").lower()
            urls=file.get("signedUrl") or []
            if isinstance(urls,dict):urls=[urls]
            url=next((x.get("url") for x in urls if isinstance(x,dict) and safe_public_https(x.get("url"))),None)
            if not url:continue
            score=(4 if fmt=="tsv" else 0)+(
                 3 if re.search(r"protein|proteome|quant|abundan",name) else 0)-(
                 3 if re.search(r"peptid|psm|spectrum",name) else 0)
            scored.append((score,url))
        scored.sort(key=lambda t:t[0],reverse=True)
        previews=[]
        # Not just the first 40 raw files: six targeted assembly files.
        for _,url in scored[:6]:
            header,status=sample_header(url)
            metrics={"status":status}
            if header is not None:
                matches=[r for a,r in index.items() if a in header]
                matched_sample_ids={r["sample"] for r in matches}
                matched_case_ids={r["case"] for r in matches}
                # PDC aliquot submitter IDs may not be the file's column
                # namespace. Independently test all six documented PDC
                # case/sample/aliquot aliases. No fuzzy/cross-patient joining.
                namespaces={}
                for key in ("aliquot_submitter_id","aliquot_id",
                            "sample_submitter_id","sample_id",
                            "case_submitter_id","case_id"):
                    names={r.get(key) for r in index.values() if r.get(key)}
                    namespaces[key]=len(names & header)
                # Conservative case-insensitive and trailing-whitespace-only
                # comparison to detect trivial quoting/case changes.
                hupper={k.strip().upper() for k in header}
                namespaces_ci={}
                for key in namespaces:
                    names={r.get(key).upper() for r in index.values() if r.get(key)}
                    namespaces_ci[key]=len(names & hupper)
                metrics.update({
                    "nonfirst_header_unique_token_n":len(header),
                    "exact_PDC_case_sample_aliquot_alias_namespace_match_n":namespaces,
                    "casefold_PDC_alias_namespace_match_n":namespaces_ci,
                    "exact_matched_PDC_aliquot_submitter_IDs_n":len(matches),
                    "exact_matched_PDC_unique_samples_n":len(matched_sample_ids),
                    "exact_matched_PDC_unique_cases_n":len(matched_case_ids),
                    "exact_matched_Tumor_aliquot_n":sum(r["kind"]=="Tumor" for r in matches),
                    "exact_matched_Normal_aliquot_n":sum(r["kind"]=="Normal" for r in matches),
                    "exact_matched_Other_aliquot_n":sum(r["kind"]=="Other" for r in matches),
                    "sample_case_or_aliquot_IDs_exported":False})
            previews.append(metrics)
        result.update({
            "status":"SOURCE_ASSEMBLY_HEADER_VS_PDC_ALIQUOT_METADATA_CHECKED",
            "PDC_source_unique_aliquot_n":len(index),
            "PDC_assembly_file_metadata_n":len(files),
            "num_assembly_headers_probed_n":len(previews),
            "per_candidate_aliquot_ID_overlap_counts":previews,
            "full_matrix_numeric_values_and_gene_annotations_verified":False})
    except Exception as e:
        result["failure_class"]=type(e).__name__
    return result

def review(v):
    if (v.get("schema")!="hcc-stage30-quant-header-vs-pdc-aliquot-source-v1" or
        {s.get("pdc_study_id") for s in v.get("studies",[])}!=set(PDC)):
        raise ValueError("PDC source target mismatch")
    for s in v["studies"]:
        if s.get("status") not in ("SOURCE_ACCESS_HOLD","SOURCE_ASSEMBLY_HEADER_VS_PDC_ALIQUOT_METADATA_CHECKED"):
            raise ValueError("invalid cohort access")
        for flag in ("gene_aliquot_full_matrix_verified","signed_urls_and_raw_identifiers_exported",
                     "RNA_case_identity_verified","full_matrix_numeric_values_and_gene_annotations_verified"):
            if s.get(flag) is True:raise ValueError("unsupported file or RNA match claim")
        if s.get("status")=="SOURCE_ASSEMBLY_HEADER_VS_PDC_ALIQUOT_METADATA_CHECKED":
            if s.get("PDC_source_unique_aliquot_n")!=331 or not 0<=s.get("num_assembly_headers_probed_n",-1)<=6:
                raise ValueError("changed original study specimen count or candidate cap")
            for m in s["per_candidate_aliquot_ID_overlap_counts"]:
                for k in ("exact_matched_PDC_aliquot_submitter_IDs_n",
                          "exact_matched_PDC_unique_samples_n",
                          "exact_matched_PDC_unique_cases_n"):
                    n=m.get(k,0)
                    if type(n) is not int or not 0<=n<=331:raise ValueError("source count impossible")
                if m.get("sample_case_or_aliquot_IDs_exported") is True:raise ValueError("patient ID leak")
    return {"schema":"hcc-stage30-header-vs-biospecimen-postflight-v1",
        "max_matched_proteome_PDC_aliquot_headers_n":max((m.get("exact_matched_PDC_aliquot_submitter_IDs_n",0) for m in v["studies"][0].get("per_candidate_aliquot_ID_overlap_counts",[])),default=0),
        "max_matched_phospho_PDC_aliquot_headers_n":max((m.get("exact_matched_PDC_aliquot_submitter_IDs_n",0) for m in v["studies"][1].get("per_candidate_aliquot_ID_overlap_counts",[])),default=0),
        "G30_B":"PARTIAL_HEADER_JOIN_EVIDENCE_BUT_FULL_NUMERICAL_MATRIX_HOLD",
        "G30_C":"HOLD_NODE_RNA_DATA",
        "G30_D":"STOP_NO_INDEPENDENT_MATCHED_RNA_AND_PROTEIN",
        "G30_E":"NOT_STARTED",
        "external_cohort_gene_effects_reproduced":False}

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True);p.add_argument("--review-out",required=True)
    a=p.parse_args()
    x={"schema":"hcc-stage30-quant-header-vs-pdc-aliquot-source-v1",
       "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
       "max_signed_assembly_urls_read_per_study":6,
       "max_source_bytes_per_header":16384,
       "studies":[study(pid) for pid in PDC]}
    out=pathlib.Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(x,indent=2)+"\n")
    y=review(x)
    out=pathlib.Path(a.review_out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(y,indent=2)+"\n")
    print(json.dumps(y,indent=2))
    print("PDC_QUANT_ALIQUOT_LINK_COUNTS",json.dumps(x["studies"]))
if __name__=="__main__":main()
