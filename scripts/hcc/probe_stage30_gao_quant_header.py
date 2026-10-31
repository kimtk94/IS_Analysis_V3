#!/usr/bin/env python3
"""Bounded non-persistent official PDC Protein Assembly table-header inspection.

Never store URL, file name/ID, header tokens, genes, samples, patient or
aliquot names, or molecular values. At most 2 metadata-selected candidates per
study; GET Range 0..8191 only, count delimiters/numeric tokens in memory.
Do NOT infer matched quantitative HCC cohort from a table-like preview.
"""
import argparse,datetime as dt,hashlib,json,pathlib,re,urllib.parse,urllib.request
from probe_stage30_gao_catalog import PDC,pdc_query,STUDY_FIELDS

MAX_HEADER=8192
ALLOWED_HOSTS=("pdc.cancer.gov","pdc-docs.cancer.gov",
               "amazonaws.com","cloudfront.net","storage.googleapis.com",
               "cancer.gov")

def safe_public_https(url):
    if not isinstance(url,str) or len(url)>6000:return False
    try:
        p=urllib.parse.urlsplit(url)
        host=(p.hostname or "").lower()
        return (p.scheme=="https" and bool(host) and
                not p.username and not p.password and
                all(x not in host for x in ("localhost","127.0.0.1","169.254.")) and
                any(host==d or host.endswith("."+d) for d in ALLOWED_HOSTS))
    except ValueError:return False

def inspect_bytes(body):
    if body.startswith(b"\x1f\x8b"):return {"status":"GZIP_HEADER_ONLY","bytes":len(body)}
    text=body.decode("utf-8","replace")
    lines=[s for s in text.splitlines() if s and not s.startswith("#")][:6]
    if len(lines)<2:return {"status":"NOT_ENOUGH_TABLE_LINES","bytes":len(body)}
    candidates=("\t",",",";")
    sep=max(candidates,key=lambda x:lines[0].count(x))
    cols=lines[0].count(sep)+1
    counts=[len(x.split(sep)) for x in lines]
    rectangular=sum(n==cols for n in counts)
    toks=[]
    for row in lines[1:]:
        toks.extend(row.split(sep)[1:])
    numeric=0
    for x in toks:
        try:
            z=float(x.strip().strip('"'))
            if z==z and abs(z)<1e100:numeric+=1
        except (TypeError,ValueError,OverflowError):pass
    return {"status":"TEXT_TABLE_PREVIEWED",
        "source_header_line_sha256":hashlib.sha256(lines[0].encode()).hexdigest(),
        "preview_line_count":len(lines),
        "approx_header_columns":cols,
        "header_and_data_lines_rectangular_n":rectangular,
        "source_preview_cell_n":len(toks),
        "finite_numeric_cells_observed_n":numeric,
        "numeric_fraction_in_first_few_rows":round(numeric/len(toks),4) if toks else 0,
        "delimiter_kind":{"\t":"tab",",":"comma",";":"semicolon"}[sep],
        "header_or_sample_identifiers_exported":False,
        "sample_by_gene_matrix_verified":False}

def browse_one(url):
    if not safe_public_https(url):return {"status":"DOWNLOAD_HOST_NOT_ALLOWLISTED"}
    try:
        req=urllib.request.Request(url,headers={
            "User-Agent":"HCC-Stage30-Bounded-Protein-Assembly-Header/1.0",
            "Range":"bytes=0-8191"})
        with urllib.request.urlopen(req,timeout=23) as res:
            raw=res.read(MAX_HEADER)
            code=res.status
        x=inspect_bytes(raw)
        x["http_status"]=code
        x["bytes_read"]=len(raw)
        return x
    except Exception as err:
        # Never expose signed URL or urllib error string.
        return {"status":"FILE_HEADER_ACCESS_HOLD","error_class":type(err).__name__}

def study(pid):
    result={"pdc_study_id":pid,"status":"METADATA_ACCESS_HOLD",
        "signed_urls_or_filename_exported":False,
        "patient_gene_or_aliquot_identifiers_exported":False,
        "full_file_downloaded":False,
        "matched_RNA_protein_cohort_verified":False}
    try:
        d,_=pdc_query('{study(pdc_study_id:"'+pid+'"){'+STUDY_FIELDS+'}}')
        st=d.get("study")
        if isinstance(st,list):
            if len(st)!=1:raise ValueError("ambiguous study")
            st=st[0]
        if not isinstance(st,dict) or st.get("pdc_study_id")!=pid:
            raise ValueError("study not found")
        uid=st.get("study_id")
        if not isinstance(uid,str) or not re.fullmatch(r"[0-9a-fA-F-]{32,36}",uid):
            raise ValueError("missing UUID")
        d,_=pdc_query('{filesPerStudy(study_id:"'+uid+'",data_category:"Protein Assembly",offset:0,limit:30){file_name data_category file_type file_format signedUrl {url}}}')
        rows=d.get("filesPerStudy")
        if not isinstance(rows,list) or len(rows)>30 or any(x.get("data_category")!="Protein Assembly" for x in rows):
            raise ValueError("filtered source unexpected")
        scored=[]
        for row in rows:
            name=str(row.get("file_name") or "")
            typ=str(row.get("file_format") or "").lower()
            # Prioritize text/tsv protein assembly summaries, not PSM raw.
            score=(3 if typ in ("tsv","txt") else 0)
            score+=4 if re.search(r"protein|proteome|quant|abundan",name,re.I) else 0
            score-=4 if re.search(r"peptide|peptid|psm|mzml",name,re.I) else 0
            urls=row.get("signedUrl") or []
            if isinstance(urls,dict):urls=[urls]
            u=next((a.get("url") for a in urls if isinstance(a,dict) and safe_public_https(a.get("url"))),None)
            if u:scored.append((score,u))
        scored.sort(key=lambda x:x[0],reverse=True)
        chosen=scored[:2]
        previews=[browse_one(url) for _,url in chosen]
        result.update({"status":"READ_ONLY_ASSEMBLY_PREVIEW_COMPLETED",
            "assembly_file_metadata_n":len(rows),
            "signed_urls_allowed_host_present_n":len(scored),
            "candidate_header_tests_n":len(previews),
            "candidate_header_findings":previews,
            "has_tabular_numeric_header_evidence":any(x.get("status")=="TEXT_TABLE_PREVIEWED" and
                x.get("approx_header_columns",0)>2 and x.get("finite_numeric_cells_observed_n",0)>0 for x in previews),
            "gene_by_aliquot_matrix_and_tumor_normal_verified":False})
    except Exception as e:
        result["failure_category"]=type(e).__name__
    return result

def review(v):
    if v.get("schema")!="hcc-stage30-pdc-protein-assembly-bounded-header-v1" or
      {x.get("pdc_study_id") for x in v.get("studies",[])}!=set(PDC):
        raise ValueError("source contract")
    for x in v["studies"]:
        for f in ("signed_urls_or_filename_exported","patient_gene_or_aliquot_identifiers_exported",
                  "full_file_downloaded","matched_RNA_protein_cohort_verified",
                  "gene_by_aliquot_matrix_and_tumor_normal_verified"):
            if x.get(f) is True:raise ValueError("scientific or source secret leakage")
        if x.get("status") not in ("METADATA_ACCESS_HOLD","READ_ONLY_ASSEMBLY_PREVIEW_COMPLETED"):
            raise ValueError("invalid access code")
        if x.get("status")=="READ_ONLY_ASSEMBLY_PREVIEW_COMPLETED":
            if not 0<=x["candidate_header_tests_n"]<=2:raise ValueError("too many file probes")
            for p in x.get("candidate_header_findings",[]):
                if p.get("header_or_sample_identifiers_exported") is True:
                    raise ValueError("source header identifier leak")
                if p.get("sample_by_gene_matrix_verified") is True:
                    raise ValueError("header probe overstated matrix evidence")
    return {"schema":"hcc-stage30-bounded-header-postflight-v1",
        "candidate_studies_with_numeric_header_evidence_n":
            sum(x.get("has_tabular_numeric_header_evidence") is True for x in v["studies"]),
        "G30_B":"HOLD_COMPLETE_GENE_BY_ALIQUOT_MATRIX_NOT_VALIDATED",
        "G30_C":"HOLD_NODE_RNA_EXPRESSION_FILE_NOT_VALIDATED",
        "G30_D":"STOP_NO_SAME_PATIENT_TUMOR_RNA_PROTEOME_JOIN",
        "G30_E":"NOT_STARTED",
        "biological_replication_eligible":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True);p.add_argument("--review-out",required=True)
    a=p.parse_args()
    v={"schema":"hcc-stage30-pdc-protein-assembly-bounded-header-v1",
       "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
       "max_bytes_read_per_candidate":MAX_HEADER,
       "max_download_header_candidates_per_study":2,
       "studies":[study(k) for k in PDC]}
    out=pathlib.Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(v,indent=2)+"\n")
    rt=review(v)
    p=pathlib.Path(a.review_out);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(rt,indent=2)+"\n")
    print(json.dumps(rt,indent=2))
    print("HEADER_PREVIEW",json.dumps({x["pdc_study_id"]:{
        "status":x["status"],"assembly":x.get("assembly_file_metadata_n"),
        "signedUrl_candidates":x.get("signed_urls_allowed_host_present_n"),
        "headers_tested":x.get("candidate_header_tests_n"),
        "results":x.get("candidate_header_findings",[]),
        "failure":x.get("failure_category")} for x in v["studies"]}))
if __name__=="__main__":main()
