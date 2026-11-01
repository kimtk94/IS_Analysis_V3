#!/usr/bin/env python3
"""PDC000198/PDC000199 Protein Assembly-only read-only file inventory.

The first page of 6K raw mass spectra/PSMs does not sample the eight/nine
Protein Assembly outputs. Query specifically by data_category via documented
PDC filesPerStudy with no signed URLs. Report aggregate file suffix, format,
name token classes and count, never file/patient identifiers.
"""
import argparse,collections,datetime as dt,json,pathlib,re
from probe_stage30_gao_catalog import PDC,pdc_query,STUDY_FIELDS

def norm_suffix(filename):
    name=(filename or "").lower()
    for ext in (".tsv.gz",".txt.gz",".csv.gz",".xlsx",".tsv",".txt",".csv",".mzid",".gz",".xml",".pdf",".rds",".zip"):
        if name.endswith(ext):return ext
    return "other"

def classify(filename):
    s=str(filename or "").lower()
    matches=[]
    for name,pattern in (
        ("protein_quant_keyword",r"protein|proteome|abundan|quant|intensit|normalized|normalised"),
        ("gene_quant_keyword",r"gene|ensembl|uniprot"),
        ("peptide_keyword",r"peptid"),
        ("summary_keyword",r"summary|result|report|assembly"),
        ("matrix_keyword",r"matrix|table|expression"),
        ("search_engine_keyword",r"maxquant|mascot|fragpipe|proteome.discoverer|pd_"),
    ):
        if re.search(pattern,s):matches.append(name)
    return matches or ["unclassified"]

def read(pid):
    answer={"pdc_study_id":pid,"status":"ACCESS_HOLD",
        "source_file_names_and_identifiers_not_exported":True,
        "signed_download_urls_not_exported":True,
        "quantitative_matrix_downloaded":False,
        "gene_aliquot_mapping_validated":False}
    try:
        x,_=pdc_query('{study(pdc_study_id:"'+pid+'"){'+STUDY_FIELDS+'}}')
        study=x.get("study")
        if isinstance(study,list):
            if len(study)!=1:raise ValueError("ambiguous study")
            study=study[0]
        if not isinstance(study,dict) or study.get("pdc_study_id")!=pid:
            raise ValueError("source study missing")
        uid=study.get("study_id")
        if not isinstance(uid,str) or not re.fullmatch(r"[0-9a-fA-F-]{32,36}",uid):
            raise ValueError("missing study version")
        query='{ filesPerStudy(study_id:"'+uid+'",data_category:"Protein Assembly",offset:0,limit:30) { file_name data_category file_type file_format file_size md5sum } }'
        v,sha=pdc_query(query)
        rows=v.get("filesPerStudy")
        if not isinstance(rows,list):raise ValueError("unexpected PDC file schema")
        if len(rows)>30:raise ValueError("too many files returned")
        # Require that filter really worked; fail if server ignores filter.
        if any(x.get("data_category")!="Protein Assembly" for x in rows):
            raise ValueError("source category filter ignored")
        suffixes=collections.Counter();formats=collections.Counter()
        filetypes=collections.Counter();lex=collections.Counter()
        observed_md5=0
        for r in rows:
            suffixes[norm_suffix(r.get("file_name"))]+=1
            formats[str(r.get("file_format") or "")[:80]]+=1
            filetypes[str(r.get("file_type") or "")[:80]]+=1
            if re.fullmatch(r"[0-9a-fA-F]{32}",str(r.get("md5sum") or "")):
                observed_md5+=1
            for v in classify(r.get("file_name")):lex[v]+=1
        answer.update({"status":"VERIFIED_ASSEMBLY_FILE_METADATA",
           "assembly_file_count_in_response":len(rows),
           "assembly_response_sha256":sha,
           "assembly_file_name_suffixes":dict(sorted(suffixes.items())),
           "assembly_file_formats":dict(sorted(formats.items())),
           "assembly_file_types":dict(sorted(filetypes.items())),
           "assembly_filename_token_classes":dict(sorted(lex.items())),
           "assembly_files_with_recorded_md5_n":observed_md5,
           "quantitative_gene_by_sample_content_verified":False})
    except Exception as e:
        answer["failure_type"]=type(e).__name__
    return answer

def review(v):
    if v.get("schema")!="hcc-stage30-gao-assembly-metadata-v1":
        raise ValueError("assembly schema drift")
    rows=v.get("studies")
    if not isinstance(rows,list) or {x.get("pdc_study_id") for x in rows}!=set(PDC):
        raise ValueError("unexpected assembly PDC targets")
    for x in rows:
        if x.get("status") not in ("VERIFIED_ASSEMBLY_FILE_METADATA","ACCESS_HOLD"):
            raise ValueError("unsupported source access claim")
        for k in ("source_file_names_and_identifiers_not_exported","signed_download_urls_not_exported"):
            if x.get(k) is not True:raise ValueError("source private identifier emitted")
        for k in ("quantitative_matrix_downloaded","gene_aliquot_mapping_validated","quantitative_gene_by_sample_content_verified"):
            if x.get(k) is True:raise ValueError("unproved matrix availability")
        if x["status"]=="VERIFIED_ASSEMBLY_FILE_METADATA":
            n=x.get("assembly_file_count_in_response")
            if type(n) is not int or not 0<=n<=30:raise ValueError("invalid candidate count")
            if type(x.get("assembly_files_with_recorded_md5_n")) is not int or not 0<=x["assembly_files_with_recorded_md5_n"]<=n:
                raise ValueError("invalid digest count")
    return {"schema":"hcc-stage30-gao-assembly-source-postflight-v1",
        "assembly_metadata_verified_study_n":sum(x["status"]=="VERIFIED_ASSEMBLY_FILE_METADATA" for x in rows),
        "source_file_candidates_not_actual_downloaded_quantitation":True,
        "G30_B":"HOLD_PROTEIN_ABUNDANCE_QUANT_MATRIX_NOT_ACTUALLY_VERIFIED",
        "G30_C":"HOLD_NODE_RNA_QUANT_MATRIX_NOT_VERIFIED",
        "G30_D":"STOP_NO_REAL_CASE_SAMPLE_ALIQUOT_JOIN",
        "G30_E":"NOT_STARTED",
        "clinical_or_biological_replication_asserted":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True);p.add_argument("--review-out",required=True)
    a=p.parse_args()
    data={"schema":"hcc-stage30-gao-assembly-metadata-v1",
          "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
          "requested_data_category":"Protein Assembly",
          "studies":[read(s) for s in PDC]}
    dest=pathlib.Path(a.out);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(data,indent=2)+"\n")
    post=review(data)
    target=pathlib.Path(a.review_out);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(post,indent=2)+"\n")
    print(json.dumps(post,indent=2))
    print("PDC_ASSEMBLY",json.dumps({x["pdc_study_id"]:{k:x.get(k) for k in
        ("status","failure_type","assembly_file_count_in_response",
         "assembly_file_name_suffixes","assembly_file_formats","assembly_file_types",
         "assembly_filename_token_classes","assembly_files_with_recorded_md5_n")}
          for x in data["studies"]},separators=(",",":")))
if __name__=="__main__":main()
