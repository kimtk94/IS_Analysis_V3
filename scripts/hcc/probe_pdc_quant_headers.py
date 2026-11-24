#!/usr/bin/env python3
"""Inspect short public quantitative-file headers (never bulk proteomics or restricted RNA).

PDC signed download URLs are used in memory and never printed/saved. HTTP Range
is requested and reads are bounded; full files must be separately authorized.
Results describe schemas ONLY and cannot promote multi-omic pairing gates.
"""
import argparse
import datetime as dt
import json
import pathlib
import urllib.parse
import urllib.request

from probe_quant_matrix_manifest import request_graphql

PDC_FILES = {
    "PDC000198": (
        "Zhou_Liver_Cancer_Proteome.tmt11.tsv",
        "Zhou_Liver_Cancer_Proteome.sample.txt",
    ),
    "PDC000199": (
        "Zhou_Liver_Cancer_Phosphoproteome.phosphosite.tmt11.tsv",
        "Zhou_Liver_Cancer_Phosphoproteome.sample.txt",
    )
}
READ_LIMIT = 32768
ALLOWED_HOSTS = ("pdcdatastore.s3.amazonaws.com", "pdc.cancer.gov",
                 "pdc-cancer-gov.s3.amazonaws.com",
                 "d3iwtkuvwz4jtf.cloudfront.net")


def check_download_host(url):
    parsed = urllib.parse.urlsplit(url)
    if parsed.scheme != "https" or parsed.hostname not in ALLOWED_HOSTS:
        raise ValueError("Unexpected download host: blocked")
    return True


def sniff_tsv(raw):
    # First 2 lines only. Do not expose patient/sample identifiers.
    string = raw.decode("utf-8-sig", errors="replace")
    lines = string.splitlines()
    if not lines:
        return {"status":"EMPTY_RESPONSE"}
    first=lines[0]
    parts=first.split("\t")
    candidate_fields=("gene","protein","symbol","uniprot","accession",
                      "peptide","site","phospho","plex","tmt","sample","aliquot")
    # Output fixed ontology tokens, NOT original columns: a sample ID may be a TSV header.
    tags=sorted({key for label in parts[:25] for key in candidate_fields
                 if key in label.lower()})
    return {
        "status":"HEADER_SNIFFED",
        "bytes_read":len(raw),
        "header_delimiter_tabs": first.count("\t"),
        "header_columns_in_prefix":len(parts),
        "header_complete_within_range": bool("\n" in string),
        "header_biological_field_categories":tags,
        "second_row_column_count_in_prefix":len(lines[1].split("\t")) if len(lines)>1 else None,
        "no_sample_ids_persisted":True
    }


def get_file_link(study, filename):
    # No signed URL is ever serialized. Only metadata with file size and format.
    safe_filename=json.dumps(filename)
    q=('{filesPerStudy(pdc_study_id:"'+study+'" file_name:'+safe_filename+
       ' offset:0 limit:10){file_id file_name file_size file_format data_category signedUrl{url}}}')
    data=request_graphql(q).get("filesPerStudy") or []
    if isinstance(data, dict):
        data = data.get("files", []) if "files" in data else [data]
    if not isinstance(data, list):
        raise ValueError("Unexpected filesPerStudy structure")
    selected=[f for f in data if isinstance(f, dict) and f.get("file_name")==filename]
    if len(selected)!=1:
        raise ValueError(f"expected one exact file, got {len(selected)}")
    rec=selected[0]
    signed=rec.get("signedUrl") or []
    if isinstance(signed, dict):
        signed=[signed]
    elif isinstance(signed, str):
        signed=[{"url":signed}]
    if not isinstance(signed, list):
        raise ValueError("Unexpected signedUrl structure")
    link=next((x.get("url") for x in signed if isinstance(x, dict) and x.get("url")),None)
    return rec,link


def bounded_prefix(url):
    check_download_host(url)
    request=urllib.request.Request(url,headers={
        "User-Agent":"HCC-Research-Metadata-Probe/3.0",
        "Range": f"bytes=0-{READ_LIMIT-1}",
    })
    with urllib.request.urlopen(request,timeout=18) as response:
        data=response.read(READ_LIMIT+1)
        code=response.status
        partial=response.headers.get("Content-Range")
    if len(data)>READ_LIMIT:
        raise ValueError("server returned more than allowed prefix")
    return data,code,bool(partial)


def probe():
    records=[]
    for study,files in PDC_FILES.items():
        for filename in files:
            result={"study":study,"file_name":filename,"status":"UNVERIFIED",
                    "full_quantitative_matrix_verified":False,
                    "rna_protein_patient_join_verified":False}
            try:
                file,link=get_file_link(study,filename)
                result.update({"file_size":file.get("file_size"),
                               "file_format":file.get("file_format"),
                               "data_category":file.get("data_category")})
                if not link:
                    result["status"]="FILE_METADATA_ONLY_NO_SIGNED_LINK"
                else:
                    # Safe hostname only, never retain URL path/query/signature.
                    result["signed_url_hostname"]=urllib.parse.urlsplit(link).hostname
                    payload,http,partial=bounded_prefix(link)
                    result["status"]="BOUNDED_TSV_PREFIX_FETCHED"
                    result["http_status"]=http
                    result["content_range_acknowledged"]=partial
                    result["schema"]=sniff_tsv(payload)
            except Exception as exc:
                # Never serialize exception messages: signed links can contain credentials.
                result["error_type"]=type(exc).__name__
            records.append(result)
    return records


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",default="results/hcc/feasibility/pdc_quant_headers.json")
    args=parser.parse_args()
    report={"checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
            "scope":"bounded 32 KiB prefixes of 4 PUBLIC processed files",
            "records":probe(),"expression_join_status":"NOT_VERIFIED"}
    target=pathlib.Path(args.output)
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
