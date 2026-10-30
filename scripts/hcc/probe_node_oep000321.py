#!/usr/bin/env python3
"""Probe public NODE OEP000321 landing page without scraping restricted files.

Records only response metadata and counts of candidate public file references.
Does not infer public accessibility of the processed RSEM matrix from a landing
page, from third-party papers, or from a JavaScript application shell.
"""
import argparse
import datetime as dt
import json
import pathlib
import re
import urllib.error
import urllib.request

URL = "https://www.biosino.org/node/project/detail/OEP000321"
MAX_BYTES = 262144

def probe():
    result={"accession":"OEP000321","url":URL,"status":"UNVERIFIED",
            "public_processed_rna_expression_matrix_url_verified":False,
            "downloaded_expression_matrix":False,
            "rna_protein_matched_patient_count":None}
    try:
        request=urllib.request.Request(URL,headers={"User-Agent":"HCC-Open-Data-Provenance/1.0","Accept":"text/html"})
        with urllib.request.urlopen(request,timeout=16) as response:
            data=response.read(MAX_BYTES+1)
            if len(data)>MAX_BYTES:
                raise ValueError("Exceeded HTML inspection size limit")
            result["http_status"]=response.status
            result["content_type"]=response.headers.get("content-type","").split(";")[0]
            result["html_bytes_inspected"]=len(data)
            result["final_host"]=urllib.request.urlparse(response.url).hostname if hasattr(urllib.request,"urlparse") else "unrecorded"
        html=data.decode("utf-8",errors="replace")
        result["html_mentions_accession"]=bool(re.search(r"OEP000321",html,re.I))
        # Source shell can contain generic links; these are hints only.
        result["link_like_download_token_count"]=len(re.findall(r"(?:download|filelist|\\.tsv|\\.txt|\\.gz|\\.csv)",html,re.I))
        result["status"]="PUBLIC_LANDING_HTML_ONLY"
    except urllib.error.HTTPError as e:
        result["status"]="HTTP_ERROR"
        result["http_status"]=e.code
    except urllib.error.URLError:
        result["status"]="NETWORK_UNREACHABLE"
    except Exception as e:
        result["status"]="UNVERIFIED"
        result["error_type"]=type(e).__name__
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",default="results/hcc/feasibility/node_oep000321_access.json")
    args=p.parse_args()
    report={"checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
            "NODE":probe(),
            "secondary_literature_hint":"Some independent papers report downloading an RSEM gene expression matrix from OEP000321; no direct file URL or sample join has been verified.",
            "V02":"BLOCKED"}
    path=pathlib.Path(args.out)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))

if __name__=="__main__":
    main()
