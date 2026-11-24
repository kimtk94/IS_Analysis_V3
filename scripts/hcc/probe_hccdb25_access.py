#!/usr/bin/env python3
"""Public HCCDB25 (OEP000321) fallback RNA metadata/access checks.

HCCDB v2.0 archives Gao 2019 OEP000321 transcriptome as HCCDB25,
as supported by HCCDB v2.0 documentation and later publications.
Probe only canonical site roots. No inference of a working expression
download link from a web shell or unrelated supplement.
"""
import argparse
import datetime as dt
import json
import pathlib
import re
import urllib.error
import urllib.parse
import urllib.request

SITES = {
  "hccdb_v2_canonical":"https://lifeome.net/database/hccdb2",
  "hccdb_portal_https":"https://lifeome.net:809/",
  "hccdb_portal_http":"http://lifeome.net:809/",
  "hccdb_legacy_root":"http://lifeome.net/database/hccdb2",
}
LIMIT = 160 * 1024

def inspect_landing(url):
    rec={"source":url,"status":"UNVERIFIED",
         "downloaded_rna_expression":False,
         "verified_hccdb25_expression_file":False,
         "joined_pdc_aliquot_ids":False}
    try:
        req=urllib.request.Request(url,headers={"User-Agent":"HCC-Public-RNA-Landing-Access/1.0",
                                              "Accept":"text/html,*/*",
                                              "Range":f"bytes=0-{LIMIT-1}"})
        with urllib.request.urlopen(req,timeout=11) as response:
            host=urllib.parse.urlsplit(response.url).hostname
            if host not in ("lifeome.net","www.lifeome.net"):
                raise ValueError("unexpected cross-domain redirect")
            content=response.read(LIMIT+1)
            if len(content)>LIMIT:
                raise ValueError("landing exceeds bounded check size")
            rec["http_status"]=response.status
            rec["final_host"]=host
            rec["content_type"]=response.headers.get("content-type","").split(";")[0]
            rec["received_bytes"]=len(content)
        sample=content.decode("utf-8",errors="replace")
        rec["mentions_HCCDB25"]=bool(re.search(r"HCCDB25",sample,re.I))
        rec["mentions_OEP000321"]=bool(re.search(r"OEP000321",sample,re.I))
        rec["generic_download_keyword_count"]=len(re.findall(r"(?:download|repository|expression)",sample,re.I))
        rec["status"]="LANDING_HTML_ONLY"
    except urllib.error.HTTPError as e:
        rec["status"]="HTTP_ERROR"
        rec["http_status"]=e.code
    except urllib.error.URLError as e:
        rec["status"]="NETWORK_ERROR"
        rec["error_type"]=type(e.reason).__name__
    except Exception as e:
        rec["status"]="UNVERIFIED"
        rec["error_type"]=type(e).__name__
    return rec


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",default="results/hcc/feasibility/hccdb25_access.json")
    args=p.parse_args()
    report={"checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
            "archive_alias":"OEP000321 -> HCCDB25",
            "archive_alias_evidence":"https://academic.oup.com/gpb/article/22/1/qzae011/7600437",
            "documented_dataset":"HCCDB25",
            "landing_checks":{name:inspect_landing(url) for name,url in SITES.items()},
            "matrix_download_status":"NOT_VERIFIED",
            "pdc_rna_protein_join":"NOT_VERIFIED"}
    dest=pathlib.Path(args.out)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
