#!/usr/bin/env python3
"""Discover HCCDB25 public web-portal download API without scraping user data.

Only the HTML shell and bounded public JS asset prefixes on lifeome.net:809 are
read. This audits the existence of public download routes, not the data payload.
Do not treat strings in JavaScript as verified downloadable RNA matrices.
"""
import argparse
import datetime as dt
from html.parser import HTMLParser
import json
import pathlib
import re
import urllib.error
import urllib.parse
import urllib.request

ROOT="http://lifeome.net:809/"
HTML_CAP=64*1024
JS_CAP=700*1024
MAX_ASSETS=12

class Scripts(HTMLParser):
    def __init__(self):
        super().__init__()
        self.items=[]
    def handle_starttag(self,tag,attrs):
        d=dict(attrs)
        if tag=="script" and d.get("src"):
            self.items.append(d["src"])
        if tag=="link" and d.get("rel") in ("preload","modulepreload") and d.get("as")=="script" and d.get("href"):
            self.items.append(d["href"])


def safe_url(path):
    url=urllib.parse.urljoin(ROOT,path)
    parsed=urllib.parse.urlsplit(url)
    return url if parsed.hostname=="lifeome.net" and parsed.scheme in ("http","https") and parsed.port in (80,809,None) else None


def small_get(url,cap,timeout=10):
    req=urllib.request.Request(url,headers={"User-Agent":"HCCDB25-Public-Asset-Probe/1.0",
                                            "Accept":"text/html,application/javascript,*/*",
                                            "Range":f"bytes=0-{cap-1}"})
    with urllib.request.urlopen(req,timeout=timeout) as response:
        host=urllib.parse.urlsplit(response.url).hostname
        if host!="lifeome.net":
            raise ValueError("unexpected redirect host")
        payload=response.read(cap+1)
        meta={"status":response.status,
              "content_type":response.headers.get("content-type","").split(";")[0],
              "response_bytes_prefix":min(len(payload),cap),
              "truncated":len(payload)>cap}
    return payload[:cap],meta


def summarize_script_text(data):
    txt=data.decode("utf-8",errors="replace")
    # Only indicate token presence/counts, no potentially sensitive URLs
    motifs=("hccdb25","repository","download","expression","/api/","axios","fetch(",
            "buke","bulk","datasource","dataset")
    return {"token_count":{s:len(re.findall(re.escape(s),txt,re.I)) for s in motifs},
            "candidate_api_route_count":len(set(re.findall(r'["\x27](/(?:api|download|data)/[a-zA-Z0-9/_-]{3,100})["\x27]',txt))),
            "does_not_verify_download":True}


def probe():
    result={"status":"UNVERIFIED","public_hccdb25_rna_matrix_verified":False}
    try:
        data,http=small_get(ROOT,HTML_CAP)
        result["landing"]=http
        html=data.decode("utf-8",errors="replace")
        scripts=Scripts()
        scripts.feed(html)
        paths=list(dict.fromkeys(scripts.items))
        result["script_asset_count"]=len(paths)
        result["inspected"]=[]
        for original in paths[:MAX_ASSETS]:
            url=safe_url(original)
            rec={"asset_origin_accepted":bool(url)}
            if url:
                rec["path_suffix_js"]=urllib.parse.urlsplit(url).path.endswith(".js")
                try:
                    js,meta=small_get(url,JS_CAP)
                    rec.update(meta)
                    rec.update(summarize_script_text(js))
                    rec["status"]="PUBLIC_JS_PREFIX_SCANNED"
                except Exception as exc:
                    rec["status"]="NOT_SCANNED"
                    rec["error_type"]=type(exc).__name__
            result["inspected"].append(rec)
        result["status"]="PUBLIC_FRONTEND_ASSETS_INVENTORIED"
    except Exception as exc:
        result["error_type"]=type(exc).__name__
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",default="results/hcc/feasibility/hccdb25_frontend_routes.json")
    args=p.parse_args()
    report={"utc":dt.datetime.now(dt.timezone.utc).isoformat(),
            "target":"HCCDB25 alternate archive of OEP000321",
            "public_url":ROOT,
            "frontend":probe(),
            "verified_rna_file_url":False,
            "V02":"BLOCKED"}
    out=pathlib.Path(args.out)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
