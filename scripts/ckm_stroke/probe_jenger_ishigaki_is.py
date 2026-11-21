#!/usr/bin/env python3
"""Probe published Ishigaki BBJ original AIS source without obtaining private JWT.
Data provenance and redirects must be checked before downloading full summary stats.
"""
import argparse,json,urllib.request,urllib.error,urllib.parse,ssl
from pathlib import Path
URLS=[
 ("JENGER_published_http","http://jenger.riken.jp:8080/download/Ischemic_stroke"),
 ("JENGER_same_https","https://jenger.riken.jp:8080/download/Ischemic_stroke"),
 ("JENGER_home","https://jenger.riken.jp/en/"),
 ("NBDC_hum0014","https://humandbs-production.ddbj.nig.ac.jp/research/hum0014/v32"),
 ("OpenGWAS_2019_metadata","https://opengwas.io/datasets/bbj-a-129"),
]
def probe(url):
 data={"requested":url}
 for method in ("HEAD","GET"):
  try:
   headers={"User-Agent":"TPMI43321-2019-BBJ-original-access/1.0",
            "Range":"bytes=0-1023"} if method=="GET" else {"User-Agent":"TPMI43321-2019-BBJ-original-access/1.0"}
   req=urllib.request.Request(url,headers=headers,method=method)
   with urllib.request.urlopen(req,timeout=22) as r:
    code=r.status
    buf=r.read(2048) if method=="GET" else b""
    final=r.geturl()
    data[method]={"status":"PASS_HTTP","http":code,"url":final,
            "content_length":r.headers.get("Content-Length"),
            "content_type":r.headers.get("Content-Type"),
            "content_range":r.headers.get("Content-Range"),
            "sample_text":buf[:300].decode("utf-8","replace") if "text" in r.headers.get("Content-Type","") else None}
  except Exception as exc:
   data[method]={"status":"HOLD_SOURCE_ACCESS",
                 "error_type":type(exc).__name__,"error":str(exc)[:180]}
 return data
def main():
 p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path);a=p.parse_args()
 rows=[{"source":k,**probe(url)} for k,url in URLS]
 state={"schema":"tpmi_BBJ_Ishigaki_2019_source_probe_v1",
        "target":"GRCh37:19:15573379:T:C","sources":rows,
        "source_original_association_observed":False,
        "independence_verified":False,"scientific_status":"HOLD"}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(state,indent=2,ensure_ascii=False)+"\n")
 print(json.dumps({"sources":rows,"scientific":"HOLD"},ensure_ascii=False))
if __name__=="__main__":main()
