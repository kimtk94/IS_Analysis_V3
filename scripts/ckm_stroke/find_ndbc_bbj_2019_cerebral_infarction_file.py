#!/usr/bin/env python3
"""Identify ischemic cerebral-infarction source in NBDC hum0014-v17 40 disease GWAS.
Only NBDC dataset metadata is trusted. No association inference or raw data download.
"""
from __future__ import annotations
import argparse, concurrent.futures, html, json, re, urllib.request
from pathlib import Path
IDS=[90,91,89,107,98,92,99,94,95,97,100,93,110,101,126,104,105,103,106,109,108,111,115,113,96,114,116,117,119,118,123,121,120,112,124,122,125,102,128,127]
BASE="https://humandbs-production.ddbj.nig.ac.jp/en/dataset/NHA"
def check(n):
    ident=f"NHA{n:06d}";url=BASE+f"{n:06d}"
    req=urllib.request.Request(url,headers={"User-Agent":"IS_Analysis_V3 BBJ academic audit/1.0"})
    try:
        with urllib.request.urlopen(req,timeout=30) as r:
            content=r.read(500000).decode("utf-8",errors="replace")
            if r.status!=200 or r.geturl()!=url: raise ValueError("redirect or non-200")
        text=html.unescape(re.sub(r"<[^>]+>"," ",content))
        text=" ".join(text.split())
        # Disease names in common Methods text are NOT dataset-specific evidence.
        label=re.search(r"(?:Unrestricted-access files linked to this dataset)(.{0,800}?)(?:Analysis method)",text,re.I)
        label=label.group(1) if label else ""
        disease=re.search(r"(?:Disease\s+)(.{0,130}?)(?:Sample description|gDNA extracted)",text,re.I)
        disease=disease.group(1) if disease else ""
        hit=bool(re.search(r"cerebral infarction|ischemic stroke|ischaemic stroke|I639",label+" "+disease,re.I))
        return {"id":ident,"status":"PASS_HTTP","dataset_specific_label":label[:350],
                "dataset_disease":disease[:160],"is_cerebral_infarction":hit,
                "url":url}
    except Exception as e:
        return {"id":ident,"status":"HOLD_HTTP","error":str(e)[:200],"url":url}
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
        rows=list(pool.map(check,IDS))
    hits=[r for r in rows if r.get("is_cerebral_infarction")]
    out={"schema":"ndbc_bbj_2019_cerebral_infarction_dataset_audit_v1","version":"hum0014.v17",
         "target":"GRCh37:19:15573379:T:C","tested_datasets":len(rows),
         "successful_metadata_responses":sum(r["status"]=="PASS_HTTP" for r in rows),
         "matches":hits,
         "nonresponsive_datasets":[r for r in rows if r["status"]!="PASS_HTTP"],
         "dataset_identification_status":"PASS_UNIQUE_ICD10_I639" if len(hits)==1 else "HOLD_AMBIGUOUS_OR_ABSENT",
         "original_effect_observed":False,"scientific_status":"HOLD"}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({k:v for k,v in out.items() if k not in ("nonresponsive_datasets",)},ensure_ascii=False))
if __name__=="__main__":main()
