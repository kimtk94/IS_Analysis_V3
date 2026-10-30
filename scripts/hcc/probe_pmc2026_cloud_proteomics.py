#!/usr/bin/env python3
"""2026 NCBI PMC Cloud article-media manifest and bounded XLSX verification.

Per NCBI 2026 official PMC Cloud README, article versions are stored at
s3://pmc-oa-opendata/PMC<id>.<version> and metadata at
metadata/PMC<id>.<version>.json. Read metadata first; resolve media_urls
*from NCBI metadata*, never guess an S3 object or assume old PMC HTTP works.

No patient identifiers, gene values, raw/signed URLs or full workbooks
are persisted. Public 2023 HCC supplementary mmc3.xlsx remains only a
candidate until a validated OOXML workbook is actually downloaded.
"""
import argparse
import collections
import datetime as dt
import hashlib
import io
import json
import pathlib
import re
import urllib.error
import urllib.parse
import urllib.request
import zipfile
import xml.etree.ElementTree as ET

BASE="https://pmc-oa-opendata.s3.amazonaws.com/"
PMCID="PMC10783603"
DOI="10.1016/j.xcrm.2023.101315"
TARGET="mmc3.xlsx"
MAX_META=2*1024*1024
MAX_XLSX=24*1024*1024
VERSIONS=(1,2,3)
NS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"

def official(url):
    p=urllib.parse.urlsplit(url)
    return (p.scheme=="https" and p.hostname=="pmc-oa-opendata.s3.amazonaws.com"
            and p.path.startswith(("/metadata/"+PMCID+".","/"+PMCID+"."))
            and not p.username and not p.password)

def bounded(url,cap,timeout=12):
    if not official(url):raise ValueError("non-official cloud URL")
    req=urllib.request.Request(url,headers={
        "User-Agent":"HCC-PMC-2026-OpenData-Verification/1.0",
        "Accept":"application/json,application/octet-stream,*/*"})
    with urllib.request.urlopen(req,timeout=timeout) as res:
        if not official(res.url):raise ValueError("redirect outside official cloud")
        b=res.read(cap+1)
        if len(b)>cap:raise ValueError("size limit exceeded")
        return b,res.status,res.headers.get("Content-Type","").split(";")[0]

def parse_article_metadata(data):
    meta=json.loads(data)
    if str(meta.get("pmcid","")).upper()!=PMCID:
        raise ValueError("PMCID mismatch")
    if meta.get("doi") and str(meta["doi"]).lower()!=DOI:
        raise ValueError("DOI mismatch")
    raw=meta.get("media_urls") or []
    if not isinstance(raw,list):raise ValueError("media_urls type invalid")
    candidate_urls=[]
    xlsx_filenames=[]
    entry_shapes=collections.Counter()
    link_hosts=collections.Counter()
    public_suffixes=collections.Counter()
    dictionary_fields=set()
    for entry in raw:
        entry_shapes[type(entry).__name__]+=1
        if isinstance(entry,str):url=entry
        elif isinstance(entry,dict):
            dictionary_fields.update(str(x)[:60] for x in entry.keys())
            url=entry.get("url") or entry.get("s3_url") or entry.get("media_url") or ""
        else:continue
        link_hosts[urllib.parse.urlsplit(url).hostname or "NO_HOST"]+=1
        suffix=urllib.parse.unquote(urllib.parse.urlsplit(url).path).lower().rsplit(".",1)
        if len(suffix)==2 and len(suffix[-1])<=12:
            public_suffixes[suffix[-1]]+=1
        if not official(url):continue
        filename=urllib.parse.unquote(urllib.parse.urlsplit(url).path).rsplit("/",1)[-1]
        if re.fullmatch(r"[a-zA-Z0-9_.-]{1,120}\.xlsx", filename, re.I):
            xlsx_filenames.append(filename)
        if re.search(r"(?:^|[-_.])mmc3\.xlsx$",filename,re.I):
            candidate_urls.append(url)
    return {
        "license_code":str(meta.get("license_code") or "UNKNOWN")[:32],
        "open_access":meta.get("is_pmc_openaccess"),
        "article_version":meta.get("version"),
        "media_url_count":len(raw),
        "media_entry_type_counts":dict(entry_shapes),
        "media_entry_public_keys":sorted(dictionary_fields),
        "media_url_host_counts":dict(link_hosts),
        "media_file_extension_counts":dict(public_suffixes),
        "official_xlsx_media_filenames":sorted(set(xlsx_filenames)),
        "matching_file_count":len(candidate_urls),
        "candidate_urls_in_memory":candidate_urls,
    }

def xlsx_only_manifest(raw):
    if not zipfile.is_zipfile(io.BytesIO(raw)):
        raise ValueError("not XLSX")
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        if "xl/workbook.xml" not in z.namelist():raise ValueError("no workbook")
        if z.getinfo("xl/workbook.xml").file_size>MAX_META:raise ValueError("oversize xml")
        t=ET.fromstring(z.read("xl/workbook.xml"))
        sh=t.find(NS+"sheets")
        names=[s.get("name") for s in sh] if sh is not None else []
        return {"xlsx_verified":True,"sheet_count":len(names),"sheet_names":names,
                "patient_expression_values_or_raw_ids_logged":False}

def audit():
    r={"pmcid":PMCID,"doi":DOI,
       "target":TARGET,"status":"UNVERIFIED",
       "official_source":"https://pmc-oa-opendata.s3.amazonaws.com/README.txt",
       "metadata_attempts":[],"download_verified":False,
       "same_patient_rna_protein_join":False}
    for version in VERSIONS:
        meta_url=BASE+"metadata/"+PMCID+"."+str(version)+".json"
        attempt={"version":version}
        try:
            b,status,kind=bounded(meta_url,MAX_META)
            info=parse_article_metadata(b)
            candidates=info.pop("candidate_urls_in_memory")
            attempt.update({"http":status,"content_type":kind,**info})
            r["metadata_attempts"].append(attempt)
            if not candidates:continue
            for url in candidates:
                try:
                    raw,code,ctype=bounded(url,MAX_XLSX,timeout=22)
                    out=xlsx_only_manifest(raw)
                    r.update({"status":"S3_XLSX_MANIFEST_VERIFIED",
                              "download_verified":True,
                              "source_article_version":version,
                              "download_bytes":len(raw),
                              "sha256":hashlib.sha256(raw).hexdigest(),
                              "content_type":ctype,
                              "manifest":out})
                    return r
                except urllib.error.HTTPError as exc:
                    attempt["media_http_error_code"]=exc.code
                except Exception as exc:
                    attempt["media_error_type"]=type(exc).__name__
        except urllib.error.HTTPError as exc:
            attempt.update({"http":exc.code,"metadata_status":"HTTP_ERROR"})
            r["metadata_attempts"].append(attempt)
        except Exception as exc:
            attempt.update({"metadata_status":"UNVERIFIED","error_type":type(exc).__name__})
            r["metadata_attempts"].append(attempt)
    r["status"]="PMC_2026_CLOUD_NOT_VERIFIED"
    return r

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",default="results/hcc/feasibility/pmc2026_cloud_mmc3_audit.json")
    args=p.parse_args()
    report={"checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
            "audit":audit(),"V02":"BLOCKED"}
    path=pathlib.Path(args.out)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2,ensure_ascii=False))

if __name__=="__main__":main()
