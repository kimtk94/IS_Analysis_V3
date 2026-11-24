#!/usr/bin/env python3
"""Audit *publisher-provided* Cheng 2023 HCC supplementary tables.

Download only publicly accessible publisher-hosted supplement S1-S11, if
available without authorization. Read OOXML manifest/dimensions only, never
output patient IDs or gene tables. No signed URLs in logs or artifacts.
"""
import argparse
import datetime as dt
import hashlib
import io
import json
import pathlib
import re
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from html import unescape

ARTICLE = "https://academic.oup.com/bib/article/24/2/bbac601/7025463"
MINIMAL = "https://oup.silverchair-cdn.com/article-minimal/7025463"
CDN_PREFIX = "https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bib/24/2/10.1093_bib_bbac601/1/"
FILE = "table_s1-s11_bbac601.xlsx"
CAP = 24*1024*1024
XML_CAP = 5*1024*1024
NS={"ss":"http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r":"http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
HEADERS={"User-Agent":"HCC-Public-Publisher-Supplement-Audit/1.0",
         "Accept":"application/octet-stream,application/vnd.openxmlformats-officedocument.spreadsheetml.sheet,text/html,*/*"}

def accept_url(url):
    p=urllib.parse.urlsplit(url)
    return (p.scheme=="https" and p.hostname in ("oup.silverchair-cdn.com","academic.oup.com")
            and "/bib/" in p.path or
            p.scheme=="https" and p.hostname=="oup.silverchair-cdn.com"
            and p.path.endswith("/"+FILE))

def read_bounded(url, limit, timeout=18):
    if not accept_url(url):
        raise ValueError("URL outside publisher allowlist")
    req=urllib.request.Request(url,headers=HEADERS)
    with urllib.request.urlopen(req,timeout=timeout) as response:
        final=response.url
        if not accept_url(final):
            raise ValueError("Unexpected publisher redirect")
        body=response.read(limit+1)
        if len(body)>limit:
            raise ValueError("Download exceeded bounded size")
        return body, response.status, response.headers.get("Content-Type","").split(";")[0]

def discover_links(html):
    text=unescape(html.decode("utf-8",errors="replace"))
    # Search only exact supplementary basename; do NOT report temporary signed URLs.
    candidates=[]
    for m in re.finditer(r'(https?://[^\s"<>]+?'+re.escape(FILE)+r'(?:\?[^\s"<>]*)?)',text,re.I):
        u=m.group(1).replace("&amp;","&")
        if accept_url(u):candidates.append(u)
    for m in re.finditer(r'(?:href|data-url)\s*=\s*["\x27]([^"\x27]*'+re.escape(FILE)+r'[^"\x27]*)["\x27]',text,re.I):
        u=urllib.parse.urljoin(MINIMAL,unescape(m.group(1)))
        if accept_url(u):candidates.append(u)
    # Prefix without signature is a fallback but may be denied by CDN.
    candidates.append(CDN_PREFIX+FILE)
    return list(dict.fromkeys(candidates))[:6]

def xlsx_manifest(data):
    if not zipfile.is_zipfile(io.BytesIO(data)):
        raise ValueError("Response not a valid XLSX zip")
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        parts=z.namelist()
        if "xl/workbook.xml" not in parts:
            raise ValueError("missing workbook.xml")
        if z.getinfo("xl/workbook.xml").file_size>XML_CAP:
            raise ValueError("oversize workbook xml")
        wb=ET.fromstring(z.read("xl/workbook.xml"))
        sheets=[]
        for el in wb.findall("ss:sheets/ss:sheet",NS):
            sheets.append({"name":el.get("name"),"visible_state":el.get("state","visible")})
        sheetparts=sorted([p for p in parts if re.fullmatch(r"xl/worksheets/sheet\d+\.xml",p)])
        dimensions=[]
        for sp in sheetparts:
            info=z.getinfo(sp)
            if info.file_size>XML_CAP:
                dimensions.append({"sheet_xml":sp.rsplit("/",1)[-1],"status":"XML_TOO_LARGE","uncompressed_bytes":info.file_size})
                continue
            root=ET.fromstring(z.read(sp))
            dim=root.find("ss:dimension",NS)
            rows=root.findall("ss:sheetData/ss:row",NS)
            dimensions.append({"sheet_xml":sp.rsplit("/",1)[-1],
                               "dimension":dim.get("ref") if dim is not None else None,
                               "row_elements":len(rows),"xml_uncompressed_bytes":info.file_size})
        return {"sheet_count":len(sheets),"sheets":sheets,"xml_sheet_schemas":dimensions,
                "contains_shared_strings":"xl/sharedStrings.xml" in parts,
                "xlsx_signature_verified":True,
                "patient_ids_exposed":False}

def audit():
    report={"article":ARTICLE,"supplement_basename":FILE,
            "download_confirmed":False,"status":"UNVERIFIED"}
    links=[]
    for url in (MINIMAL,ARTICLE):
        try:
            html,status,ctype=read_bounded(url,3*1024*1024,timeout=15)
            report.setdefault("landing_checks",[]).append({"publisher_host":urllib.parse.urlsplit(url).hostname,
                "http_status":status,"body_bytes":len(html),"supplement_name_mentioned":FILE.encode() in html})
            links+=discover_links(html)
        except Exception as exc:
            report.setdefault("landing_checks",[]).append({"publisher_host":urllib.parse.urlsplit(url).hostname,
                "error_type":type(exc).__name__})
    links+=discover_links(b"")
    links=list(dict.fromkeys(links))[:8]
    report["publisher_candidates"]=len(links)
    attempts=[]
    for url in links:
        row={"publisher_host":urllib.parse.urlsplit(url).hostname}
        try:
            data,status,ctype=read_bounded(url,CAP,timeout=20)
            row["http_status"]=status
            row["content_type"]=ctype
            row["download_bytes"]=len(data)
            manifest=xlsx_manifest(data)
            row["zip_and_workbook_valid"]=True
            report["status"]="PUBLISHER_SUPPLEMENT_WORKBOOK_VERIFIED"
            report["download_confirmed"]=True
            report["sha256"]=hashlib.sha256(data).hexdigest()
            report["workbook_manifest"]=manifest
            attempts.append(row)
            break
        except urllib.error.HTTPError as exc:
            row["http_status"]=exc.code
            row["error_type"]="HTTPError"
        except Exception as exc:
            row["error_type"]=type(exc).__name__
        attempts.append(row)
    report["download_attempts"]=attempts
    report["V02_PAIRED_RNA_PROTEIN_DATA_JOIN"]="BLOCKED"
    report["published_408_gene_result_replicated"]=False
    report["note"]="Supplementary-table metadata, if verified, is not Gao RNA or PDC protein patient-level quantification."
    return report

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",default="results/hcc/feasibility/cheng_supplement_manifest.json")
    args=p.parse_args()
    result={"checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),"audit":audit()}
    out=pathlib.Path(args.out)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
