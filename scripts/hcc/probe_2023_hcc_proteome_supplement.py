#!/usr/bin/env python3
"""Bounded public 2023 HCC 152-patient proteomic XLSX access audit.

Study: Cell Reports Medicine 2023 10.1016/j.xcrm.2023.101315
Supplementary Table S2: mmc3.xlsx (reported 16.2 MB; 152 HCC cases).
This checks file bytes, OOXML ZIP integrity, sheets and metadata only.
It does NOT claim actual RNA–Protein pairing, any biological result or
patient-level data access. No patient IDs, gene values, tokens or raw signed
links are persisted.
"""
import argparse
import datetime as dt
import hashlib
import io
import json
import pathlib
import urllib.error
import urllib.parse
import urllib.request
import zipfile
import xml.etree.ElementTree as ET

CANDIDATES=(
  "https://pmc.ncbi.nlm.nih.gov/articles/instance/10783603/bin/mmc3.xlsx",
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC10783603/bin/mmc3.xlsx",
)
CAP=23*1024*1024
NS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"

def fetch(url,cap=CAP):
    p=urllib.parse.urlsplit(url)
    if p.scheme!="https" or p.hostname not in ("pmc.ncbi.nlm.nih.gov","www.ncbi.nlm.nih.gov"):
        raise ValueError("Not a public NCBI PMC URL")
    req=urllib.request.Request(url,headers={
       "User-Agent":"HCC-Research-Open-Supplement-Verification/1.0",
       "Accept":"application/octet-stream,application/vnd.openxmlformats-officedocument.spreadsheetml.sheet,*/*"})
    with urllib.request.urlopen(req,timeout=17) as resp:
        if urllib.parse.urlsplit(resp.url).hostname not in (
            "pmc.ncbi.nlm.nih.gov","www.ncbi.nlm.nih.gov","cdn.ncbi.nlm.nih.gov"):
            raise ValueError("Unexpected off-repository redirect")
        body=resp.read(cap+1)
        if len(body)>cap:raise ValueError("oversize supplementary download")
        return body,resp.status,resp.headers.get("Content-Type","").split(";")[0]

def workbook_metadata(data):
    if not zipfile.is_zipfile(io.BytesIO(data)):raise ValueError("not XLSX ZIP")
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        if "xl/workbook.xml" not in z.namelist():raise ValueError("workbook missing")
        if z.getinfo("xl/workbook.xml").file_size>1024*1024:raise ValueError("workbook manifest oversize")
        root=ET.fromstring(z.read("xl/workbook.xml"))
        sheets=root.find(NS+"sheets")
        names=[el.attrib.get("name","") for el in sheets] if sheets is not None else []
        worksheets=[p for p in z.namelist() if p.startswith("xl/worksheets/sheet") and p.endswith(".xml")]
        return {"valid_ooxml":True,"sheet_count":len(names),"sheet_names":names,
                "sheet_xml_parts":len(worksheets),
                "sheet_uncompressed_bytes":[z.getinfo(p).file_size for p in worksheets],
                "values_or_patient_ids_logged":False}

def audit():
    o={"study_doi":"10.1016/j.xcrm.2023.101315",
       "supplement":"mmc3.xlsx","publication_claim":152,
       "source_type":"publisher-archived publicly listed supplement",
       "download_verified":False,"attempts":[],"status":"HOLD"}
    for u in CANDIDATES:
        item={"archive_host":urllib.parse.urlsplit(u).hostname}
        try:
            b,http,ctype=fetch(u)
            item.update({"http_status":http,"content_type":ctype,"response_bytes":len(b)})
            manifest=workbook_metadata(b)
            o["status"]="XLSX_MANIFEST_VERIFIED"
            o["download_verified"]=True
            o["sha256"]=hashlib.sha256(b).hexdigest()
            o["manifest"]=manifest
            o["attempts"].append(item)
            break
        except urllib.error.HTTPError as exc:
            item.update({"http_status":exc.code,"error_type":"HTTPError"})
        except Exception as exc:
            item["error_type"]=type(exc).__name__
        o["attempts"].append(item)
    o["paired_rna_protein_n_independently_verified"]=None
    o["source_access_not_equal_scientific_go"]="TRUE"
    return o

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",default="results/hcc/feasibility/2023_152patient_proteome_mmc3_manifest.json")
    a=p.parse_args()
    record={"checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),"result":audit()}
    dest=pathlib.Path(a.out);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(record,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(record,indent=2,ensure_ascii=False))
if __name__=="__main__":
    main()
