#!/usr/bin/env python3
"""Bounded actual-archive discovery for original-imported GTEx V8 vs reprocessed QTD000161.

This is NOT a GTEx-v8-vs-EBI molecular-QTL statistical comparison.
Reports official FTP asset directories, exact metadata headers and gates.
No paid requester-pays GCS API, no full genome-wide download, no SSH.
"""
import argparse
import hashlib
from html.parser import HTMLParser
import json
import re
import urllib.parse
import urllib.request
import zlib
from pathlib import Path

ROOT="https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/"
IMPORTED=ROOT+"imported/"
QTL_DIR=ROOT+"sumstats/QTS000015/QTD000161/"
QTL_URL=QTL_DIR+"QTD000161.all.tsv.gz"
MAX_INDEX=600_000
MAX_PREFIX=65_536

class IndexParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.hrefs=[]
    def handle_starttag(self,tag,attrs):
        if tag=="a":
            href=dict(attrs).get("href")
            if href: self.hrefs.append(href)

def verified_url(url):
    parsed=urllib.parse.urlsplit(url)
    base=urllib.parse.urlsplit(ROOT)
    if parsed.scheme!="https" or parsed.hostname!=base.hostname or not parsed.path.startswith(base.path):
        raise ValueError("FTP URL must stay in official eQTL folder")
    if ".." in urllib.parse.unquote(parsed.path).split("/"):
        raise ValueError("Parent traversal prohibited")
    return url

def fetch(url,limit,*,range_header=False,fetcher=None):
    verified_url(url)
    if fetcher is not None:
        data=fetcher(url)
        if len(data)>limit: raise ValueError("Synthetic or source response larger than byte budget")
        return data
    req=urllib.request.Request(url,headers={
        "User-Agent":"IS-GTEx-Original-Imported-Lineage/1.0",
        **({"Range":f"bytes=0-{limit-1}"} if range_header else {})})
    with urllib.request.urlopen(req,timeout=22) as res:
        verified_url(res.url)
        data=res.read(limit+1)
    if len(data)>limit:
        raise ValueError("Official FTP response exceeded fixed byte cap")
    return data

def assets(raw):
    p=IndexParser()
    p.feed(raw.decode("utf-8","replace"))
    result=[]
    for href in p.hrefs:
        name=urllib.parse.unquote(urllib.parse.urlsplit(href).path).split("/")[-1]
        if not name: # trailing slash path; retain full directory name
            name=urllib.parse.unquote(urllib.parse.urlsplit(href).path).rstrip("/").split("/")[-1]
        if re.fullmatch(r"[A-Za-z0-9._-]{1,220}",name) and name not in (".",".."):
            item={"name":name,"directory":href.endswith("/")}
            if item not in result:result.append(item)
    return result

def gzip_header(prefix):
    if not prefix.startswith(b"\x1f\x8b"):
        return None
    dec=zlib.decompressobj(16+zlib.MAX_WBITS)
    chunk=dec.decompress(prefix,32768)
    first=chunk.splitlines()
    if not first:return None
    return first[0].decode("utf-8-sig","replace").split("\t")

def audit(fetcher=None):
    result={"schema":"is-original-GTEx-v8-vs-EBI-QTD000161-official-asset-discovery-v1",
        "source":"EMBL-EBI eQTL Catalogue",
        "qtl_study":"QTS000015","qtl_dataset":"QTD000161",
        "tissue":"brain_cerebellar_hemisphere","donors_from_official_catalogue":175,
        "original_GTEx_v8_FastQTL_effects_compared":False,
        "GTEx_v8_original_imported_retrieved":False,
        "GTEx_v8_original_N_nominal_allpairs_retrieved":False,
        "release6_QTD_original_bytes_sha256_verified":False,
        "scientific_status":"HOLD"}
    for title,url in (("reprocessed_QTD_dataset",QTL_DIR),("GTEx_v8_imported_root",IMPORTED)):
        try:
            raw=fetch(url,MAX_INDEX,fetcher=fetcher)
            listed=assets(raw)
            result[title]={
                "directory_url":url,"directory_sha256":hashlib.sha256(raw).hexdigest(),
                "assets":listed[:100],
                "files_discovered":len(listed),
                "gate":"OFFICIAL_INDEX_DISCOVERED" if listed else "HOLD_EMPTY_INDEX",
            }
        except (OSError,ValueError,UnicodeError) as exc:
            result[title]={"gate":"HOLD_SOURCE_INDEX_INACCESSIBLE",
                            "failure_class":type(exc).__name__}
    try:
        raw=fetch(QTL_URL,MAX_PREFIX,range_header=True,fetcher=fetcher)
        header=gzip_header(raw)
        result["reprocessed_QTD_gzip_source_prefix"]={
            "file_url":QTL_URL,
            "prefix_sha256_not_full_file":hashlib.sha256(raw).hexdigest(),
            "header":header,
            "gene_and_variant_fields_observed":bool(header and
                any("gene" in name.lower() or "molecular_trait" in name.lower() for name in header) and
                any("variant" in name.lower() or "chromosome" in name.lower() for name in header)),
            "gate":"SOURCE_QTD_HEADER_FOUND_NOT_RELEASE_PIN" if header else "HOLD_UNREADABLE_GZIP_HEADER",
        }
    except (OSError,ValueError,zlib.error,UnicodeError) as exc:
        result["reprocessed_QTD_gzip_source_prefix"]={
            "gate":"HOLD_REPROCESSED_QTD_SOURCE_UNAVAILABLE",
            "failure_class":type(exc).__name__}
    folder=result.get("GTEx_v8_imported_root",{}).get("assets",[])
    directories=[i["name"] for i in folder if i["directory"]]
    # Never guess filename from a tissue label; inspect only an observed link.
    exact_original_folder=[name for name in directories if name=="GTEx_V8"]
    if len(exact_original_folder)==1:
        folder_url=IMPORTED+urllib.parse.quote(exact_original_folder[0],safe="._-")+"/"
        try:
            raw=fetch(folder_url,MAX_INDEX,fetcher=fetcher)
            listed=assets(raw)
            result["GTEx_v8_imported_child_index"]={
                "url":folder_url,"assets":listed[:100],
                "directory_sha256":hashlib.sha256(raw).hexdigest(),
                "gate":"OFFICIAL_IMPORTED_CHILD_INDEX_FOUND" if listed else "HOLD_EMPTY_INDEX",
            }
        except (OSError,ValueError) as exc:
            result["GTEx_v8_imported_child_index"]={
                "gate":"HOLD_IMPORTED_CHILD_INACCESSIBLE",
                "failure_class":type(exc).__name__}
    result["source_discovery_gate"]=(
        "PUBLIC_SOURCE_INDEXES_DISCOVERED_ORIGINAL_QTL_COMPARISON_PENDING"
        if result.get("reprocessed_QTD_dataset",{}).get("gate")=="OFFICIAL_INDEX_DISCOVERED"
        and result.get("GTEx_v8_imported_root",{}).get("gate")=="OFFICIAL_INDEX_DISCOVERED"
        else "HOLD_SOURCE_INDEX_INCOMPLETE")
    return result

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args()
    result=audit()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("IS_GTEX_IMPORTED_DISCOVERY="+json.dumps(result,sort_keys=True))
if __name__=="__main__":
    main()
