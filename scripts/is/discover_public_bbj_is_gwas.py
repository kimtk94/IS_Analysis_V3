#!/usr/bin/env python3
"""Public IS source download discovery on a network-enabled GitHub runner.
Metadata/header inspection only: no bulk GWAS download, no tokens, no SSH.
"""
import argparse
import html.parser
import json
import urllib.request
from urllib.parse import urljoin,urlparse

BASE="https://pheweb.jp/downloads"
class Rows(html.parser.HTMLParser):
    def __init__(self):
        super().__init__()
        self.active=False
        self.cells=[]
        self.row=[]
        self.results=[]
        self.link=None
        self.words=[]
    def handle_starttag(self,tag,attrs):
        if tag=="tr":
            self.row=[]
        if tag in ("td","th"):
            self.active=True
            self.words=[]
            self.link=None
        if tag=="a" and self.active:
            self.link=dict(attrs).get("href")
    def handle_data(self,data):
        if self.active:self.words.append(data)
    def handle_endtag(self,tag):
        if tag in ("td","th") and self.active:
            self.row.append({"text":" ".join(" ".join(self.words).split()),"href":self.link})
            self.active=False
        if tag=="tr" and self.row:
            self.results.append(self.row)

def source():
    request=urllib.request.Request(BASE,headers={"User-Agent":"IS-Research-Public-Provenance-Audit/1.0"})
    with urllib.request.urlopen(request,timeout=25) as f:
        if f.status!=200:raise RuntimeError("PheWeb downloads index non-200")
        content=f.read(3_000_000)
        if len(content)>=3_000_000:raise ValueError("HTML source too large")
    parser=Rows()
    parser.feed(content.decode("utf-8"))
    matches=[]
    for row in parser.results:
        cells=[x["text"] for x in row]
        if "Ischemic stroke" not in cells or "BBJ" not in cells:continue
        if "22,664" not in cells or "152,022" not in cells:continue
        links=[x["href"] for x in row if x["href"]]
        for link in links:
            url=urljoin(BASE,link)
            parts=urlparse(url)
            if parts.scheme!="https" or parts.hostname not in {"humandbs.dbcls.jp","pheweb.jp","www.pheweb.jp"}:
                continue
            matches.append({"study":"Sakaue2021_BBJ_IS","phenotype":"Ischemic stroke",
                            "n_cases":22664,"n_controls":152022,"url":url})
    return matches

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--json-path",required=True)
    args=p.parse_args()
    try:
        matches=source()
        data={"schema":"is-public-gwas-download-discovery-v1",
              "status":"LINK_DISCOVERED" if matches else "NOT_FOUND",
              "source_page":BASE,"matches":matches,
              "file_downloaded":False,"official_file_identity_verified":False,
              "note":"BBJ 2021 alternative, not presumed same release as archived BBJ GWAS."}
    except Exception as ex:
        data={"schema":"is-public-gwas-download-discovery-v1",
              "status":"REQUEST_FAILED","source_page":BASE,
              "error":type(ex).__name__+":"+str(ex)[:200],
              "matches":[],"file_downloaded":False,
              "official_file_identity_verified":False}
    from pathlib import Path
    dest=Path(args.json_path);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(data,ensure_ascii=False))
if __name__=="__main__":main()
