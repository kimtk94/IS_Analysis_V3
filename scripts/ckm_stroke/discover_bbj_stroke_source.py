#!/usr/bin/env python3
"""Discover official publicly downloadable BBJ ischemic-stroke sumstats without guessing URLs.
Catalogs are different BBJ study versions: no assumption of independent replication.
Only checks headers and download locations, not per-SNP effect.
"""
import argparse,json,urllib.request,urllib.parse,urllib.error
from html.parser import HTMLParser
from pathlib import Path

CATALOG="https://www.pheweb.jp/downloads"
OPEN="https://opengwas.io/datasets/bbj-a-129"
class Links(HTMLParser):
    def __init__(self): super().__init__();self.rows=[];self.in_row=False;self.text=[];self.links=[];self.href=None
    def handle_starttag(self,tag,attrs):
        at=dict(attrs)
        if tag=="tr":self.in_row=True;self.text=[];self.links=[]
        if self.in_row and tag=="a":self.href=at.get("href")
    def handle_data(self,data):
        if self.in_row:self.text.append(data)
    def handle_endtag(self,tag):
        if tag=="a":self.href=None
        if tag=="tr" and self.in_row:
            t=" ".join(" ".join(self.text).split())
            if "Ischemic stroke" in t:self.rows.append({"row":t[:500],"links":self.links[:]})
            self.in_row=False
    def handle_startendtag(self,tag,attrs):pass
    def handle_entityref(self,name):self.handle_data("&"+name+";")
    # capture anchor information when text is encountered
    def handle_data(self,data):
        if self.in_row:
            self.text.append(data)
            if self.href and data.strip(): self.links.append({"label":data.strip(),"url":urllib.parse.urljoin(CATALOG,self.href)})
def get(url):
    req=urllib.request.Request(url,headers={"User-Agent":"TPMI43321-BBJ-public-source-discovery/1.0"})
    with urllib.request.urlopen(req,timeout=45) as r:
        return r.read(2500000).decode("utf-8",errors="replace"),r.geturl(),r.status
def probe(url,method="GET"):
    req=urllib.request.Request(url,method=method,headers={"User-Agent":"TPMI43321-BBJ-target-probe/1.0","Accept":"application/json,text/plain,*/*"})
    try:
        with urllib.request.urlopen(req,timeout=30) as r:
            b=r.read(20000) if method=="GET" else b""
            return {"status":"PASS_HTTP","http":r.status,"url":r.geturl(),
                    "content_type":r.headers.get("Content-Type"),
                    "content_length":r.headers.get("Content-Length"),
                    "accept_ranges":r.headers.get("Accept-Ranges"),
                    "sample_prefix":b[:500].decode("utf-8",errors="replace")}
    except Exception as e:return {"status":"HOLD_HTTP","error_type":type(e).__name__,"error":str(e)[:180]}

def run():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    out={"schema":"bbj_is_stroke_source_discovery_v1",
         "target":"GRCh37:19:15573379:T:C","requested_opengwas_id":"bbj-a-129",
         "scientific_status":"HOLD","target_variant_covered":None,
         "catalogs_are_same_cohort":None,"original_allele_effect_verified":False}
    try:
        html,url,status=get(CATALOG)
        parser=Links();parser.feed(html)
        out["pheweb"]={"status":"PASS_HTTP","url":url,"http":status,"ischemic_stroke_rows":parser.rows,
                       "rows_found":len(parser.rows)}
    except Exception as e:out["pheweb"]={"status":"HOLD_HTTP","error":str(e)[:240]}
    try:
        html,url,status=get(OPEN)
        out["opengwas"]={"status":"PASS_HTTP","url":url,"http":status,
            "jwt_required_by_provider":True,"rsid_in_full_outcome_file":"NOT_QUERIED",
            "bbj_dataset_marker":("bbj-a-129" in html.lower())}
    except Exception as e:out["opengwas"]={"status":"HOLD_HTTP","error":str(e)[:240]}
    out["actual_pheweb_bbj_zip_head"]=probe("https://humandbs.dbcls.jp/files/hum0197/hum0197.v3.BBJ.IS.v1.zip","HEAD")
    out["pheweb_variant_api"]=probe("https://pheweb.jp/api/variant/19:15573379:T:C")
    out["pheweb_variant_page"]=probe("https://pheweb.jp/variant/19:15573379:T:C")
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(out,ensure_ascii=False))
if __name__=="__main__":run()
