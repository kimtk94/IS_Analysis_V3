#!/usr/bin/env python3
"""Stage57: prefetch public source-linked HCCDB JavaScript, then replay SPA.

The exact public HTML and each script explicitly linked by it are fetched once
into volatile CI memory with size limit, source MIME, range and SHA-256 checks.
A fresh unauthenticated Chromium browser is fed those *same source bytes*.
Only same-origin public /static/js/*.js subsequently requested by actual HCCDB
runtime may load from network; XHR/fetch/websocket/downloads and third parties
stay blocked. This DOES NOT retrieve RNA or validate patient/PDC identity.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import datetime as dt
import hashlib
from html.parser import HTMLParser
import json
import pathlib
import re
import urllib.parse
import urllib.request

ORIGIN="http://lifeome.net:809"
LANDING=ORIGIN+"/"
MAX_HTML=128*1024
MAX_JS=5*1024*1024
MAX_SCRIPT_N=20
MAX_DYNAMIC_N=25
ROUTES=[("home","#/"),("repository_buke","#/repository?type=buke"),
        ("download","#/download")]
TRIPLET=("HCCDB25_mRNA_level3.txt","HCCDB25.sample.txt","HCCDB25.patient.txt")
RAW_ROUTE_PATH=re.compile(r"""["'](/(?:api|download|repository|data)/[A-Za-z0-9_./:$-]{2,115})["']""")

class Scripts(HTMLParser):
    def __init__(self):
        super().__init__();self.paths=[]
    def handle_starttag(self,tag,attrs):
        if tag.lower()=="script":
            val=dict(attrs).get("src")
            if val:self.paths.append(val)

def safe_static_path(value):
    try:
        u=urllib.parse.urlsplit(urllib.parse.urljoin(LANDING,value))
        if (u.scheme,u.hostname,u.port or 80)!=("http","lifeome.net",809):
            return None
        if u.username or u.password or u.query or u.fragment:return None
        path=urllib.parse.unquote(u.path)
        if not re.fullmatch(r"/[A-Za-z0-9_./-]+[.]js",path):return None
        if any(x in (".","..") for x in path.split("/")):return None
        if re.search(r"/(?:api|private)/|(?:OED|OES|OEP|OEZ)\d{5,}",path,re.I):
            return None
        return path
    except (ValueError,UnicodeError):
        return None

def safe_dynamic_path(value):
    path=safe_static_path(value)
    return path if path and re.fullmatch(r"/static/js/[A-Za-z0-9_.-]+[.]js",path) else None

def declared_script_paths(html):
    if len(html)>MAX_HTML:raise ValueError("HTML_SOURCE_CAP")
    parser=Scripts()
    parser.feed(html.decode("utf-8","replace"))
    return sorted(set(p for src in parser.paths if (p:=safe_static_path(src))))[:MAX_SCRIPT_N]

def get_public(url,cap,opener=urllib.request.urlopen):
    if url!=LANDING and safe_static_path(url) is None:
        raise ValueError("PUBLIC_SOURCE_PATH_NOT_ALLOWED")
    if url!=LANDING and not url.startswith(ORIGIN+"/"):
        raise ValueError("PUBLIC_SOURCE_ORIGIN_NOT_ALLOWED")
    req=urllib.request.Request(url,headers={
        "User-Agent":"HCC-Stage57-cached-public-source/1.0",
        "Accept":"text/html,application/javascript,*/*",
        "Range":"bytes=0-"+str(cap-1)})
    with opener(req,timeout=65) as r:
        if r.geturl()!=url or r.status not in (200,206):
            raise ValueError("PUBLIC_SOURCE_REDIRECT_OR_HTTP_STATUS")
        data=r.read(cap+1)
        mime=r.headers.get("Content-Type","").split(";")[0].strip().lower()
        length=r.headers.get("Content-Length")
        cr=r.headers.get("Content-Range","")
        status=r.status
    if len(data)>cap or not data:raise ValueError("PUBLIC_SOURCE_SIZE_LIMIT")
    if status==206:
        m=re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)",cr.strip())
        complete=(bool(m) and int(m[1])==0 and
                  int(m[2])+1==len(data) and int(m[3])==len(data))
    else:
        complete=length is not None and int(length)==len(data)
    if not complete:raise ValueError("PUBLIC_SOURCE_PARTIAL_HOLD")
    if url==LANDING and mime not in ("text/html","application/xhtml+xml"):
        raise ValueError("PUBLIC_HTML_MIME_HOLD")
    if url!=LANDING and mime not in ("application/javascript","text/javascript",
                       "application/x-javascript","text/plain","application/octet-stream"):
        raise ValueError("PUBLIC_JS_MIME_HOLD")
    return data,{"sha256":hashlib.sha256(data).hexdigest(),"bytes":len(data),
                 "mime":mime,"http_status":status,"full_source":True}

def source_cache(opener=urllib.request.urlopen):
    html,meta=get_public(LANDING,MAX_HTML,opener)
    paths=declared_script_paths(html)
    if len(paths)<3:raise ValueError("PUBLIC_HTML_STATIC_LINKS_INSUFFICIENT")
    cache={};records=[]
    def retrieve(path):
        return path,get_public(ORIGIN+path,MAX_JS,opener)
    with ThreadPoolExecutor(max_workers=min(6,len(paths))) as pool:
        futures={pool.submit(retrieve,path):path for path in paths}
        for f in as_completed(futures):
            path=futures[f]
            try:
                _,(data,info)=f.result()
                cache[path]=data
                records.append({"public_path":path,"status":"SOURCE_COMPLETE",
                                "sha256":info["sha256"],"bytes":info["bytes"]})
            except Exception as exc:
                records.append({"public_path":path,"status":"SOURCE_HOLD",
                                "error_class":type(exc).__name__})
    return {"html":html,"html_metadata":meta,"paths":paths,"cache":cache,
            "source_checks":sorted(records,key=lambda v:v["public_path"])}

def js_evidence(source,path):
    text=source.decode("utf-8","replace")
    keys=("repository","download","buke","HCCDB25","__webpack_require__","webpackJsonp")
    safe=set()
    for match in RAW_ROUTE_PATH.finditer(text):
        route=match.group(1)
        if len(route)>120 or re.search(r"(?:OED|OES|OEP|OEZ)\d+|[?&#]|token|private|patient|\d{5,}",
                                       route,re.I):
            continue
        safe.add(route)
    return {"public_js_path":path,"sha256":hashlib.sha256(source).hexdigest(),
            "bytes":len(source),
            "marker_counts":{k:len(re.findall(re.escape(k),text,re.I)) for k in keys},
            "safe_literal_route_paths":sorted(safe)[:20],
            "triplet_file_name_mentions":{name:text.count(name) for name in TRIPLET},
            "downloaded_RNA_file":False,"patient_identity_join_verified":False}

def review(out):
    if out.get("rna_original_downloaded") is not False or (
        out.get("pdc_official_patient_join_verified") is not False):
        raise ValueError("BIOLOGICAL_SOURCE_GATES_CANNOT_BE_PROMOTED")
    return {"schema":"hcc-stage57-cached-source-postflight-v1",
            "original_HCCDB25_triplet_verified":False,
            "official_PDC_RNA_patient_identity_verified":False,
            "G30_B":"HOLD_PDC_OFFICIAL_ALIQUOT_MAP",
            "G30_C":"HOLD_HCCDB25_ORIGINAL_RNA_TRIPLET",
            "G30_D":"STOP_NO_SAME_BIOPSY_RNA_PROTEIN_JOIN",
            "biological_claims_allowed":False}

def static_audit(opener=urllib.request.urlopen):
    src=source_cache(opener)
    out={"schema":"hcc-stage57-source-prefetch-v1",
         "status":"SOURCE_PREFETCH_AUDITED",
         "source_html":src["html_metadata"],
         "declared_public_js_n":len(src["paths"]),
         "complete_source_js_n":len(src["cache"]),
         "all_declared_source_js_complete":len(src["paths"])==len(src["cache"]),
         "source_checks":src["source_checks"],
         "static_route_analysis":[js_evidence(data,path)
                                   for path,data in sorted(src["cache"].items())],
         "rna_original_downloaded":False,"pdc_official_patient_join_verified":False}
    out["review"]=review(out)
    return out

def browser_audit():
    out={"schema":"hcc-stage57-cached-browser-route-v1",
         "status":"SOURCE_BROWSER_HOLD",
         "rna_original_downloaded":False,"pdc_official_patient_join_verified":False,
         "route_checks":[]}
    try:
        src=source_cache()
    except Exception as exc:
        out["prefetch_error_class"]=type(exc).__name__
        out["review"]=review(out)
        return out
    out["official_declared_js_n"]=len(src["paths"])
    out["complete_prefetched_js_n"]=len(src["cache"])
    out["source_checks"]=src["source_checks"]
    try:
        from playwright.sync_api import sync_playwright
        with sync_playwright() as p:
            browser=p.chromium.launch(headless=True)
            ctx=browser.new_context(service_workers="block",accept_downloads=False,
                                    storage_state={"cookies":[],"origins":[]})
            page=ctx.new_page()
            observed={}
            responses=[]
            blocked={"backend":0,"external":0,"unexpected_script":0,"nonessential":0}
            def allowed_route(route):
                req=route.request
                rt=req.resource_type
                url=req.url
                u=urllib.parse.urlsplit(url)
                if rt in ("xhr","fetch","websocket"):
                    blocked["backend"]+=1
                    route.abort("blockedbyclient")
                    return
                if rt=="document" and (u.scheme,u.hostname,u.port or 80)==("http","lifeome.net",809) and u.path=="/":
                    route.fulfill(status=200,content_type="text/html; charset=utf-8",
                                  body=src["html"])
                    return
                if rt=="script":
                    path=safe_static_path(url)
                    if path in src["cache"]:
                        route.fulfill(status=200,content_type="application/javascript",
                                      body=src["cache"][path])
                        return
                    if safe_dynamic_path(url):
                        route.continue_()
                        return
                    blocked["unexpected_script"]+=1
                    route.abort("blockedbyclient")
                    return
                if u.hostname!="lifeome.net" or (u.port or 80)!=809:
                    blocked["external"]+=1
                else:
                    blocked["nonessential"]+=1
                route.abort("blockedbyclient")
            ctx.route("**/*",allowed_route)
            def on_response(resp):
                if resp.request.resource_type!="script":return
                path=safe_dynamic_path(resp.url)
                if path and path not in src["cache"] and resp.status==200:
                    responses.append((path,resp))
            page.on("response",on_response)
            baseline=set()
            for label,hashroute in ROUTES:
                before=set(observed.keys())
                item={"route":label,"status":"BROWSER_ROUTE_HOLD",
                      "new_dynamic_chunks_n":0,"new_dynamic_chunks":[]}
                try:
                    nav=page.goto(LANDING+hashroute,wait_until="domcontentloaded",timeout=30000)
                    item["document_status"]=nav.status if nav else 200
                except Exception as exc:
                    item["navigation_error_class"]=type(exc).__name__
                page.wait_for_timeout(17000)
                item["ready_state"]=page.evaluate("document.readyState")
                item["app_element_exists"]=page.locator("#app").count()>0
                for path,resp in list(responses):
                    if path in observed:continue
                    try:
                        headers=resp.all_headers()
                        if int(headers.get("content-length","0"))>MAX_JS:continue
                        data=resp.body()
                        if len(data)>MAX_JS:continue
                        observed[path]=js_evidence(data,path)
                    except Exception:
                        continue
                new=sorted(set(observed)-before-baseline)
                baseline.update(new)
                item["new_dynamic_chunks"]=new
                item["new_dynamic_chunks_n"]=len(new)
                item["blocked_backend_requests_cumulative"]=blocked["backend"]
                item["status"]=("BROWSER_ROUTE_OBSERVED" if item["ready_state"]=="complete"
                    else "BROWSER_DOCUMENT_LOADING_HOLD")
                out["route_checks"].append(item)
            out["dynamic_script_evidence"]=[observed[k] for k in sorted(observed)]
            out["blocked_network_counts"]=blocked
            out["status"]="CACHED_BROWSER_RUN_COMPLETE"
            ctx.close();browser.close()
    except Exception as exc:
        out["status"]="BROWSER_EXECUTION_HOLD"
        out["browser_error_class"]=type(exc).__name__
    out["review"]=review(out)
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mode",choices=("source","browser"),required=True)
    ap.add_argument("--out",required=True)
    args=ap.parse_args()
    try:
        report=static_audit() if args.mode=="source" else browser_audit()
    except Exception as exc:
        report={"schema":"hcc-stage57-source-execution-fail-closed-v1",
                "status":"SOURCE_EXECUTION_HOLD","error_class":type(exc).__name__,
                "rna_original_downloaded":False,"pdc_official_patient_join_verified":False}
        report["review"]=review(report)
    report["checked_utc"]=dt.datetime.now(dt.timezone.utc).isoformat()
    target=pathlib.Path(args.out)
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print("HCC_STAGE57_"+args.mode.upper(),json.dumps(report))
if __name__=="__main__":main()
