#!/usr/bin/env python3
"""Stage52: classify exact public HCCDB/NODE frontend data-route evidence.

This is a bounded static JS analysis, NOT an authenticated API crawler, data
downloader, or evidence of a biological/patient join. Derived route templates
may be printed ONLY when they are literal, public, query-free, anonymized paths
from a same-origin asset. No expression matrix, patient ID, sample ID, cookies,
request body, credentials, or arbitrary query strings are emitted.
"""
import argparse
import datetime as dt
import hashlib
from html.parser import HTMLParser
import json
import pathlib
import re
import urllib.parse
import urllib.request

SITES={
    "hccdb25":("http://lifeome.net:809/","lifeome.net","http",809),
    "node":("https://www.biosino.org/node/project/detail/OEP000321",
            "www.biosino.org","https",443),
}
FILES=("HCCDB25_mRNA_level3.txt","HCCDB25.sample.txt","HCCDB25.patient.txt")
CAP_HTML=65536
CAP_JS=1100000
MAX_STATIC=8
ROUTE_LITERAL=re.compile(r"""(?:(?:'|"))(/(?:node/)?(?:api|download|repository|project|data|sample|analysis|metadata)(?:/[A-Za-z0-9_./:-]{1,125})?)(?:(?:'|"))""",re.I)
CHUNK_LITERAL=re.compile(r"""(?:chunk[-_][A-Za-z0-9_-]{4,65})""")
SUSPECT_ID=re.compile(r"(?:OES|OED|OEZ|OER)\d+|[0-9a-f]{8}-[0-9a-f-]{15,}",re.I)

class ScriptParser(HTMLParser):
    def __init__(self):
        super().__init__();self.scripts=[]
    def handle_starttag(self,tag,attrs):
        if tag=="script":
            value=dict(attrs).get("src")
            if value:self.scripts.append(value)

def static_js_url(site,source):
    root,host,scheme,port=SITES[site]
    if not isinstance(source,str) or not source or source.startswith("//"):return None
    try:
        url=urllib.parse.urlsplit(urllib.parse.urljoin(root,source))
        if (url.scheme,url.hostname,url.port or (443 if scheme=="https" else 80))!=(scheme,host,port):
            return None
        if url.username or url.password or url.query or url.fragment:return None
        if not url.path.endswith(".js"):return None
        if any(x in (".","..") for x in urllib.parse.unquote(url.path).split("/")):return None
        return urllib.parse.urlunsplit((scheme,url.netloc,url.path,"",""))
    except (ValueError,UnicodeError):
        return None

def public_route_only(raw):
    try:
        if len(raw)>140 or SUSPECT_ID.search(raw) or "?" in raw or "#" in raw or ".." in raw:
            return None
        p=urllib.parse.urlsplit(raw)
        if p.scheme or p.netloc or not raw.startswith("/"):return None
        if not re.fullmatch(r"/[A-Za-z0-9_/.:+-]+",raw):return None
        return raw
    except (ValueError,UnicodeError):
        return None

def scan_js(raw):
    text=raw.decode("utf-8","replace")
    paths=sorted({x for x in (public_route_only(v) for v in ROUTE_LITERAL.findall(text)) if x})
    chunks=sorted(set(CHUNK_LITERAL.findall(text)))
    return {"candidate_public_route_literals":paths[:35],
            "candidate_public_route_literals_n":len(paths),
            "chunk_module_name_hints":chunks[:30],
            "chunk_module_name_hints_n":len(chunks),
            "source_hccdb25_mentions_n":text.lower().count("hccdb25"),
            "exact_triplet_filename_references":{name:text.count(name) for name in FILES},
            "file_payload_or_authenticated_api_verified":False,
            "unfiltered_source_excerpts_exported":False}

def read_public(url,cap,opener=urllib.request.urlopen):
    req=urllib.request.Request(url,headers={
        "User-Agent":"HCC-Stage52-public-source-route-verifier/1.0",
        "Accept":"text/html,application/javascript,*/*",
        "Range":f"bytes=0-{cap-1}"})
    with opener(req,timeout=22) as r:
        if r.geturl()!=url:raise ValueError("unexpected redirect")
        if r.status not in (200,206):raise ValueError("unexpected HTTP status")
        mime=r.headers.get("Content-Type","").split(";")[0].lower()
        size_header=r.headers.get("Content-Range","")
        data=r.read(cap+1)
        status=r.status
    is_full=False
    if status==206:
        match=re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)",size_header)
        if not match or int(match.group(1))!=0:
            raise ValueError("untrusted HTTP content range")
        total=int(match.group(3))
        if int(match.group(2))+1!=len(data) or total<len(data):
            raise ValueError("inconsistent HTTP range response")
        is_full=total==len(data)
    else:
        is_full=len(data)<=cap
    if len(data)>cap:
        is_full=False
    part=data[:cap]
    return part,{"http_status":status,"mime":mime,"inspected_bytes":len(part),
                 "source_file_complete":is_full,"sha256_source_or_prefix":hashlib.sha256(part).hexdigest()}

def audit(site,opener=urllib.request.urlopen):
    root,_,_,_=SITES[site]
    report={"site":site,"status":"SOURCE_ACCESS_HOLD",
            "rna_expression_matrix_downloaded":False,
            "official_pdc_aliquot_rna_sample_join":False}
    try:
        page,metadata=read_public(root,CAP_HTML,opener)
        if metadata["mime"] not in ("text/html","application/xhtml+xml"):
            raise ValueError("unexpected HTML media type")
        parser=ScriptParser();parser.feed(page.decode("utf-8","replace"))
        urls=[]
        for candidate in parser.scripts:
            url=static_js_url(site,candidate)
            if url and url not in urls:urls.append(url)
        report["landing"]={"html":metadata,"script_references_n":len(parser.scripts),
                           "same_origin_script_assets_n":len(urls)}
        report["scripts"]=[]
        for i,url in enumerate(urls[:MAX_STATIC]):
            # Source index allows review of reproducible public bundles without exporting URLs.
            item={"index":i+1}
            try:
                js,http=read_public(url,CAP_JS,opener)
                if http["mime"] not in ("application/javascript","text/javascript",
                                       "application/x-javascript","text/plain",
                                       "application/octet-stream"):
                    raise ValueError("JS MIME invalid")
                item.update({"status":"JS_SOURCE_PREFIX_INSPECTED","http":http,"routes":scan_js(js)})
            except Exception as exc:
                item.update({"status":"JS_SOURCE_ACCESS_HOLD","error_class":type(exc).__name__})
            report["scripts"].append(item)
        report["identified_public_routes"]=sorted({p for item in report["scripts"]
            for p in item.get("routes",{}).get("candidate_public_route_literals",[])})[:45]
        report["identified_public_routes_n"]=len(report["identified_public_routes"])
        report["triplet_static_filename_mention_counts"]={n:sum(
            item.get("routes",{}).get("exact_triplet_filename_references",{}).get(n,0)
            for item in report["scripts"]) for n in FILES}
        report["partial_js_assets_n"]=sum(not item.get("http",{}).get("source_file_complete",False)
            for item in report["scripts"] if item.get("http"))
        report["status"]="PUBLIC_STATIC_ROUTE_EVIDENCE_COLLECTED"
    except Exception as exc:
        report["error_class"]=type(exc).__name__
    return report

def decision(report):
    if report.get("rna_expression_matrix_downloaded") is not False or (
        report.get("official_pdc_aliquot_rna_sample_join") is not False):
        raise ValueError("RNA or same-biopsy identification unverified")
    return {"schema":"hcc-stage52-public-source-routes-verdict-v1",
            "site":report["site"],"route_hints_n":report.get("identified_public_routes_n",0),
            "actual_file_download_verified":False,
            "G30_B":"HOLD_OFFICIAL_PDC_ALIQUOT_ID",
            "G30_C":"HOLD_HCCDB25_NODE_RNA_ORIGINAL_FILES",
            "G30_D":"STOP_NO_OFFICIAL_INDEPENDENT_RNA_PROTEIN_PATIENT_JOIN",
            "biology_results_allowed":False}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--site",required=True,choices=sorted(SITES))
    parser.add_argument("--out",required=True)
    parser.add_argument("--review-out",required=True)
    args=parser.parse_args()
    report={"schema":"hcc-stage52-public-javascript-route-evidence-v1",
            "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
            "source":audit(args.site)}
    result=decision(report["source"])
    for file,content in ((args.out,report),(args.review_out,result)):
        path=pathlib.Path(file);path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(content,indent=2)+"\n",encoding="utf-8")
    print("HCC_STAGE52_ROUTE_SOURCE",json.dumps(report))
    print("HCC_STAGE52_REVIEW",json.dumps(result))
if __name__=="__main__":main()
