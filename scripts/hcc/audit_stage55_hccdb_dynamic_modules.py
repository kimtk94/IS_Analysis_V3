#!/usr/bin/env python3
"""Stage55 HCCDB25 official public SPA module/manifest source inspection.

Inspects only HCCDB's own landing page, its literal script refs and a manifest
referenced verbatim by HTML. NEVER invents an API endpoint, requests patient
records, downloads transcriptome matrices or treats a JS route as file access.
All source excerpts are reduced to safe structural counts/module hints.
"""
import argparse
import datetime as dt
import hashlib
from html.parser import HTMLParser
import json
import pathlib
import re
import urllib.error
import urllib.parse
import urllib.request

ROOT="http://lifeome.net:809/"
HOST="lifeome.net"
PORT=809
MAX_HTML=128*1024
MAX_JS=3*1024*1024
MAX_MANIFEST=128*1024
MAX_ASSETS=14
MAX_MANIFEST_ASSETS=6
TOTAL_CAP=15*1024*1024
TARGET_FILES=("HCCDB25_mRNA_level3.txt","HCCDB25.sample.txt","HCCDB25.patient.txt")
ROUTE_MARKERS=("/repository","/download","HCCDB25")
BULK_MARKERS=("buke","bulk","transcriptome","repository","download")
API_PATH=re.compile(r"""['"]((?:/|https?://lifeome[.]net:809/)(?:api|[A-Za-z0-9_.-]+/api|download|repository)/[A-Za-z0-9_./{}:$-]{1,130})['"]""",re.I)
DYNAMIC_MARKER=re.compile(r"""(?:chunk[-_][A-Za-z0-9_-]{3,70}|webpackChunkName|__webpack_require__|\.e\(|\.u\(|import\()""")
JS_PATH=re.compile(r"""(?<![A-Za-z0-9])(?:/static/js/|static/js/)[A-Za-z0-9_./-]+[.]js""")
SENSITIVE=re.compile(r"(?:OES|OED|OEZ|OEP)[0-9]{5,}|[?&#]|token|secret|bearer|/private/",re.I)

class HTMLSource(HTMLParser):
    def __init__(self):
        super().__init__();self.scripts=[];self.manifests=[]
    def handle_starttag(self,tag,attrs):
        d=dict(attrs)
        if tag=="script" and d.get("src"):self.scripts.append(d["src"])
        if tag=="link" and d.get("rel")=="manifest" and d.get("href"):
            self.manifests.append(d["href"])

def allowed_url(src,extension):
    if not isinstance(src,str) or not src or src.startswith("//"):
        return None
    try:
        u=urllib.parse.urlsplit(urllib.parse.urljoin(ROOT,src))
        if (u.scheme,u.hostname,u.port or 80)!=("http",HOST,PORT):
            return None
        if u.username or u.password or u.query or u.fragment:
            return None
        if any(x in (".","..") for x in urllib.parse.unquote(u.path).split("/")):
            return None
        if not u.path.lower().endswith(extension):return None
        return urllib.parse.urlunsplit((u.scheme,u.netloc,u.path,"",""))
    except (ValueError,UnicodeError):
        return None

def get_limited(url,limit,opener=urllib.request.urlopen):
    if url!=ROOT and not (allowed_url(url,".js")==url or allowed_url(url,".json")==url):
        raise ValueError("URL_NOT_SOURCE_LINK_ALLOWED")
    req=urllib.request.Request(url,headers={"User-Agent":"HCC-Stage55-public-HCCDB-frontend/1.0",
        "Accept":"text/html,application/javascript,application/json,*/*",
        "Range":f"bytes=0-{limit-1}"})
    with opener(req,timeout=18) as r:
        if r.geturl()!=url or r.status not in (200,206):
            raise ValueError("UNEXPECTED_SOURCE_REDIRECT_OR_STATUS")
        data=r.read(limit+1)
        mime=r.headers.get("Content-Type","").split(";")[0].strip().lower()
        range_header=r.headers.get("Content-Range","")
        length_header=r.headers.get("Content-Length","")
        status=r.status
    truncated=len(data)>limit
    body=data[:limit]
    total=None
    complete=False
    if status==206:
        m=re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)",range_header)
        if not m or int(m.group(1))!=0 or int(m.group(2))+1!=len(body):
            raise ValueError("UNTRUSTED_CONTENT_RANGE")
        total=int(m.group(3))
        complete=total==len(body) and not truncated
    else:
        try:total=int(length_header)
        except (ValueError,TypeError):total=None
        complete=total==len(body) and not truncated
    return body,{"http_status":status,"source_mime":mime,
        "bytes_inspected":len(body),"source_bytes_total":total,
        "complete_source":complete,
        "source_or_prefix_sha256":hashlib.sha256(body).hexdigest()}

def public_static_routes(blob):
    text=blob.decode("utf-8","replace")
    result={}
    for term in ROUTE_MARKERS:
        offsets=[m.start() for m in re.finditer(re.escape(term),text,re.I)]
        windows=[]
        for at in offsets[:15]:
            snippet=text[max(0,at-400):at+450]
            module_hints=sorted(set(re.findall(r"chunk[-_][A-Za-z0-9_-]{3,60}",snippet)))[:12]
            # This explicitly records local routing context only, not code or user data.
            windows.append({"offset":at,
                "lazy_import_structure_near_route":bool(re.search(
                    r"(?:\.e\(|\.u\(|import\(|webpackChunkName|chunk[-_])",snippet)),
                "module_hints_near_route":module_hints,
                "api_word_near_route":bool(re.search(r"\bapi\b|axios|fetch\(",snippet,re.I)),
                "file_extension_near_route":bool(re.search(
                    r"\.(?:txt|tsv|gct|gz|rds)\b",snippet,re.I))})
        result[term]={"literal_occurrences_n":len(offsets),"first_source_locations":windows}
    paths=[]
    for path in API_PATH.findall(text):
        if SENSITIVE.search(path):continue
        if len(path)>145:continue
        paths.append(path)
    explicit_js=sorted({x for x in JS_PATH.findall(text) if not SENSITIVE.search(x)})
    return {"route_context":result,
            "source_literal_api_path_candidates_n":len(set(paths)),
            "source_literal_api_path_samples":sorted(set(paths))[:20],
            "explicit_static_js_asset_paths":explicit_js[:20],
            "explicit_static_js_asset_paths_n":len(explicit_js),
            "generic_dynamic_module_signals_n":len(DYNAMIC_MARKER.findall(text)),
            "HCCDB25_exact_text_occurrences_n":text.lower().count("hccdb25"),
            "expected_triplet_exact_filename_occurrences":{name:text.count(name) for name in TARGET_FILES},
            "public_download_file_verified":False,
            "official_RNA_patient_sample_mapping_verified":False}

def safe_manifest_entries(payload):
    """Only explicit .js file entries from a manifest linked by official HTML."""
    if not isinstance(payload,dict):return []
    values=[]
    for k,v in payload.items():
        if k in ("files","assets","entrypoints") and isinstance(v,dict):
            values.extend(x for x in v.values() if isinstance(x,str))
        if k in ("files","assets","entrypoints","precache") and isinstance(v,list):
            for x in v:
                if isinstance(x,str):values.append(x)
                elif isinstance(x,dict):
                    for kk in ("url","src","file"):
                        if isinstance(x.get(kk),str):values.append(x[kk])
    return list(dict.fromkeys(values))[:1000]

def audit(opener=urllib.request.urlopen):
    out={"schema":"hcc-stage55-hccdb25-public-frontend-v1",
        "status":"OFFICIAL_FRONTEND_SOURCE_HOLD",
        "dataset":"HCCDB25",
        "RNA_original_matrix_downloaded":False,
        "official_NODE_PDC_same_biopsy_crosswalk":False,
        "access_control_bypassed":False}
    try:
        html,meta=get_limited(ROOT,MAX_HTML,opener)
        if meta["source_mime"] not in ("text/html","application/xhtml+xml"):
            raise ValueError("LANDING_MIME_UNEXPECTED")
        parser=HTMLSource();parser.feed(html.decode("utf-8","replace"))
        all_assets=list(dict.fromkeys(u for src in parser.scripts
                                      if (u:=allowed_url(src,".js"))))
        manifest_reports=[]
        for link in parser.manifests[:2]:
            url=allowed_url(link,".json")
            if not url:continue
            try:
                raw,network=get_limited(url,MAX_MANIFEST,opener)
                if not network["complete_source"]:
                    raise ValueError("MANIFEST_TRUNCATED")
                parsed=json.loads(raw)
                proposed=[allowed_url(x,".js") for x in safe_manifest_entries(parsed)]
                validated=[x for x in proposed if x]
                for x in validated[:MAX_MANIFEST_ASSETS]:
                    if x not in all_assets:all_assets.append(x)
                manifest_reports.append({"status":"PUBLIC_LINKED_MANIFEST_READ",
                    "http":network,"explicit_static_js_entries_n":len(validated)})
            except Exception as exc:
                manifest_reports.append({"status":"PUBLIC_MANIFEST_HOLD",
                                         "error_class":type(exc).__name__})
        out["landing"]={"http":meta,"literal_scripts_n":len(parser.scripts),
            "same_origin_scripts_n":len([x for x in parser.scripts if allowed_url(x,".js")]),
            "manifest_links_n":len(parser.manifests)}
        out["linked_public_manifest_checks"]=manifest_reports
        out["asset_checks"]=[]
        remaining=TOTAL_CAP
        for i,url in enumerate(all_assets[:MAX_ASSETS]):
            if remaining<1:break
            cap=min(MAX_JS,remaining)
            record={"origin":"literal_landing_or_manifest","position":i+1,
                    "status":"SOURCE_ACCESS_HOLD"}
            try:
                js,network=get_limited(url,cap,opener)
                remaining-=len(js)
                if network["source_mime"] not in ("text/javascript","application/javascript",
                        "application/x-javascript","text/plain","application/octet-stream"):
                    raise ValueError("SOURCE_JS_MIME_UNEXPECTED")
                record.update({"status":"SOURCE_JS_ANALYZED",
                               "http":network,"evidence":public_static_routes(js)})
            except Exception as exc:
                record["error_class"]=type(exc).__name__
            out["asset_checks"].append(record)
        out["route_occurrences_total"]={route:sum(
            asset.get("evidence",{}).get("route_context",{}).get(route,{})
                 .get("literal_occurrences_n",0) for asset in out["asset_checks"])
            for route in ROUTE_MARKERS}
        out["full_js_sources_n"]=sum(x.get("http",{}).get("complete_source",False)
                                  for x in out["asset_checks"])
        out["partial_js_sources_n"]=sum(x.get("status")=="SOURCE_JS_ANALYZED"
                      and not x.get("http",{}).get("complete_source",False)
                      for x in out["asset_checks"])
        out["download_api_templates_n"]=sum(x.get("evidence",{})
                    .get("source_literal_api_path_candidates_n",0)
                    for x in out["asset_checks"])
        out["status"]="PUBLIC_HCCDB_FRONTEND_SOURCE_INVENTORY_DONE"
    except Exception as exc:
        out["error_class"]=type(exc).__name__
    return out

def review(out):
    if out.get("RNA_original_matrix_downloaded") is not False or (
       out.get("official_NODE_PDC_same_biopsy_crosswalk") is not False):
        raise ValueError("UNSUPPORTED_BIOLOGICAL_GO")
    return {"schema":"hcc-stage55-public-frontend-postflight-v1",
            "site":"HCCDB25",
            "verified_actual_RNA_triplet_file_content":False,
            "official_source_endpoint_identified":False,
            "G30_B":"HOLD_PDC_OFFICIAL_ALIQUOT_MAP",
            "G30_C":"HOLD_HCCDB25_RNA_SOURCE_TRIPLET",
            "G30_D":"STOP_NO_SOURCE_AUTHORIZED_PATIENT_JOIN",
            "biological_claims_allowed":False}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True)
    args=ap.parse_args()
    source=audit()
    item={"schema":"hcc-stage55-hccdb25-route-inventory-v1",
          "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
          "source":source,"review":review(source)}
    path=pathlib.Path(args.out);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(item,indent=2)+"\n",encoding="utf-8")
    print("HCC_STAGE55_FRONTEND",json.dumps(item))

if __name__=="__main__":main()
