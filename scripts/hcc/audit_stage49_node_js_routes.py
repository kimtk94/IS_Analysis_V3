#!/usr/bin/env python3
"""HCC Stage49: read-only public NODE frontend asset route census.

This detects whether the OEP000321 SPA reveals PUBLIC metadata query routes.
It does not call undocumented APIs, bypass login, download RNA or link patients.
Only bounded aggregate counts and public-source SHA256 digests are emitted.
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

ROOT="https://www.biosino.org"
PROJECT=ROOT+"/node/project/detail/OEP000321"
MAX_HTML=128*1024
MAX_JS=512*1024
MAX_SCRIPTS=8
MAX_TOTAL_JS=3*1024*1024
ALLOWED_TYPES={"application/javascript","text/javascript","application/x-javascript",
               "text/plain","application/octet-stream"}
SIGNALS=("project/detail","project","sample","experiment","download",
         "analysis","metadata","/api/","axios","fetch(","OED","OES","OEZ","RSEM")
API_LITERAL=re.compile(r"['\x22](/(?:node/)?(?:api|project|data|download)/[A-Za-z0-9_/\-{}:$]{2,130})['\x22]",re.I)
SOURCE_ID=re.compile(r"\b(?:OEP|OES|OED|OEZ|OER)\d{5,}\b",re.I)

class ScriptParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.scripts=[]
    def handle_starttag(self,tag,attrs):
        if tag.lower()=="script":
            d=dict(attrs)
            if d.get("src"):self.scripts.append(d["src"])

def asset_url(src):
    """Allow only static same-origin JS; reject credentials, traversal and query."""
    if not isinstance(src,str) or not src or src.startswith("//"):
        return None
    try:
        p=urllib.parse.urlsplit(urllib.parse.urljoin(PROJECT,src))
        if p.scheme!="https" or p.hostname!="www.biosino.org" or p.port not in (None,443):
            return None
        if p.username or p.password or p.query or p.fragment:
            return None
        if not p.path.lower().endswith(".js"):
            return None
        if any(x in (".","..") for x in urllib.parse.unquote(p.path).split("/")):
            return None
        return urllib.parse.urlunsplit((p.scheme,p.netloc,p.path,"",""))
    except (ValueError,UnicodeError):
        return None

def bounded_get(url,cap,opener=urllib.request.urlopen):
    if url!=PROJECT and (asset_url(url)!=url):
        raise ValueError("fetch URL not allowlisted")
    request=urllib.request.Request(url,headers={"User-Agent":"HCC-Stage49-NODE-public-asset-audit/1.0",
        "Range":f"bytes=0-{cap-1}","Accept":"text/html,application/javascript,*/*"})
    with opener(request,timeout=22) as r:
        final=r.geturl()
        if final!=url:
            raise ValueError("redirected source not allowed")
        status=r.status
        mime=r.headers.get("Content-Type","").split(";")[0].strip().lower()
        data=r.read(cap+1)
    if status not in (200,206):
        raise ValueError("unexpected HTTP status")
    if len(data)>cap:
        data=data[:cap]
        truncated=True
    else:
        truncated=False
    return data,{"http_status":status,"mime_type":mime,"bytes_inspected":len(data),
                 "truncated_prefix":truncated,
                 "prefix_sha256":hashlib.sha256(data).hexdigest()}

def inspect_js(raw):
    s=raw.decode("utf-8","replace")
    literals={m.group(1) for m in API_LITERAL.finditer(s)}
    # No literal strings, paths, tokens, sample IDs or endpoint URLs leave CI.
    return {"signal_counts":{k:len(re.findall(re.escape(k),s,re.I)) for k in SIGNALS},
            "public_api_path_literal_candidates_n":len(literals),
            "source_accession_like_tokens_n":len(SOURCE_ID.findall(s)),
            "OEP000321_literal_present":"OEP000321" in s,
            "executable_api_metadata_call_verified":False,
            "source_identifiers_or_endpoints_exported":False}

def audit(opener=urllib.request.urlopen):
    out={"project":"OEP000321","source":"https://www.biosino.org/node/",
         "status":"SOURCE_ACCESS_HOLD","rna_original_expression_file_verified":False,
         "official_sample_manifest_verified":False,"same_patient_PDC_RNA_join_verified":False}
    try:
        raw,meta=bounded_get(PROJECT,MAX_HTML,opener)
        if meta["mime_type"] not in ("text/html","application/xhtml+xml"):
            raise ValueError("project did not return HTML")
        parser=ScriptParser()
        parser.feed(raw.decode("utf-8","replace"))
        seen=set()
        paths=[]
        blocked=0
        for src in parser.scripts:
            u=asset_url(src)
            if not u:
                blocked+=1
            elif u not in seen:
                paths.append(u)
                seen.add(u)
        out["landing"]={"http":meta,
            "script_elements_n":len(parser.scripts),
            "allowed_unique_static_assets_n":len(paths),
            "disallowed_script_assets_n":blocked,
            "project_accession_literal_in_html":"OEP000321" in raw.decode("utf-8","replace")}
        out["assets"]=[]
        budget=MAX_TOTAL_JS
        for url in paths[:MAX_SCRIPTS]:
            if budget<1:break
            size=min(MAX_JS,budget)
            budget-=size
            try:
                js,network=bounded_get(url,size,opener)
                if network["mime_type"] not in ALLOWED_TYPES:
                    raise ValueError("JS content type not supported")
                record={"status":"PUBLIC_JS_PREFIX_AUDITED","http":network,
                        "source":inspect_js(js)}
            except Exception as exc:
                record={"status":"ASSET_ACCESS_HOLD","error_class":type(exc).__name__}
            out["assets"].append(record)
        out["status"]="PUBLIC_JS_ROUTE_CENSUS_COMPLETED"
        out["candidate_api_path_literals_total_n"]=sum(
            a.get("source",{}).get("public_api_path_literal_candidates_n",0) for a in out["assets"])
        out["metadata_api_from_source_called"]=False
    except Exception as exc:
        out["error_class"]=type(exc).__name__
    return out

def review(v):
    if v.get("rna_original_expression_file_verified") is not False or (
        v.get("same_patient_PDC_RNA_join_verified") is not False):
        raise ValueError("unsupported RNA/proteome claim")
    if any(a.get("source",{}).get("source_identifiers_or_endpoints_exported") is not False
           for a in v.get("assets",[]) if a.get("source")):
        raise ValueError("source metadata disclosure")
    return {"schema":"hcc-stage49-node-js-route-postflight-v1",
            "public_source_route_inventory_completed":v.get("status")=="PUBLIC_JS_ROUTE_CENSUS_COMPLETED",
            "possible_public_api_path_literals_n":v.get("candidate_api_path_literals_total_n",0),
            "verified_original_processed_rna_matrix":False,
            "G30_B":"HOLD_OFFICIAL_PDC_ALIQUOT_ID_JOIN",
            "G30_C":"HOLD_NODE_OEP000321_ORIGINAL_RNA_AND_SAMPLE_MANIFEST",
            "G30_D":"STOP_NO_SAME_BIOPSY_RNA_PROTEOME_JOIN",
            "biological_claims_allowed":False}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True);ap.add_argument("--review-out",required=True)
    args=ap.parse_args()
    result={"schema":"hcc-stage49-node-public-frontend-route-census-v1",
            "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
            "source":audit()}
    verdict=review(result["source"])
    for target,obj in ((args.out,result),(args.review_out,verdict)):
        dest=pathlib.Path(target);dest.parent.mkdir(parents=True,exist_ok=True)
        dest.write_text(json.dumps(obj,indent=2)+"\n",encoding="utf-8")
    print("HCC_STAGE49_NODE_PUBLIC_SOURCE",json.dumps(result))
    print("HCC_STAGE49_REVIEW",json.dumps(verdict))

if __name__=="__main__":main()
