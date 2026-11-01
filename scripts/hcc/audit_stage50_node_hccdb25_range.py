#!/usr/bin/env python3
"""Stage50: bounded NODE / HCCDB25 public frontend range & route evidence audit.

Evidence scope: HTML + explicitly linked same-origin static JS only, up to a hard
byte limit per asset. No undocumented API requests, no login bypass, no source
sample names, raw expression data, per-patient IDs or candidate URLs exported.
A JavaScript endpoint literal alone is not a verified downloadable matrix.
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
    "node":{"url":"https://www.biosino.org/node/project/detail/OEP000321",
            "host":"www.biosino.org","scheme":"https","port":443},
    "hccdb25":{"url":"http://lifeome.net:809/",
               "host":"lifeome.net","scheme":"http","port":809},
}
MAX_HTML=131072
CHUNK=262144
MAX_ASSET_BYTES=3145728
MAX_ASSETS=6
MAX_TOTAL=8*1024*1024
API_LITERAL=re.compile(r"""['"](/(?:api|node|data|download|repository|project|sample)/[A-Za-z0-9_\-/{}:$]{2,125})['"]""",re.I)
HCCDB=re.compile(r"HCCDB25",re.I)
EXPRESSION=re.compile(r"RSEM|TPM|FPKM|expression|matrix|count",re.I)
SOURCE_TOKEN=re.compile(r"\b(?:OEP|OED|OES|OEZ)[0-9]{5,}\b",re.I)
CONTENT_RANGE=re.compile(r"^bytes (\d+)-(\d+)/(\d+|\*)$")
MIME_JS={"application/javascript","text/javascript","application/x-javascript",
         "text/plain","application/octet-stream","application/json"}
class ScriptLinks(HTMLParser):
    def __init__(self):
        super().__init__();self.urls=[]
    def handle_starttag(self,tag,attrs):
        d=dict(attrs)
        if tag.lower()=="script" and d.get("src"):
            self.urls.append(d["src"])

def same_origin_asset(value,site):
    if not isinstance(value,str) or not value or value.startswith("//"):
        return None
    try:
        url=urllib.parse.urljoin(site["url"],value)
        u=urllib.parse.urlsplit(url)
        if u.hostname!=site["host"] or u.scheme!=site["scheme"]:
            return None
        if (u.port or (443 if u.scheme=="https" else 80))!=site["port"]:
            return None
        if u.username or u.password or u.query or u.fragment:
            return None
        if not u.path.lower().endswith(".js"):
            return None
        parts=urllib.parse.unquote(u.path).split("/")
        if any(part in (".","..") for part in parts):
            return None
        return urllib.parse.urlunsplit((u.scheme,u.netloc,u.path,"",""))
    except (ValueError,UnicodeError):
        return None

def read_range(url,start,amount,opener=urllib.request.urlopen):
    if amount<1 or amount>CHUNK or start<0:
        raise ValueError("range request outside hard bounds")
    request=urllib.request.Request(url,headers={
        "Range":f"bytes={start}-{start+amount-1}",
        "Accept":"text/html,application/javascript,*/*",
        "User-Agent":"HCC-Stage50-public-ranges/1.0"})
    with opener(request,timeout=22) as resp:
        status=resp.status
        if resp.geturl()!=url:
            raise ValueError("redirect changed original allowed source")
        headers=resp.headers
        mime=headers.get("Content-Type","").split(";")[0].lower().strip()
        body=resp.read(amount+1)
    if status not in (200,206):
        raise ValueError("unexpected HTTP status")
    if len(body)>amount:
        if status==200 and start==0:
            # Range ignored: don't leak a huge body; record incomplete source.
            return body[:amount],None,False,mime,"RANGE_IGNORED_PARTIAL"
        raise ValueError("content exceeded requested range")
    if status==206:
        m=CONTENT_RANGE.fullmatch(headers.get("Content-Range","").strip())
        if not m or int(m.group(1))!=start or int(m.group(2))!=start+len(body)-1:
            raise ValueError("untrusted Content-Range")
        if m.group(3)=="*":raise ValueError("unknown Content-Range total")
        total=int(m.group(3))
        if total<start+len(body):raise ValueError("impossible total length")
        return body,total,start+len(body)>=total,mime,"RANGE_CONFIRMED"
    if start!=0:
        raise ValueError("server ignored nonzero range")
    # Full 200 allowed only when Content-Length proves no unseen suffix.
    try: total=int(headers.get("Content-Length",""))
    except (ValueError,TypeError):total=None
    complete=total is not None and total==len(body)
    return body,total,complete,mime,"FULL_200" if complete else "UNBOUNDED_200_PARTIAL"

def read_asset(url,limit=MAX_ASSET_BYTES,opener=urllib.request.urlopen):
    if limit<1 or limit>MAX_ASSET_BYTES:raise ValueError("invalid byte cap")
    digest=hashlib.sha256()
    received=bytearray()
    total=None
    outcomes=[]
    offset=0
    complete=False
    while offset<limit:
        n=min(CHUNK,limit-offset)
        block,size,finished,mime,state=read_range(url,offset,n,opener)
        if not block:raise ValueError("empty source range")
        if total is not None and size is not None and total!=size:
            raise ValueError("source size changed during range verification")
        total=size if size is not None else total
        if offset==0 and mime not in MIME_JS and mime not in ("text/html","application/xhtml+xml"):
            raise ValueError("unsupported source media type")
        digest.update(block)
        received.extend(block)
        offset+=len(block)
        outcomes.append(state)
        if finished:
            complete=True
            break
        if state in ("RANGE_IGNORED_PARTIAL","UNBOUNDED_200_PARTIAL"):
            break
        if len(block)<n:
            break
    return bytes(received),{
        "inspected_bytes":len(received),
        "claimed_total_bytes":total,
        "full_source_read_verified":complete,
        "bounded_source_prefix_only":not complete,
        "source_sha256_if_complete_else_prefix_sha256":digest.hexdigest(),
        "range_response_classes":dict(sorted({x:outcomes.count(x) for x in set(outcomes)}.items())),
    }

def summarize_js(raw):
    txt=raw.decode("utf-8","replace")
    return {
        "source_HCCDB25_literal_mentions_n":len(HCCDB.findall(txt)),
        "OEP000321_literal_mentions_n":txt.count("OEP000321"),
        "generic_expression_keyword_mentions_n":len(EXPRESSION.findall(txt)),
        "distinct_literal_api_paths_n":len(set(API_LITERAL.findall(txt))),
        "accession_like_token_mentions_n":len(SOURCE_TOKEN.findall(txt)),
        "request_invocation_keywords_n":sum(len(re.findall(re.escape(s),txt,re.I))
                                          for s in ("axios","fetch(","XMLHttpRequest")),
        "public_endpoint_runtime_verified":False,
        "expression_matrix_content_verified":False,
        "source_identifiers_or_endpoint_literals_exported":False,
    }

def inspect(site_name,opener=urllib.request.urlopen):
    site=SITES[site_name]
    result={"site":site_name,"status":"PUBLIC_SITE_HOLD",
            "RNA_expression_matrix_retrieved":False,
            "NODE_to_PDC_case_sample_aliquot_join_verified":False,
            "source_urls_identifiers_or_numeric_values_exported":False}
    try:
        raw,meta=read_asset(site["url"],MAX_HTML,opener)
        html=raw.decode("utf-8","replace")
        if not html.lstrip().lower().startswith(("<!doctype","<html")):
            raise ValueError("unexpected landing content")
        parser=ScriptLinks();parser.feed(html)
        urls=[];seen=set();rejected=0
        for source in parser.urls:
            normalized=same_origin_asset(source,site)
            if not normalized:rejected+=1
            elif normalized not in seen:
                urls.append(normalized);seen.add(normalized)
        result["landing"]={"http":meta,"script_tags_n":len(parser.urls),
                           "same_origin_script_assets_n":len(urls),
                           "rejected_script_sources_n":rejected,
                           "project_accession_embedded_in_html":"OEP000321" in html}
        summary=[]
        consumed=0
        for url in urls[:MAX_ASSETS]:
            budget=min(MAX_ASSET_BYTES,MAX_TOTAL-consumed)
            if budget<1:break
            try:
                js,http=read_asset(url,budget,opener)
                consumed+=len(js)
                summary.append({"status":"BOUNDED_PUBLIC_JS_AUDITED",
                                "network":http,"indicators":summarize_js(js)})
            except Exception as exc:
                summary.append({"status":"ASSET_ACCESS_HOLD","error_class":type(exc).__name__})
        result["assets"]=summary
        result["total_js_bytes_inspected"]=consumed
        result["script_asset_count_audited"]=len(summary)
        result["full_static_source_files_verified_n"]=sum(
            x.get("network",{}).get("full_source_read_verified",False) for x in summary)
        result["partial_static_source_files_n"]=sum(
            x.get("network",{}).get("bounded_source_prefix_only",False) for x in summary)
        result["explicit_literal_api_path_candidates_n"]=sum(
            x.get("indicators",{}).get("distinct_literal_api_paths_n",0) for x in summary)
        result["source_HCCDB25_literal_mentions_n"]=sum(
            x.get("indicators",{}).get("source_HCCDB25_literal_mentions_n",0) for x in summary)
        result["status"]="PUBLIC_STATIC_SOURCE_AUDITED"
    except Exception as exc:
        result["error_class"]=type(exc).__name__
    return result

def review(item):
    if item.get("RNA_expression_matrix_retrieved") is not False or (
        item.get("NODE_to_PDC_case_sample_aliquot_join_verified") is not False):
        raise ValueError("unverified multiomics claim")
    if any(x.get("indicators",{}).get("source_identifiers_or_endpoint_literals_exported") is not False
           for x in item.get("assets",[]) if x.get("indicators")):
        raise ValueError("unapproved source identifiers in report")
    return {
        "schema":"hcc-stage50-public-source-postflight-v1",
        "site":item["site"],
        "source_js_audit_completed":item.get("status")=="PUBLIC_STATIC_SOURCE_AUDITED",
        "full_js_sources_checked":item.get("full_static_source_files_verified_n",0),
        "partial_js_sources_checked":item.get("partial_static_source_files_n",0),
        "verified_processed_RNA_expression_matrix":False,
        "G30_B":"HOLD_PDC_OFFICIAL_ALIQUOT_JOIN",
        "G30_C":"HOLD_NODE_HCCDB25_ORIGINAL_RNA_MATRIX_AND_SAMPLE_METADATA",
        "G30_D":"STOP_NO_AUTHORIZED_PATIENT_MATCHED_RNA_PROTEIN",
        "biological_or_clinical_claims_allowed":False,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--site",choices=sorted(SITES),required=True)
    p.add_argument("--out",required=True)
    p.add_argument("--review-out",required=True)
    args=p.parse_args()
    source=inspect(args.site)
    item={"schema":"hcc-stage50-node-hccdb25-range-audit-v1",
          "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
          "source":source}
    verdict=review(source)
    for file,value in ((args.out,item),(args.review_out,verdict)):
        path=pathlib.Path(file);path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(value,indent=2)+"\n",encoding="utf-8")
    print("HCC_STAGE50_SOURCE",json.dumps(item))
    print("HCC_STAGE50_REVIEW",json.dumps(verdict))

if __name__=="__main__":main()
