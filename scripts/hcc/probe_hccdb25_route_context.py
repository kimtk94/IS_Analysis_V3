#!/usr/bin/env python3
"""Find the actual HCCDB2 SPA modules for /repository and /download.

Reads a public ~9 KB app.js bundle only. Captures bounded, printable JS excerpts
near known public ROUTE labels (not data/API credentials). Does not access patient
records, invoke downloads, or claim a validated HCCDB25 expression file.
"""
import argparse
import datetime as dt
import json
import pathlib
import re
from probe_hccdb25_assets import small_get, safe_url

APP="/static/js/app.3f5b30f0.js"
MAX_ASSET=64*1024
MARKERS=("/repository","/download","repository","Repository","buke","chunk-")
MAX_WINDOWS=18
CONTEXT=180


def context_windows(js):
    txt=js.decode("utf-8","replace")
    out=[]
    for term in MARKERS:
        for match in list(re.finditer(re.escape(term),txt,re.I))[:5]:
            start=max(0,match.start()-CONTEXT)
            stop=min(len(txt),match.end()+CONTEXT)
            segment=txt[start:stop]
            # Only allow bounded legible metadata; avoid leaking arbitrary query strings,
            # location identifiers, or tokens. Public paths and static chunk names only.
            static_chunks=list(dict.fromkeys(
                re.findall(r"chunk[-_][a-z0-9_-]{2,55}",segment,re.I)))[:10]
            module_tokens=list(dict.fromkeys(
                re.findall(r"(?:webpackChunkName|component|import|require|children|path|name)",segment,re.I)))[:10]
            out.append({
                "marker":term,
                "offset_in_public_bundle":match.start(),
                "nearby_chunk_names":static_chunks,
                "nearby_structure_tokens":module_tokens,
                "context_has_module_loader":bool(re.search(r"(?:\\.e\\(|\\.bind\\(|import\\(|webpackChunkName)",segment)),
                "context_has_direct_data_file_extension":bool(re.search(r"\\.(?:csv|tsv|xlsx|rds|gz|gct)\\b",segment,re.I)),
                "context_has_api_url_literal":bool(re.search(r"['\"]https?://[^'\"]+/api/",segment)),
            })
            if len(out)>=MAX_WINDOWS:
                return out
    return out


def probe():
    result={
        "status":"UNVERIFIED",
        "hccdb25_expression_file_verified":False,
        "public_matrix_download_performed":False,
        "pdc_case_join_verified":False,
    }
    url=safe_url(APP)
    if not url:
        result["status"]="UNSAFE_HOST"
        return result
    try:
        raw,net=small_get(url,MAX_ASSET,timeout=12)
        result["status"]="PUBLIC_ROUTE_MODULE_INSPECTED"
        result["app_js_path"]=APP
        result["bounded_bytes"]=len(raw)
        result["network"]=net
        result["route_windows"]=context_windows(raw)
        result["route_modules_observed"]=len(result["route_windows"])>0
    except Exception as exc:
        result["error_type"]=type(exc).__name__
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",default="results/hcc/feasibility/hccdb25_route_context.json")
    args=p.parse_args()
    result={
        "timestamp_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
        "official_literature_route":"http://lifeome.net:809/#/repository?type=buke",
        "route_analysis":probe(),
        "V02":"BLOCKED",
        "note":"Static SPA route identification is not evidence that HCCDB25 expression matrix or patient identifiers have been downloaded.",
    }
    out=pathlib.Path(args.out)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
