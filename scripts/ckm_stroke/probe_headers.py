#!/usr/bin/env python3
"""Read bounded compressed prefixes; no full GWAS download or science claim."""
import argparse
import csv
import datetime as dt
import json
import os
from pathlib import Path
import zlib
from urllib.request import urlopen, Request
from urllib.parse import urlparse

HOSTS={"ckdgen.imbi.uni-freiburg.de","csg.sph.umich.edu"}
LIMIT=131072

def extract_first_lines(payload):
    if payload.startswith(b"\x1f\x8b"):
        dec=zlib.decompressobj(16+zlib.MAX_WBITS)
        buf=dec.decompress(payload, max_length=1_000_000)
        encoding="gzip"
    else:
        buf=payload
        encoding="plain"
    text=buf.decode("utf-8",errors="replace")
    lines=[x[:500] for x in text.splitlines()[:8]]
    return {"encoding":encoding,"first_lines":lines,
            "complete_first_line":("\n" in text)}

def inspect(url, opener=urlopen):
    hostname=urlparse(url).hostname
    if hostname not in HOSTS or urlparse(url).scheme!="https":
        raise ValueError("untrusted source host")
    req=Request(url,headers={"Range":f"bytes=0-{LIMIT-1}",
                             "Accept-Encoding":"identity","User-Agent":"CKM-Stroke-Bounded-Probe/1.0"})
    try:
        with opener(req,timeout=28) as res:
            final=urlparse(res.geturl())
            if final.hostname not in HOSTS:
                return {"status":"HOLD","reason":"untrusted redirect"}
            data=res.read(LIMIT)
            result=extract_first_lines(data)
            result.update({"status":"HEADER_OBSERVED" if result["complete_first_line"] else "HOLD",
                           "http_status":res.status,"bytes_read":len(data),
                           "content_range":res.headers.get("Content-Range",""),
                           "content_length":res.headers.get("Content-Length",""),
                           "accept_ranges":res.headers.get("Accept-Ranges","")})
            return result
    except (OSError,ValueError,zlib.error) as exc:
        return {"status":"SOURCE_UNAVAILABLE","reason":str(exc)[:240]}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--sources",default="config/ckm_stroke/public_sources.tsv")
    parser.add_argument("--output",default="artifacts/ckm_stroke/header_probe.json")
    args=parser.parse_args()
    with open(args.sources,newline="",encoding="utf-8") as fh:
        rows=[x for x in csv.DictReader(fh,delimiter="\t") if x["kind"]=="file"]
    results=[]
    for row in rows:
        result=inspect(row["url"])
        results.append({"trait":row["trait"],"ancestry":row["ancestry"],"source":row["source"],
                        "url":row["url"],**result})
        print(row["trait"],row["ancestry"],result["status"],result.get("first_lines",[])[:2],flush=True)
    report={"schema":"ckm-header-probe-v1",
            "timestamp_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
            "source_commit":os.getenv("GITHUB_SHA",""),
            "probes":results, "scientific_status":"NOT_RUN",
            "note":"Bounded first bytes only. Not a full gzip/MD5/schema validation."}
    out=Path(args.output)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    return 0

if __name__=="__main__":
    raise SystemExit(main())
