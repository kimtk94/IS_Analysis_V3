#!/usr/bin/env python3
"""GIGASTROKE EUR 4 subtype publisher-sidecar and bounded source header check.

Do not download full outcomes, infer rsIDs from positions, or claim harmonization.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import urllib.request
from urllib.parse import urlsplit
import zlib

from scripts.ckm_stroke.verify_gigastroke_raw import (
    load_manifest, source_urls, get_text, check_sidecars,
)

MAX_READ=8192
MAX_DECOMPRESSED=32768


def header_probe(row, opener=urllib.request.urlopen, text_fetcher=None):
    raw,meta,sums=source_urls(row)
    fetch_text = text_fetcher or (lambda uri:get_text(uri,opener))
    check_sidecars(row,fetch_text(meta),fetch_text(sums))
    req=urllib.request.Request(raw,headers={
        "User-Agent":"CKM-Stroke-gigastroke-bounded-header/1.0",
        "Range":f"bytes=0-{MAX_READ-1}"
    })
    with opener(req,timeout=60) as response:
        target=urlsplit(response.geturl())
        if response.status not in (200,206) or target.hostname!="ftp.ebi.ac.uk" or target.scheme!="https":
            raise ValueError("nonpublisher source redirect/status")
        buf=response.read(MAX_READ)
        if len(buf)<24 or buf[:2]!=b"\x1f\x8b":
            raise ValueError("publisher source not a gzip header prefix")
        decompressor=zlib.decompressobj(wbits=16+zlib.MAX_WBITS)
        output=decompressor.decompress(buf,MAX_DECOMPRESSED)
        if b"\n" not in output:
            raise ValueError("first header line missing from bounded prefix")
        first=output.split(b"\n",1)[0].decode("utf-8-sig","strict").strip("\r")
        if not first or len(first)>16384:
            raise ValueError("no bounded tabular header found")
        cols=first.split("\t")
        if len(cols)<4 or len(cols)>100:
            raise ValueError("outcome header width outside expected bounds")
        return {
            "phenotype":row["phenotype"],"accession":row["accession"],
            "reported_trait":row["reported_trait"],"build":row["build"],
            "original_published_md5":row["published_raw_md5"],
            "official_yaml_and_md5_sidecar":"PASS",
            "gzip_header_columns":cols,
            "gzip_prefix_bytes_read":len(buf),
            "full_raw_md5_verified":False,
            "header_only_status":"PASS_PUBLISHER_HEADER_ONLY",
            "scientific_status":"HOLD",
        }


def audit(manifest:Path,opener=urllib.request.urlopen):
    out={"schema":"ckm-gigastroke-eur-four-bounded-outcome-header-v1",
         "scientific_status":"HOLD","ready_for_ldsc":False,
         "source_full_raw_verified":False,"results":[]}
    entries=load_manifest(manifest)
    for row in entries:
        try:
            out["results"].append(header_probe(row,opener=opener))
        except Exception as exc:
            out["results"].append({
                "phenotype":row["phenotype"],"accession":row["accession"],
                "header_only_status":"HOLD_PUBLISHER_HEADER_ACCESS",
                "error_type":type(exc).__name__,"error":str(exc)[:240],
                "scientific_status":"HOLD",
            })
    out["technical_status"]=(
        "PASS_4_PUBLISHER_HEADERS_ONLY"
        if all(x["header_only_status"]=="PASS_PUBLISHER_HEADER_ONLY" for x in out["results"])
        else "HOLD_ONE_OR_MORE_OUTCOME_HEADERS"
    )
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",type=Path,default=Path("config/ckm_stroke/gigastroke_eur_official.tsv"))
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    report=audit(args.manifest)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":report["technical_status"],
                      "outcomes":[{"phenotype":r["phenotype"],"status":r["header_only_status"],
                                   "columns":r.get("gzip_header_columns")} for r in report["results"]]},
                     ensure_ascii=False))
    if report["technical_status"]!="PASS_4_PUBLISHER_HEADERS_ONLY":
        raise SystemExit(1)


if __name__=="__main__":
    main()
