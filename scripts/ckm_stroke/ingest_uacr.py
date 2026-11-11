#!/usr/bin/env python3
"""Bounded full-source CKDGen UACR pilot: official MD5, gzip CRC, standardization.

Public GWAS statistics are temporary on the runner; only audit JSON is retained.
Never treats source data integrity as proof of valid downstream MR or LDSC.
"""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import tempfile
import time
from urllib.error import HTTPError, URLError
from urllib.request import urlopen, Request

from scripts.ckm_stroke.normalize_sumstats import normalize

UACR = {
    "url": "https://ckdgen.imbi.uni-freiburg.de/files/Teumer2019/formatted_20180517-UACR_overall-EA-nstud_18-SumMac_400.tbl.rsid.gz",
    "filename": "formatted_20180517-UACR_overall-EA-nstud_18-SumMac_400.tbl.rsid.gz",
    "md5": "51270db6886892199975b343452f1d90",
    "md5_source": "https://ckdgen.imbi.uni-freiburg.de/files/Teumer2019/MD5SUM",
}
BLOCK = 2*1024*1024


def hash_file(path, algorithm):
    h=hashlib.new(algorithm)
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(BLOCK),b""):
            h.update(chunk)
    return h.hexdigest()


def stream_to_file(response, target, cap_bytes):
    size=0
    with target.open("wb") as out:
        while True:
            block=response.read(BLOCK)
            if not block: break
            size+=len(block)
            if size>cap_bytes:
                raise ValueError("source exceeds configured hard size limit")
            out.write(block)
    return size


def ingest(root, max_mb=650, opener=urlopen):
    root=Path(root)
    root.mkdir(parents=True, exist_ok=True)
    max_bytes=max_mb*1024*1024
    target=root/UACR["filename"]
    normalized=root/"UACR_EUR_GRCh37.normalized.tsv.gz"
    head=Request(UACR["url"],method="HEAD",headers={"User-Agent":"CKM-Stroke-P1-Integrity/1.0"})
    try:
        with opener(head,timeout=45) as resp:
            length=resp.headers.get("Content-Length","")
            if length and int(length)>max_bytes:
                return {"status":"HOLD","scientific_status":"NOT_RUN",
                        "reason":"declared source exceeds cap","declared_bytes":int(length)}
    except (HTTPError,URLError,OSError,ValueError) as exc:
        # Some hosts reject HEAD even when GET works. Keep a bounded GET cap.
        length=""
        head_error=str(exc)[:200]
    else:
        head_error=""
    last_error=""
    for attempt in range(1,4):
        try:
            request=Request(UACR["url"],headers={"User-Agent":"CKM-Stroke-P1-Integrity/1.0"})
            with opener(request,timeout=120) as resp:
                if resp.status!=200:
                    raise OSError("source download returned HTTP "+str(resp.status))
                nbytes=stream_to_file(resp,target,max_bytes)
            break
        except (HTTPError,URLError,OSError,ValueError) as exc:
            last_error=str(exc)[:200]
            if target.exists(): target.unlink()
            if attempt==3:
                return {"status":"HOLD","scientific_status":"NOT_RUN",
                        "reason":"download failed or exceeded size cap","details":last_error}
            time.sleep(attempt)
    if length and int(length)!=nbytes:
        return {"status":"FAIL","scientific_status":"NOT_RUN",
                "reason":"Content-Length mismatch","declared":int(length),"downloaded":nbytes}
    digest_md5=hash_file(target,"md5")
    digest_sha=hash_file(target,"sha256")
    if digest_md5.lower()!=UACR["md5"]:
        return {"status":"FAIL","scientific_status":"NOT_RUN","reason":"original MD5 mismatch",
                "expected_md5":UACR["md5"],"observed_md5":digest_md5,
                "input_sha256":digest_sha}
    try:
        results=normalize(target,normalized,"ckdgen2019",digest_sha)
    except (OSError,ValueError,EOFError) as exc:
        return {"status":"FAIL","scientific_status":"NOT_RUN",
                "reason":"normalization or gzip integrity failed","details":str(exc)[:300]}
    return {"status":"PASS" if results["status"] in ("PASS","HOLD") else "FAIL",
            "scientific_status":"HOLD","source_verified":True,
            "source_md5_verified":True,"gzip_stream_read_to_end":True,
            "source_bytes":nbytes,"official_md5":UACR["md5"],
            "observed_md5":digest_md5,"original_sha256":digest_sha,
            "normalization_status":results["status"],
            "input_rows":results["input_rows"],"normalized_rows":results["output_rows"],
            "excluded":results["excluded"],
            "excluded_total":results["excluded_total"],
            "excluded_examples":results["excluded_examples"],
            "malformed_field_count_histogram":results["malformed_field_count_histogram"],
            "recovered_missing_rsid":results["recovered_missing_rsid"],
            "recovered_examples":results["recovered_examples"],
            "normalized_sha256":results["output_sha256"],
            "note":"Original source verified against consortium MD5; full source parsed. Not yet allele harmonized, independent LDSC or GWAS duplicate QC."}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--max-mb",type=int,default=650)
    parser.add_argument("--out",default="artifacts/ckm_stroke/uacr_full_integrity.json")
    opts=parser.parse_args()
    try:
        if opts.max_mb<1 or opts.max_mb>750: raise ValueError("invalid max-mb")
        with tempfile.TemporaryDirectory(prefix="ckm_stroke_p1_") as d:
            result=ingest(d,opts.max_mb)
    except (OSError,ValueError,OverflowError) as exc:
        result={"status":"FAIL","scientific_status":"NOT_RUN","reason":str(exc)}
    result.update({"source":"CKDGen_Teumer2019","trait":"UACR","ancestry":"EUR",
                   "url":UACR["url"],"official_checksum_source":UACR["md5_source"],
                   "timestamp_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
                   "source_commit":os.getenv("GITHUB_SHA","")})
    dest=Path(opts.out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False),flush=True)
    summ=os.getenv("GITHUB_STEP_SUMMARY")
    if summ:
        with open(summ,"a",encoding="utf-8") as fh:
            fh.write("\n## CKM UACR genuine source integrity pilot\n\n")
            for key in ("status","scientific_status","source_md5_verified",
                        "source_bytes","input_rows","normalized_rows","normalization_status"):
                fh.write("- "+key+": "+str(result.get(key,"unavailable"))+"\n")
    return 1 if result["status"]=="FAIL" else 0

if __name__=="__main__":
    raise SystemExit(main())
