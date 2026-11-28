#!/usr/bin/env python3
"""Bounded parallel range transfer of one public, publisher-MD5-verified GWAS gzip."""
import argparse
import concurrent.futures
import hashlib
import os
import subprocess
from pathlib import Path
SIZE=720566637
MD5="aa0056c580cc28ac77140c62dd87e079"
URL="https://zenodo.org/api/records/16976585/files/GWAS_summary_statistics.txt.gz/content"
CHUNK=32*1024*1024

def fetch(args):
    i,start,end,root=args
    p=root/f"part-{i:03d}"
    hdr=root/f"part-{i:03d}.headers"
    size=end-start+1
    for attempt in range(2):
        p.unlink(missing_ok=True)
        hdr.unlink(missing_ok=True)
        try:
            # Exact 206 and Content-Range required; abort on a server sending a full file.
            r=subprocess.run(["curl","-fLsS","--connect-timeout","15","--max-time","95",
                "--max-filesize",str(size),"--range",f"{start}-{end}","--dump-header",str(hdr),
                "--header","Accept-Encoding: identity","--output",str(p),"--write-out","%{http_code}",
                URL],capture_output=True,text=True,timeout=110)
            if r.returncode or r.stdout.strip()!="206":raise ValueError("HTTP/range transfer failure "+r.stderr[:180])
            expected=f"content-range: bytes {start}-{end}/{SIZE}"
            if expected not in hdr.read_text(errors="replace").lower().splitlines():
                raise ValueError("Unexpected Content-Range")
            if p.stat().st_size!=size:raise ValueError("Incomplete segment")
            return i,p
        except (OSError,ValueError,subprocess.TimeoutExpired) as exc:
            print(f"chunk {i} attempt {attempt+1}: {exc}",flush=True)
    raise RuntimeError(f"Chunk {i} failed bounded retries")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    root=a.output.parent/"koges_parts";root.mkdir(exist_ok=True)
    spans=[(i,start,min(SIZE,start+CHUNK)-1,root) for i,start in enumerate(range(0,SIZE,CHUNK))]
    parts={}
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:
        for i,p in ex.map(fetch,spans):
            parts[i]=p
            print("Verified byte range",i,p.stat().st_size,flush=True)
    temp=a.output.with_suffix(".partial")
    md5=hashlib.md5()
    with temp.open("wb") as out:
        for i in range(len(spans)):
            with parts[i].open("rb") as f:
                while b:=f.read(1024*1024):
                    out.write(b);md5.update(b)
    if temp.stat().st_size!=SIZE or md5.hexdigest()!=MD5:
        raise ValueError("Assembled publisher MD5 or size mismatch")
    os.replace(temp,a.output)
    print("Verified assembled full publisher file",SIZE,md5.hexdigest(),flush=True)
if __name__=="__main__":main()
