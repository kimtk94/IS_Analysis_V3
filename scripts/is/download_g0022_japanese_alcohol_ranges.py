#!/usr/bin/env python3
"""Reliable HTTP 206 range acquisition from Zenodo's public RECORD endpoint.

The Zenodo /api/.../content endpoint can terminate a large transfer and does
NOT support ranged recovery. The public /records/.../files/... endpoint supports
byte-range GET (HTTP 206). Each block must match exact Content-Range+length;
assemble with original metadata MD5. Default is plan-only; never promote a
partial or incorrect file. Data go outside Git, in separate evidence staging.
"""
import argparse, hashlib, json, os, re, time
from pathlib import Path
from urllib.parse import quote
import requests

ROOT=Path("/srv/is-analysis/data/is/gwas_exposure/japanese_alcohol_2024")
SOURCES={
 "alcohol":"1_Alcohol_intake_Unstratified.tsv.gz",
 "drinking":"5_Drinking_Unstratified.tsv.gz",
}
BLOCK_SIZE=4*1024*1024

def url_for(name):
    if name not in SOURCES.values():
        raise ValueError("Official filename not allowed")
    return "https://zenodo.org/records/10038152/files/"+quote(name)+"?download=1"

def metadata(root):
    x=json.loads((root/"ZENODO_10038152_METADATA.json").read_text())
    if str(x.get("id"))!="10038152":raise ValueError("Wrong official Zenodo accession")
    return {row["key"]:row for row in x["files"]}

def md5(path):
    digest=hashlib.md5()
    with Path(path).open("rb") as f:
        for blob in iter(lambda:f.read(2**20),b""): digest.update(blob)
    return digest.hexdigest()

def fetch_range(session,url,start,end,total,output,max_attempts=4):
    expected=end-start+1
    if output.is_file() and output.stat().st_size==expected:
        return "ALREADY_STAGED"
    if output.exists():
        output.unlink() # only within predetermined untracked staging
    last=None
    for attempt in range(max_attempts):
        temp=Path(str(output)+".part")
        try:
            with session.get(url,headers={"Range":f"bytes={start}-{end}"},
                             stream=True,timeout=(12,50)) as response:
                if response.status_code!=206:
                    raise ValueError(f"HTTP {response.status_code}, requires exact 206")
                hdr=response.headers.get("Content-Range","")
                m=re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)",hdr)
                if not m or tuple(map(int,m.groups()))!=(start,end,total):
                    raise ValueError("Range response mismatch "+hdr)
                num=response.headers.get("Content-Length")
                if num and int(num)!=expected:raise ValueError("Content-Length mismatch")
                n=0
                with temp.open("wb") as out:
                    for part in response.iter_content(chunk_size=2**18):
                        if not part:continue
                        n+=len(part)
                        if n>expected:raise ValueError("Range payload exceeds requested length")
                        out.write(part)
                if n!=expected:raise IOError(f"Range truncated {n}/{expected}")
                temp.replace(output)
                return "NEW"
        except (requests.RequestException,OSError,ValueError) as e:
            last=e
            if temp.exists():temp.unlink()
            if attempt<max_attempts-1:time.sleep(min(3,attempt+1))
    raise RuntimeError(f"Unable to retrieve byte range {start}-{end}: {last}")

def acquire(source,root,block_size=BLOCK_SIZE,session=None):
    refs=metadata(root)
    name=SOURCES[source]
    row=refs[name]
    expected_md5=row["checksum"].removeprefix("md5:")
    size=int(row["size"])
    if len(expected_md5)!=32 or size<1:
        raise ValueError("Unverifiable official source metadata")
    output=root/name
    if output.is_file() and output.stat().st_size==size and md5(output)==expected_md5:
        return {"source":source,"file":name,"status":"MD5_VERIFIED_REUSED","bytes":size}
    if output.exists():
        raise RuntimeError("A nonmatching complete-source candidate exists; manual inspection required")
    stage=root/(name+".ranges")
    stage.mkdir(parents=True,exist_ok=True)
    if session is None:session=requests.Session()
    url=url_for(name)
    stages=[]
    for start in range(0,size,block_size):
        end=min(size-1,start+block_size-1)
        part=stage/(f"r{start:012d}-{end:012d}.bin")
        stages.append((start,end,part))
        mode=fetch_range(session,url,start,end,size,part)
        print(source,"RANGE",start,end,mode,flush=True)
    temp=root/(name+".assembled.part")
    with temp.open("wb") as out:
        for start,end,part in stages:
            if part.stat().st_size!=end-start+1:
                raise RuntimeError("Staged block has invalid byte count")
            with part.open("rb") as inp:
                for blob in iter(lambda:inp.read(2**20),b""):out.write(blob)
    if temp.stat().st_size!=size:
        raise RuntimeError("Final assembled byte count invalid")
    got=md5(temp)
    if got!=expected_md5:
        raise RuntimeError(f"Source MD5 mismatch got {got} expected {expected_md5}")
    temp.replace(output)
    result={"source":source,"file":name,"status":"OFFICIAL_FULL_SOURCE_MD5_VERIFIED",
            "bytes":size,"expected_md5":expected_md5,
            "source_url":url,"blocks":len(stages)}
    (root/(name+".source_verified.json")).write_text(json.dumps(result,indent=2))
    print("SUCCESS",json.dumps(result),flush=True)
    return result

def cli():
    parser=argparse.ArgumentParser()
    parser.add_argument("--root",type=Path,default=ROOT)
    parser.add_argument("--file",choices=(*SOURCES,"both"),default="both")
    parser.add_argument("--block-mib",type=int,default=4)
    parser.add_argument("--execute",action="store_true")
    a=parser.parse_args()
    if a.block_mib<1 or a.block_mib>8:parser.error("Choose 1–8 MiB blocks")
    targets=list(SOURCES) if a.file=="both" else [a.file]
    refs=metadata(a.root)
    if not a.execute:
        print(json.dumps({"mode":"PLAN_ONLY","strategy":"EXACT_HTTP206_CHUNK_REASSEMBLY_FULL_MD5",
          "block_size_mib":a.block_mib,
          "sources":[{"type":x,"filename":SOURCES[x],
              "bytes":refs[SOURCES[x]]["size"],"md5":refs[SOURCES[x]]["checksum"]} for x in targets]},indent=2))
        return
    for kind in targets:
        acquire(kind,a.root,a.block_mib*2**20)

if __name__=="__main__":cli()
