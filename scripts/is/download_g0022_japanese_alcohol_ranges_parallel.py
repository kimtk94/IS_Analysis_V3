#!/usr/bin/env python3
"""Threaded exact HTTP206 source acquisition. Full MD5 required before promotion.

Use ONLY after no other downloader process is using the same source staging.
Concurrency ≤8 and bounded requests; range chunks reusable on interruption.
"""
import argparse,concurrent.futures,hashlib,json,time
from pathlib import Path
import requests
from download_g0022_japanese_alcohol_ranges import (
    ROOT,SOURCES,BLOCK_SIZE,metadata,md5,url_for,fetch_range)

def fetch_one(job):
    name,url,start,end,total,path=job
    session=requests.Session()
    try:
        status=fetch_range(session,url,start,end,total,path,max_attempts=4)
        return start,end,status
    finally:session.close()

def acquire_parallel(source,root,workers=6,block_size=BLOCK_SIZE):
    refs=metadata(root);name=SOURCES[source];m=refs[name]
    size=int(m["size"]);check=m["checksum"].removeprefix("md5:")
    output=root/name
    if output.exists() and output.stat().st_size==size and md5(output)==check:
        print("ALREADY_MD5_VERIFIED",name,flush=True);return
    if output.exists():raise RuntimeError("Unexpected existing full file; manual reconciliation required")
    stage=root/(name+".ranges");stage.mkdir(exist_ok=True,parents=True)
    chunks=[];joblist=[]
    for start in range(0,size,block_size):
        end=min(size-1,start+block_size-1)
        part=stage/f"r{start:012d}-{end:012d}.bin"
        chunks.append((start,end,part))
        if not (part.is_file() and part.stat().st_size==end-start+1):
            joblist.append((name,url_for(name),start,end,size,part))
    print("START",source,"total_chunks",len(chunks),
          "remaining_chunks",len(joblist),"workers",workers,flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures=[pool.submit(fetch_one,job) for job in joblist]
        done=0;errors=[]
        for f in concurrent.futures.as_completed(futures):
            try:
                start,end,status=f.result()
                done+=1
                if done==1 or done%6==0 or done==len(futures):
                    print("PROGRESS",source,done,"/",len(joblist),
                          "last_range",start,end,status,flush=True)
            except Exception as e:
                errors.append(str(e))
                print("ERROR",source,str(e)[:500],flush=True)
        if errors:
            raise RuntimeError("Range retrieval incomplete; "+str(errors[:3]))
    candidate=root/(name+".assembled.part")
    with candidate.open("wb") as out:
        for start,end,part in chunks:
            if part.stat().st_size!=end-start+1:
                raise ValueError("Incomplete source part: "+part.name)
            with part.open("rb") as inp:
                for block in iter(lambda:inp.read(2**20),b""):out.write(block)
    if candidate.stat().st_size!=size:raise ValueError("Assembled file length mismatch")
    actual=md5(candidate)
    if actual!=check:raise ValueError("Source MD5 mismatch, never promoted")
    candidate.replace(output)
    verification={"file":name,"source_record":"https://zenodo.org/records/10038152",
        "size":size,"expected_md5":check,"actual_md5":actual,
        "num_chunks":len(chunks),"http_ranges_verified":True,
        "status":"OFFICIAL_COMPLETE_MD5_VERIFIED"}
    (root/(name+".source_verified.json")).write_text(json.dumps(verification,indent=2))
    print("OFFICIAL_MD5_PASS",json.dumps(verification),flush=True)
    return verification

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--file",choices=(*SOURCES,"both"),default="both")
    ap.add_argument("--root",type=Path,default=ROOT)
    ap.add_argument("--workers",type=int,default=6)
    ap.add_argument("--block-mib",type=int,default=4)
    ap.add_argument("--execute",action="store_true")
    opts=ap.parse_args()
    if opts.workers<1 or opts.workers>8 or opts.block_mib<1 or opts.block_mib>8:
        ap.error("Supported: 1–8 workers, 1–8 MiB chunks")
    chosen=list(SOURCES) if opts.file=="both" else [opts.file]
    if not opts.execute:
        print(json.dumps({"status":"PLAN_ONLY","targets":chosen,
          "workers":opts.workers,"block_mib":opts.block_mib,
          "md5_gate":"ENFORCED"}))
    else:
        for x in chosen:acquire_parallel(x,opts.root,opts.workers,opts.block_mib*2**20)
