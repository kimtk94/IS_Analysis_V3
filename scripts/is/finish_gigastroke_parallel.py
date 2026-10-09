#!/usr/bin/env python3
"""Finish an interrupted MD5-gated GWAS download using checked concurrent HTTP ranges.
Only runs with --execute. Requires single-writer owner to be stopped.
"""
import argparse, csv, hashlib, json, re, subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

MANIFEST=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/GIGASTROKE_SOURCE_VERIFIED_V3.tsv")
DEST=Path("/srv/is-analysis/data/is/reference/gigastroke/additional")
def acquire(args):
    with args.manifest.open() as f:
        hit=[r for r in csv.DictReader(f,delimiter="\t") if r["accession"]==args.accession]
    if len(hit)!=1:raise ValueError("Accession not verified")
    r=hit[0]
    if not re.fullmatch(r"[a-fA-F0-9]{32}",r["source_md5"]):raise ValueError("Invalid source checksum")
    url=r["source_url"]
    if not url.startswith("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"):
        raise ValueError("Unapproved download source")
    head=subprocess.run(["curl","-fsSIL","--max-time","20",url],
        text=True,capture_output=True,check=True)
    length_match=re.findall(r"(?im)^content-length:\s*(\d+)\s*$",head.stdout)
    ranges="accept-ranges: bytes" in head.stdout.lower()
    if not length_match or not ranges:raise RuntimeError("Server must support exact range requests")
    total=int(length_match[-1])
    raw=args.dest/f"{args.accession}_buildGRCh37.tsv.gz"
    partial=Path(str(raw)+".part")
    if raw.exists():raise FileExistsError("Verified target exists; no overwrite")
    if not partial.is_file():raise FileNotFoundError("Expected interrupted .part prefix")
    start=partial.stat().st_size
    if not 0<start<total:raise ValueError("Invalid prefix byte length")
    plan={"accession":args.accession,"total_bytes":total,"existing_prefix_bytes":start,
          "remaining_bytes":total-start,"parts":args.parts,"status":"PLAN_ONLY"}
    if not args.execute:
        print(json.dumps(plan,indent=2));return
    folder=Path(str(raw)+".ranges")
    folder.mkdir(parents=True,exist_ok=True)
    pending=total-start
    spans=[]
    for i in range(args.parts):
        lo=start+(pending*i)//args.parts
        hi=start+(pending*(i+1))//args.parts-1
        if hi>=lo:spans.append((i,lo,hi))
    def fetch(item):
        idx,lo,hi=item
        dest=folder/f"{idx:02}.bin"
        header=folder/f"{idx:02}.headers"
        expected=hi-lo+1
        if dest.is_file() and dest.stat().st_size==expected and header.is_file():
            meta=header.read_text().lower()
            if f"content-range: bytes {lo}-{hi}/{total}" in meta:
                return item
        cmd=["curl","-fsSL","--retry","3","--connect-timeout","20","--speed-time","70",
             "--speed-limit","1024","--range",f"{lo}-{hi}",
             "-D",str(header),"-o",str(dest),url]
        subprocess.run(cmd,check=True,timeout=args.timeout)
        hdr=header.read_text().lower()
        if f"content-range: bytes {lo}-{hi}/{total}" not in hdr or dest.stat().st_size!=expected:
            raise ValueError(f"HTTP range integrity failed for {idx}")
        return item
    with ThreadPoolExecutor(max_workers=args.parts) as executor:
        futures=[executor.submit(fetch,item) for item in spans]
        for future in as_completed(futures):
            print("segment_done",future.result()[0],flush=True)
    stitched=Path(str(raw)+".assembled")
    md5=hashlib.md5()
    with stitched.open("wb") as out:
        for src in [partial]+[folder/f"{i:02}.bin" for i,_,_ in spans]:
            with src.open("rb") as f:
                for buf in iter(lambda:f.read(1024*1024),b""):
                    out.write(buf);md5.update(buf)
    result=md5.hexdigest()
    if result!=r["source_md5"] or stitched.stat().st_size!=total:
        raise ValueError("MD5 or length mismatch: source intact; no promotion")
    stitched.replace(raw)
    # Temporary inputs are safe to remove only after exact-length MD5 verification.
    partial.unlink()
    for i,_,_ in spans:
        (folder/f"{i:02}.bin").unlink()
        (folder/f"{i:02}.headers").unlink()
    folder.rmdir()
    plan.update({"status":"MD5_VERIFIED_FINAL","md5":result,"final":str(raw)})
    print(json.dumps(plan,indent=2))
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--accession",required=True)
    p.add_argument("--manifest",type=Path,default=MANIFEST)
    p.add_argument("--dest",type=Path,default=DEST)
    p.add_argument("--parts",type=int,choices=[2,3,4,5,6],default=4)
    p.add_argument("--timeout",type=int,default=1800)
    p.add_argument("--execute",action="store_true")
    acquire(p.parse_args())
