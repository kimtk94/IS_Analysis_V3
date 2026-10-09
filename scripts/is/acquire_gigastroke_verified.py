#!/usr/bin/env python3
"""Checksum-gated GIGASTROKE acquisition, safe by default (plan only).

Explicit --execute and exact --accession required. Never overwrite canonical files.
"""
import argparse,csv,hashlib,json,subprocess
from pathlib import Path
PLAN=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/GIGASTROKE_SOURCE_VERIFIED_V3.tsv")
DEST=Path("/srv/is-analysis/data/is/reference/gigastroke/additional")
def digest(path):
    md5=hashlib.md5()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(2**20),b""):md5.update(block)
    return md5.hexdigest()
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",type=Path,default=PLAN)
    p.add_argument("--dest",type=Path,default=DEST)
    p.add_argument("--accession",required=True)
    p.add_argument("--execute",action="store_true")
    args=p.parse_args()
    with args.manifest.open(newline="") as f:
        selected=[r for r in csv.DictReader(f,delimiter="\t") if r["accession"]==args.accession]
    if len(selected)!=1:raise ValueError("Not exactly one accession found")
    r=selected[0]
    if not r["accession"].startswith("GCST") or len(r["source_md5"])!=32:
        raise ValueError("Invalid source")
    url=r["source_url"]
    if not url.startswith("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"):
        raise ValueError("Unapproved source URL")
    output=args.dest/f'{r["accession"]}_buildGRCh37.tsv.gz'
    existing=output.is_file() and digest(output)==r["source_md5"]
    preview={"accession":r["accession"],"phenotype":r["phenotype"],
        "ancestry":r["ancestry_verified"],"url":url,"destination":str(output),
        "expected_md5":r["source_md5"],"verified_already":existing,
        "status":"ALREADY_VERIFIED" if existing else "PLAN_ONLY"}
    if not args.execute or existing:
        print(json.dumps(preview,indent=2));return
    args.dest.mkdir(parents=True,exist_ok=True)
    partial=Path(str(output)+".part")
    subprocess.run(["curl","-fL","--retry","3","--connect-timeout","20",
        "--speed-time","90","--speed-limit","1024","-C","-","-o",str(partial),url],
        check=True)
    actual=digest(partial)
    if actual!=r["source_md5"]:
        raise ValueError("Downloaded checksum mismatch; incomplete file retained as .part")
    partial.replace(output)
    preview.update({"status":"DOWNLOADED_MD5_VERIFIED","actual_md5":actual})
    print(json.dumps(preview,indent=2))
if __name__=="__main__":main()
