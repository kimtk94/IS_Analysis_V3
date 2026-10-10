#!/usr/bin/env python3
"""Guarded official Zenodo Japanese alcohol GWAS download.

DEFAULT PLAN ONLY. Use --execute --file {alcohol,drinking}, retries source
transfers but NEVER calls a partial/size-only file validated. Stores checksum
and source metadata in untracked data workspace, not git.
"""
import argparse,hashlib,json,subprocess
from pathlib import Path
ROOT=Path("/srv/is-analysis/data/is/gwas_exposure/japanese_alcohol_2024")
FILES={"alcohol":"1_Alcohol_intake_Unstratified.tsv.gz",
       "drinking":"5_Drinking_Unstratified.tsv.gz"}
def planned(meta,root,requested):
    by={x["key"]:x for x in meta["files"]}
    plan=[]
    for item in requested:
        name=FILES[item];data=by[name]
        if not data["checksum"].startswith("md5:") or len(data["checksum"])!=36:
            raise ValueError("Missing official MD5")
        url=data["links"]["self"]
        if not url.startswith("https://zenodo.org/api/records/10038152/files/"):
            raise ValueError("Non-official URL blocked")
        plan.append({"label":item,"name":name,"size":data["size"],
          "md5":data["checksum"][4:],"url":url,
          "complete_path":str(root/name),"part_path":str(root/(name+".part")),
          "status":"PLAN_ONLY_NOT_VERIFIED"})
    return plan
def checksum_file(path):
    h=hashlib.md5()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()
def execute(plan):
    for rec in plan:
        out=Path(rec["complete_path"]);partial=Path(rec["part_path"])
        if not (out.exists() and out.stat().st_size==rec["size"] and checksum_file(out)==rec["md5"]):
            for attempt in range(1,5):
                cmd=["curl","--fail","--location","--show-error","--silent",
                     "--retry","2","--connect-timeout","20",
                     "--speed-time","30","--speed-limit","1024",
                     "--max-time","150","--continue-at","-",
                     "--output",str(partial),rec["url"]]
                response=subprocess.run(cmd)
                if response.returncode==0 and partial.is_file() and partial.stat().st_size==rec["size"]:
                    if checksum_file(partial)!=rec["md5"]:
                        raise ValueError("Checksum mismatch after full download")
                    partial.replace(out)
                    break
                print("TRANSFER_RETRY",rec["name"],attempt,
                      "partial_bytes",partial.stat().st_size if partial.exists() else 0,flush=True)
            else:
                raise RuntimeError("SOURCE_TRANSFER_INCOMPLETE "+rec["name"])
        if not(out.is_file() and out.stat().st_size==rec["size"] and checksum_file(out)==rec["md5"]):
            raise ValueError("Source file is NOT MD5 valid "+rec["name"])
        rec["status"]="SOURCE_MD5_VERIFIED"
        print("SOURCE_MD5_VERIFIED",rec["name"],out.stat().st_size,flush=True)
    return plan
def run(args):
    if not args.meta.exists():raise FileNotFoundError("Official Zenodo metadata not cached")
    m=json.loads(args.meta.read_text())
    if str(m.get("id"))!="10038152":raise ValueError("Incorrect record ID")
    args.root.mkdir(parents=True,exist_ok=True)
    labels=[args.file] if args.file!="both" else list(FILES)
    plan=planned(m,args.root,labels)
    if args.execute:execute(plan)
    out=args.root/"G0022_ALCOHOL_PUBLIC_GWAS_DOWNLOAD_STATUS.json"
    out.write_text(json.dumps({"record":"zenodo_10038152","executed":args.execute,
        "source_files":plan},indent=2))
    print(json.dumps(plan,indent=2))
    return plan
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=ROOT)
    p.add_argument("--meta",type=Path,default=ROOT/"ZENODO_10038152_METADATA.json")
    p.add_argument("--file",choices=(*FILES,"both"),default="both")
    p.add_argument("--execute",action="store_true")
    run(p.parse_args())
