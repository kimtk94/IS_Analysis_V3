#!/usr/bin/env python3
"""Bounded public RDA transfer, not a full-table Python RData conversion."""
import argparse,hashlib,json,subprocess
from pathlib import Path
BASE="https://raw.githubusercontent.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/70461ea4a55716bb00f8467860521e8a8f94b9f9/data/"
SOURCES={"MUSCLE_TRNSCRPT_DA.rda":(22711273,"0a4168c61102863487a4c4ecb95ec66bc250329e"),
"HUMAN_FEATURE_TO_GENE.rda":(20001757,"db4dd772a08ed10b1063201179386a2c113ecb55")}
def main():
    a=argparse.ArgumentParser();a.add_argument("--root",type=Path,required=True);a.add_argument("--out",type=Path,required=True);p=a.parse_args()
    p.root.mkdir(parents=True,exist_ok=True);p.out.parent.mkdir(parents=True,exist_ok=True)
    result=[]
    for f,(size,sha) in SOURCES.items():
        path=p.root/f;url=BASE+f
        print("Downloading",f,flush=True)
        subprocess.run(["curl","-fLsS","--connect-timeout","15","--max-time","90","--max-filesize",str(size),url,"-o",str(path)],check=True,timeout=100)
        data=path.read_bytes();got=hashlib.sha1(("blob "+str(len(data))).encode()+b"\0"+data).hexdigest()
        if len(data)!=size or got!=sha:raise ValueError("Publisher object identity "+f)
        result.append({"file":f,"bytes":size,"git_blob":got,"sha256":hashlib.sha256(data).hexdigest(),"url":url})
        print(json.dumps(result[-1]),flush=True)
    p.out.write_text(json.dumps(result,indent=2)+"\n")
if __name__=="__main__":main()
