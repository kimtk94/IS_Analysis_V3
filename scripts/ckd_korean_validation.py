#!/usr/bin/env python3
"""Fixed Korean outcome follow-up of the previous public CKB exploratory screen.

No CKB encrypted statistics or private individual-level data are accessed.
"""
import argparse
import collections
import gzip
import hashlib
import io
import json
import math
import re
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

PIN="ad0d9181d3d1b4d598b7c9e5a3d5be11825db8681a877069a4410c96834328c5"
EXPECTED={"rs16846391","rs2070600","rs4390169","rs17050272","rs77924615","rs10185736","rs11162405"}
RECORD="16976585"
COMP=str.maketrans("ACGT","TGCA")
PHEWEB="https://koges.leelabsg.org"

def digest(p,alg="sha256"):
    with p.open("rb") as f: return hashlib.file_digest(f,alg).hexdigest()

def get(url,cap=2000000):
    with urllib.request.urlopen(urllib.request.Request(url,headers={"User-Agent":"CKD-public-summary-research/1.0","Accept":"application/json, */*"}),timeout=25) as f:
        raw=f.read(cap+1)
    if len(raw)>cap: raise ValueError("Response exceeds bound")
    return raw

def anchors(root):
    p=root/"egfr.zip"
    if digest(p)!=PIN: raise ValueError("Prior eGFR artifact digest changed")
    with zipfile.ZipFile(p) as z:
        if z.testzip(): raise ValueError("Prior archive CRC")
        name=[n for n in z.namelist() if n.rsplit("/",1)[-1]=="assay_results.json"]
        if len(name)!=1: raise ValueError("Prior assay table missing")
        rows=json.loads(z.read(name[0]))
    out=collections.defaultdict(list)
    for r in rows:
        if r["status"]=="EXPLORATORY_WALD" and r["bonferroni_p"]<.05:
            out[r["rsid"]].append(r)
    if set(out)!=EXPECTED or sum(map(len,out.values()))!=10: raise ValueError("Frozen family changed")
    for rs,rr in out.items():
        a=rr[0]
        if any((x["chr37"],x["pos37"],x["ref37"],x["alt37"])!=(a["chr37"],a["pos37"],a["ref37"],a["alt37"]) for x in rr):
            raise ValueError("Anchor discordance")
    return dict(sorted(out.items()))

def orient(ea,oa,a,b):
    if {ea,oa} in ({"A","T"},{"C","G"}):return None
    if (ea,oa)==(a,b):return 1
    if (ea,oa)==(b,a):return -1
    return None

def wald(bx,sx,by,sy):
    if not all(math.isfinite(x) for x in (bx,sx,by,sy)) or bx==0 or min(sx,sy)<=0:raise ValueError("Nonfinite estimate")
    b=by/bx
    se=math.hypot(sy,b*sx)/abs(bx)
    return {"mr_beta":b,"mr_se":se,"mr_ci95":[b-1.959964*se,b+1.959964*se],
            "mr_p":math.erfc(abs(b/se)/math.sqrt(2))}

def summarize(rs,rr,by,sy,ea,oa,source,expected_egfr=True):
    a=rr[0]
    s=orient(a["ea37"],a["oa37"],ea,oa)
    if s is None:raise ValueError("Allele mismatch or palindrome")
    p=math.erfc(abs(by/sy)/math.sqrt(2))
    result={"rsid":rs,"outcome_source":source,"outcome_ea":ea,"outcome_oa":oa,
            "outcome_beta":by,"outcome_se":sy,"two_sided_p":p,
            "bonferroni_p_7":min(1,7*p),"egfr_direction_agrees":bool(by*s*a["outcome_beta_aligned"]>0) if expected_egfr else None,
            "assays":[]}
    for r in rr:
        sign=orient(r["ea37"],r["oa37"],ea,oa)
        if sign is None:raise ValueError("Assay allele mismatch")
        m=wald(r["beta"],r["se"],by*sign,sy)
        result["assays"].append({"platform":r["platform"],"assay_id":r["assay_id"],
            "exposure_beta":r["beta"],"exposure_se":r["se"],"previous_egfr_mr_beta":r["mr_beta"],
            "new_outcome_beta_aligned":by*sign,**m})
    return result

def write(out,obj):
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(obj,indent=2,allow_nan=False)+"\n")
    print(json.dumps(obj,indent=2,allow_nan=False))

def pheweb(root,out):
    panel=anchors(root)
    result=[]
    for rs,rr in panel.items():
        a=rr[0]
        variant=f'{a["chr37"]}-{a["pos37"]}-{a["ref37"]}-{a["alt37"]}'
        url=f"{PHEWEB}/api/variant/{variant}"
        record={"rsid":rs,"variant":variant,"url":url,"status":"HOLD"}
        try:
            raw=get(url);p=json.loads(raw)
            record["sha256"]=hashlib.sha256(raw).hexdigest()
            if not isinstance(p,dict):
                raise ValueError("No variant record")
            record["header"]=sorted(p)
            ref,alt=str(p.get("ref","")).upper(),str(p.get("alt","")).upper()
            if (str(p.get("chrom")),int(p.get("pos",0)),ref,alt)!=(str(a["chr37"]),a["pos37"],a["ref37"],a["alt37"]):
                raise ValueError("PheWeb variant identity not exact")
            reported=str(p.get("rsids",p.get("rsid","")))
            if rs not in re.findall(r"rs\d+",reported):
                raise ValueError("PheWeb rsID not verified")
            phenos=p.get("phenos",[])
            if not isinstance(phenos,list):raise ValueError("Unexpected phenos")
            record["available_renal_labels"]=[{"label":x.get("phenostring"),"code":x.get("phenocode"),"fields":sorted(x)} for x in phenos if re.search(r"creatinine|urea|kidney|glomerular|egfr|renal",str(x.get("phenostring","")),re.I)]
            selected=[]
            for label in ("Blood Urea Nitrogen","Creatinine","Estimated Glomerular Filtration Rate","eGFR"):
                hits=[x for x in phenos if str(x.get("phenostring","")).casefold()==label.casefold()]
                if len(hits)!=1:continue
                x=hits[0]
                by,sy,pval=map(float,(x["beta"],x["sebeta"],x["pval"]))
                if not all(math.isfinite(t) for t in (by,sy,pval)) or sy<=0 or not 0<=pval<=1:raise ValueError("Invalid phenotype association")
                # Site explicitly defines BETA relative to ALT. Do not use BBJ meta-analysis.
                estimate=summarize(rs,rr,by,sy,alt,ref,label,expected_egfr=label.lower() in ("egfr","estimated glomerular filtration rate"))
                selected.append({"label":label,"phenocode":x.get("phenocode"),"source_p":pval,"result":estimate})
            record.update(status="EXACT_PUBLIC_VARIANT_QUERY",supporting_traits=selected)
        except (OSError,ValueError,KeyError,TypeError,OverflowError) as e:
            record["reason"]=str(e)[:300]
        result.append(record)
    write(out,{"status":"PUBLIC_KOGES_2022_VARIANT_QUERY","primary_family":sorted(panel),
               "records":result,"limits":"KoGES 2022 is supportive; BUN/creatinine are not eGFR. Do not treat overlapping KoGES releases as independent."})

def zenodo_meta(root,out):
    raw=get("https://zenodo.org/api/records/"+RECORD,5000000)
    p=json.loads(raw)
    meta=p["metadata"]
    if str(p["id"])!=RECORD or meta.get("doi")!="10.5281/zenodo."+RECORD:raise ValueError("Record identity mismatch")
    files=[{"key":x["key"],"size":x["size"],"checksum":x.get("checksum"),"download":x.get("links",{}).get("self")} for x in p["files"]]
    result={"status":"PUBLIC_RECORD_METADATA","record_id":RECORD,"doi":meta.get("doi"),"title":meta.get("title"),
            "description":str(meta.get("description",""))[:4000],"access_right":meta.get("access_right"),
            "license":meta.get("license"),"files":files,"sha256":hashlib.sha256(raw).hexdigest()}
    write(out,result)
    # This is an official record, never use a file that is not clearly a public GWAS summary statistic.
    root.mkdir(parents=True,exist_ok=True)
    (root/"zenodo_record.json").write_bytes(raw)

def selftest():
    assert orient("G","A","A","G")==-1
    assert orient("A","T","A","T") is None
    assert abs(wald(.5,.02,.1,.01)["mr_beta"]-.2)<1e-12
    print("Synthetic tests passed")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("mode",choices=["selftest","pheweb","zenodo"])
    p.add_argument("--root",type=Path,default=Path("."))
    p.add_argument("--out",type=Path,default=Path("output.json"))
    a=p.parse_args()
    if a.mode=="selftest":selftest()
    elif a.mode=="pheweb":pheweb(a.root,a.out)
    else:zenodo_meta(a.root,a.out)
if __name__=="__main__":main()
