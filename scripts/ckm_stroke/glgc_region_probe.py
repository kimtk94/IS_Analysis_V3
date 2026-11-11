#!/usr/bin/env python3
"""Probe three official GLGC EUR 2021 indexed summary-statistics regions.

Uses the consortium's .tbi index via htslib/tabix HTTP range access.
THIS IS A PARTIAL REGION/QC DEMONSTRATION, not whole-genome completion.
"""
import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import math
import os
import re
import subprocess
from pathlib import Path

BASE="https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/"
FILE_NAMES={
    "LDL_C":"LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz",
    "HDL_C":"HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz",
    "TG":"logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz",
}
REQUIRED={"rsID","CHROM","POS_b37","REF","ALT","N","EFFECT_SIZE","SE","pvalue"}
AUTOSOMES={str(x) for x in range(1,23)}
ALLOWED={"A","C","G","T"}

def validate_rows(text,chrom,start,end):
    lines=[s for s in text.splitlines() if s.strip()]
    if not lines:
        return {"status":"HOLD","reason":"no output"}
    header=None
    for s in lines:
        if s.startswith("##"):
            continue
        if s.startswith("#"):
            s=s.lstrip("#")
        header=s
        break
    if header is None:
        return {"status":"HOLD","reason":"tabix returned no column header"}
    f=csv.DictReader(io.StringIO("\n".join(lines[lines.index(next(x for x in lines if not x.startswith("##"))):])),delimiter="\t")
    # Most GLGC files have a plain first header row, not a VCF metadata section.
    fields=set(f.fieldnames or [])
    if not REQUIRED.issubset(fields):
        return {"status":"FAIL","reason":"missing required GLGC fields",
                "columns":sorted(fields),"missing":sorted(REQUIRED-fields)}
    n=0; invalid=0; n_pal=0; many=0; examples=[]; p_invalid=0
    for row in f:
        n+=1
        try:
            chromval=row["CHROM"].removeprefix("chr")
            pos=int(row["POS_b37"])
            ref,alt=row["REF"].upper(),row["ALT"].upper()
            beta,se,p=map(float,(row["EFFECT_SIZE"],row["SE"],row["pvalue"]))
            if chromval!=chrom or not start<=pos<=end:
                raise ValueError("unexpected chromosome/region")
            if len(ref)!=1 or len(alt)!=1 or ref not in ALLOWED or alt not in ALLOWED or ref==alt:
                many+=1
                continue
            if not all(map(math.isfinite,(beta,se,p))) or se<=0 or not 0<p<=1:
                p_invalid+=1
                continue
            if frozenset((ref,alt)) in {frozenset(("A","T")),frozenset(("C","G"))}:
                n_pal+=1
            if len(examples)<2:
                examples.append({"chrom":chromval,"pos":pos,"ref":ref,"alt":alt,"beta":beta,"se":se,"p":p})
        except (ValueError,TypeError,KeyError,AttributeError):
            invalid+=1
    if n<=0: return {"status":"HOLD","reason":"region returned zero variants","header":sorted(fields)}
    return {"status":"HOLD" if invalid or p_invalid else "PASS",
            "rows_region":n,"invalid_rows":invalid,
            "non_snp_rows":many,"invalid_effect_rows":p_invalid,
            "palindromic_snp_rows":n_pal,
            "observed_fields":sorted(fields),"sample_variants":examples,
            "note":"Region-only; nothing has been deduplicated, allele-harmonized, or compared to another GWAS."}


def probe(trait,chrom="22",start=20000000,end=20050000,timeout=180):
    if trait not in FILE_NAMES: raise ValueError("unknown trait")
    if chrom not in AUTOSOMES or not 1<=start<end:
        raise ValueError("invalid coordinate")
    url=BASE+FILE_NAMES[trait]
    argv=["tabix","-h",url,f"{chrom}:{start}-{end}"]
    try:
        p=subprocess.run(argv,capture_output=True,text=True,timeout=timeout,check=False)
    except (OSError,subprocess.TimeoutExpired) as exc:
        return {"status":"SOURCE_UNAVAILABLE","reason":str(exc)[:300],"url":url}
    if p.returncode!=0:
        return {"status":"SOURCE_UNAVAILABLE","returncode":p.returncode,
                "stderr_tail":p.stderr[-1000:],"url":url,
                "advice":"Check remote .tbi index availability and HTTP byte ranges"}
    parsed=validate_rows(p.stdout,chrom,start,end)
    parsed.update({"url":url,"index_url":url+".tbi",
                   "source":"GLGC_Graham2021","trait":trait,"ancestry":"EUR",
                   "region":f"{chrom}:{start}-{end}",
                   "tabix_stdout_sha256":hashlib.sha256(p.stdout.encode()).hexdigest(),
                   "scientific_status":"NOT_RUN",
                   "scope":"INDEXED_REGION_ONLY_NOT_FULL_GWAS"})
    return parsed


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--trait",required=True,choices=sorted(FILE_NAMES))
    ap.add_argument("--chrom",default="22")
    ap.add_argument("--start",type=int,default=20000000)
    ap.add_argument("--end",type=int,default=20050000)
    ap.add_argument("--out",default="artifacts/ckm_stroke/glgc_region_probe.json")
    args=ap.parse_args()
    result=probe(args.trait,args.chrom,args.start,args.end)
    result.update({"source_commit":os.getenv("GITHUB_SHA",""),
                   "checked_at_utc":dt.datetime.now(dt.timezone.utc).isoformat()})
    out=Path(args.out)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False),flush=True)
    summary=os.getenv("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary,"a",encoding="utf-8") as handle:
            handle.write("\n### GLGC indexed EUR "+args.trait+"\n"+
                         "- Status: "+result["status"]+"\n"+
                         "- Region: "+result.get("region",f"{args.chrom}:{args.start}-{args.end}")+"\n"+
                         "- Records in region: "+str(result.get("rows_region",0))+"\n"+
                         "- Scientific status: NOT_RUN\n")
    return int(result["status"]=="FAIL")

if __name__=="__main__":
    raise SystemExit(main())
