#!/usr/bin/env python3
"""CKM-stroke research gates: technical execution never equals scientific validation."""
import argparse
import csv
import datetime
import gzip
import json
import math
import os
from pathlib import Path
import subprocess

COLS={"phenotype","role","domain","ancestry","accession","build","data_file"}
UNKNOWN={"","REVIEW_REQUIRED","TBD"}

def manifest(root):
    with (root/"config/ckm_stroke/phenotypes.tsv").open(encoding="utf-8",newline="") as fh:
        reader=csv.DictReader(fh,delimiter="\t")
        if not reader.fieldnames or not COLS.issubset(reader.fieldnames):
            raise ValueError("manifest missing columns")
        records=list(reader)
    if not records: raise ValueError("empty manifest")
    return records

def p0(root):
    rows=manifest(root)
    seen=set()
    errors=[]
    todo=[]
    exposure_domains=set()
    outcomes=set()
    for r in rows:
        key=(r["phenotype"],r["role"],r["ancestry"])
        if key in seen: errors.append("duplicate "+str(key))
        seen.add(key)
        if r["role"] not in ("exposure","outcome"): errors.append("invalid role "+str(key))
        if r["ancestry"] not in ("EUR","EAS"): errors.append("invalid ancestry "+str(key))
        if r["build"] not in ("GRCh37","GRCh38","REVIEW_REQUIRED"): errors.append("invalid build "+str(key))
        if r["role"]=="exposure":
            exposure_domains.add(r["domain"])
            if "STROKE" in r["phenotype"].upper() or "CKM_STAGE" in r["phenotype"].upper():
                errors.append("circular endpoint "+str(key))
        if r["role"]=="outcome":
            outcomes.add(r["phenotype"])
            if r["domain"]!="stroke": errors.append("outcome domain must be stroke")
        for c in ("accession","build","data_file"):
            if r[c].strip().upper() in UNKNOWN: todo.append(str(key)+":"+c)
    if not {"metabolic","kidney"}.issubset(exposure_domains): errors.append("missing domain")
    if not {"LAA","CE","SVD"}.issubset(outcomes): errors.append("missing subtype")
    return {"status":"FAIL" if errors else "PASS","scientific_status":"HOLD",
            "schema_rows":len(rows),"unresolved":todo,"errors":errors}

def gwas_qc(path,limit=0):
    opener=gzip.open if path.suffix==".gz" else open
    with opener(path,"rt",newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f,delimiter="\t")
        if not reader.fieldnames or not {"CHR","POS","EA","NEA","BETA","SE","P"}.issubset(reader.fieldnames):
            return {"status":"FAIL","error":"normalized fields CHR POS EA NEA BETA SE P required"}
        count=0
        errors=[]
        seen=set()
        partial=False
        for row in reader:
            if limit and count>=limit:
                partial=True
                break
            count+=1
            try:
                chrom=row["CHR"].removeprefix("chr")
                pos=int(row["POS"])
                ea=row["EA"].upper()
                nea=row["NEA"].upper()
                beta,se,p=map(float,(row["BETA"],row["SE"],row["P"]))
                if chrom not in {str(x) for x in range(1,23)} or pos<=0: raise ValueError("chr/pos")
                if len(ea)!=1 or len(nea)!=1 or ea not in "ACGT" or nea not in "ACGT" or ea==nea: raise ValueError("allele")
                if not all(map(math.isfinite,(beta,se,p))) or se<=0 or not 0<p<=1: raise ValueError("beta/se/p")
                key=(chrom,pos,ea,nea)
                if key in seen: raise ValueError("duplicate")
                seen.add(key)
            except (KeyError,TypeError,ValueError) as e:
                if len(errors)<20: errors.append("row "+str(count)+" "+str(e))
        return {"status":"FAIL" if errors or count==0 else "HOLD" if partial else "PASS",
                "rows":count,"partial":partial,"errors":errors}

def p1(root,max_rows):
    missing=[]
    outputs={}
    for r in manifest(root):
        label=r["phenotype"]+"_"+r["ancestry"]
        file=r["data_file"]
        if file.strip().upper() in UNKNOWN or not (root/file).is_file():
            missing.append(label)
        else:
            outputs[label]=gwas_qc(root/file,max_rows)
    broken=[name for name,x in outputs.items() if x["status"]=="FAIL"]
    incomplete=[name for name,x in outputs.items() if x["status"]=="HOLD"]
    return {"status":"FAIL" if broken else "HOLD" if missing or incomplete else "PASS",
            "scientific_status":"HOLD","files":outputs,"missing":missing,
            "broken":broken,"incomplete":incomplete,
            "warning":"external checksum, clinical ascertainment, genome build, sample overlap not approved"}

def p2_to_p5(stage,root,execute,timeout):
    jobs=json.loads((root/"config/ckm_stroke/tasks.json").read_text(encoding="utf-8"))
    task=jobs.get(stage,{})
    if not all(isinstance(task.get(x),list) and task[x] for x in ("argv","inputs","outputs")):
        return {"status":"HOLD","scientific_status":"NOT_RUN","reason":"no approved task command/input/output"}
    missing=[x for x in task["inputs"] if not (root/x).is_file()]
    if missing: return {"status":"HOLD","scientific_status":"NOT_RUN","missing":missing}
    if not execute: return {"status":"HOLD","scientific_status":"NOT_RUN","reason":"execution not enabled"}
    args=task["argv"]
    if not all(isinstance(x,str) for x in args):
        return {"status":"FAIL","scientific_status":"NOT_RUN","reason":"argv must be strings"}
    try:
        p=subprocess.run(args,cwd=root,timeout=timeout,capture_output=True,text=True,check=False)
    except (OSError,subprocess.TimeoutExpired) as e:
        return {"status":"FAIL","scientific_status":"NOT_RUN","reason":str(e)}
    ok=p.returncode==0 and all((root/x).is_file() and (root/x).stat().st_size>0 for x in task["outputs"])
    return {"status":"PASS" if ok else "FAIL","scientific_status":"REVIEW_REQUIRED",
            "returncode":p.returncode,"stdout_tail":p.stdout[-1000:],"stderr_tail":p.stderr[-1000:]}

def p6(root):
    f=root/"data/ckm_stroke/replication.csv"
    if not f.is_file(): return {"status":"HOLD","scientific_status":"NOT_RUN","missing":[str(f)]}
    with f.open(encoding="utf-8",newline="") as fh:
        reader=csv.DictReader(fh)
        fields={"locus","beta_eur","se_eur","beta_eas","se_eas","harmonized","independent"}
        if not reader.fieldnames or not fields.issubset(reader.fieldnames):
            return {"status":"FAIL","scientific_status":"NOT_RUN","reason":"missing replication columns"}
        results=[]
        try:
            for r in reader:
                if r["harmonized"].lower()!="yes" or r["independent"].lower()!="yes":
                    return {"status":"HOLD","scientific_status":"NOT_RUN","reason":"harmonization or independence not established"}
                a,sa,b,sb=map(float,(r["beta_eur"],r["se_eur"],r["beta_eas"],r["se_eas"]))
                if not all(map(math.isfinite,(a,sa,b,sb))) or sa<=0 or sb<=0: raise ValueError("invalid estimate")
                z=(a-b)/math.hypot(sa,sb)
                results.append({"locus":r["locus"],"same_direction":a*b>0,
                                "difference_p":math.erfc(abs(z)/math.sqrt(2))})
        except ValueError as e:
            return {"status":"FAIL","scientific_status":"NOT_RUN","reason":str(e)}
    return {"status":"PASS" if results else "FAIL","scientific_status":"REVIEW_REQUIRED","results":results}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--stage",required=True,choices=["p0","p1","p2","p3","p4","p5","p6"])
    parser.add_argument("--root",default=".")
    parser.add_argument("--max-rows",type=int,default=0)
    parser.add_argument("--execute",action="store_true")
    parser.add_argument("--strict",action="store_true")
    parser.add_argument("--out",default="")
    a=parser.parse_args()
    root=Path(a.root).resolve()
    try:
        if a.stage=="p0": result=p0(root)
        elif a.stage=="p1": result=p1(root,a.max_rows)
        elif a.stage=="p6": result=p6(root)
        else: result=p2_to_p5(a.stage,root,a.execute,1500)
    except (OSError,ValueError,json.JSONDecodeError) as e:
        result={"status":"FAIL","scientific_status":"NOT_RUN","error":str(e)}
    result.update({"stage":a.stage,"run_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   "source_commit":os.getenv("GITHUB_SHA","")})
    path=Path(a.out or "artifacts/ckm_stroke/"+a.stage+".json")
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False))
    return 1 if result["status"]=="FAIL" or a.strict and result["status"]!="PASS" else 0

if __name__=="__main__":
    raise SystemExit(main())
