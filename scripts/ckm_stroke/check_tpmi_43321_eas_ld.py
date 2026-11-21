#!/usr/bin/env python3
"""TPMI 433.21 rare lead: fail-closed GRCh37 Ensembl 1000G EAS LD evidence scan.

Never equate a neighboring SNP, high D', or a strong LD proxy with outcome
replication. This script does NOT claim that any proxy appears in GIGASTROKE.
"""
import argparse
import json
import math
import os
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

QUERY="rs187759247"
BUILD="GRCh37"
TARGET="19:15573379:T:C"
# Predefined 1000 Genomes phase 3 East Asian reference subpopulations.
POPS=("CHB","CHS","JPT","CDX","KHV")
API="https://grch37.rest.ensembl.org"
MIN_R2=0.80

def validate(items, population):
    if not isinstance(items,list):
        raise ValueError("non-list Ensembl LD response")
    kept=[]
    for item in items:
        if not isinstance(item,dict):
            raise ValueError("non-dict LD row")
        a,b=item.get("variation1"),item.get("variation2")
        if QUERY not in (a,b):
            raise ValueError("response does not involve target rsID")
        other=b if a==QUERY else a
        r2=item.get("r2")
        try:
            r2=float(r2)
        except (ValueError,TypeError):
            raise ValueError("missing/invalid LD r2")
        if not math.isfinite(r2) or r2<0 or r2>1:
            raise ValueError("out-of-range LD r2")
        if other==QUERY or not isinstance(other,str) or not other.startswith("rs"):
            continue
        if r2>=MIN_R2:
            kept.append({"proxy_rsid":other,"population":population,"r2":r2,
                         "d_prime":item.get("d_prime"),
                         "phase_known":False,"outcome_alleles_matched":False})
    return kept

def fetch(pop):
    code="1000GENOMES:phase_3:"+pop
    url=f"{API}/ld/human/{QUERY}/{urllib.parse.quote(code,safe=':')}?r2={MIN_R2}"
    req=urllib.request.Request(url,headers={"Accept":"application/json",
                           "User-Agent":"TPMI-43321-EAS-LD-audit/1.0"})
    with urllib.request.urlopen(req,timeout=35) as response:
        if response.status!=200 or response.geturl().split("/")[2]!="grch37.rest.ensembl.org":
            raise ValueError("LD endpoint redirected or changed")
        data=response.read(2000001)
    if len(data)>2000000:raise ValueError("Ensembl LD response exceeds cap")
    return validate(json.loads(data.decode("utf-8")),pop)

def run():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    rows=[];states={}
    for pop in POPS:
        try:
            result=fetch(pop)
            rows.extend(result)
            states[pop]={"status":"PASS_SOURCE_QUERIED","n_high_ld":len(result)}
        except Exception as exc:
            states[pop]={"status":"HOLD_REFERENCE_UNAVAILABLE",
                         "error_type":type(exc).__name__,
                         "error_message":str(exc)[:240]}
    ranking={}
    for row in rows:
        obj=ranking.setdefault(row["proxy_rsid"],{"rsid":row["proxy_rsid"],
                    "populations":[],"r2_max":0.0,"r2_min":1.0,
                    "phase_known":False,"gigastroke_variant_matched":False})
        obj["populations"].append(row["population"])
        obj["r2_max"]=max(obj["r2_max"],row["r2"])
        obj["r2_min"]=min(obj["r2_min"],row["r2"])
    candidates=sorted(ranking.values(),key=lambda o:(-len(o["populations"]),-o["r2_max"],o["rsid"]))
    out={"schema":"ckm-tpmi-43321-eas-ld-proxy-gate-v1",
         "query":QUERY,"source_build":BUILD,"target":TARGET,
         "reference":"Ensembl GRCh37 REST / 1000 Genomes phase 3",
         "minimum_r2":MIN_R2,"populations":states,
         "high_ld_reference_pairs":len(rows),"distinct_high_ld_rsids":len(candidates),
         "candidates":candidates[:200],
         "candidates_truncated":len(candidates)>200,
         "reference_evidence_status":("REFERENCE_PROXIES_FOUND_NOT_HARMONIZED" if candidates else
           "HOLD_NO_PROVEN_PROXY"),
         "scientific_status":"HOLD",
         "note":"LD proxy candidates are not GIGASTROKE outcome matches; EAS phase and outcome locus/alleles, reference frequencies, and independence must be confirmed."}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":out["reference_evidence_status"],"populations":states,
                      "candidates":len(candidates),"science":"HOLD"}))
if __name__=="__main__":
    run()
