#!/usr/bin/env python3
"""P5e selection/collider-risk diagnostic from public GSE313736 cell metadata.

No individual donor counts/IDs in output. This is descriptive and not proof that
selected BeamA cells are identical in underlying state between DEX and Vehicle.
"""
import argparse
import csv
import gzip
import json
import math
import statistics
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from ocular_p1_integrity import BASE,download

def metric(values):
    v=sorted(float(x) for x in values)
    if not v:return {"n":0}
    return {"n":len(v),"min":round(v[0],5),
       "median":round(statistics.median(v),5),
       "max":round(v[-1],5)}
def sign_test(d):
    npos=sum(x>1e-12 for x in d)
    nneg=sum(x< -1e-12 for x in d)
    n=npos+nneg
    if n==0:return {"n":0,"positive":0,"negative":0,"two_sided_p":1.}
    lo=min(npos,nneg)
    p=min(1.,2*sum(math.comb(n,i) for i in range(lo+1))/(2**n))
    return {"n":n,"positive":npos,"negative":nneg,"two_sided_p":round(p,8)}
def evaluate(fp):
    stats=defaultdict(Counter)
    conf=defaultdict(list)
    by_donor=defaultdict(Counter)
    pools=defaultdict(Counter)
    cross=Counter()
    n=0
    with gzip.open(fp,"rt",encoding="utf8") as handle:
        rd=csv.DictReader(handle,delimiter="\t")
        required={"barcode","DEX Treatment","donor_id_demuxlet","droplet_type",
            "celltypist_chen2025_fine","celltypist_chen2025_fine_score",
            "scpred_vanzyl2022","library_id"}
        if not required.issubset(rd.fieldnames or []):
            raise ValueError("CELL_ANNOTATION_COLUMNS_MISSING")
        for r in rd:
            n+=1
            if r["droplet_type"]!="SNG":continue
            a=r["DEX Treatment"].strip().lower()
            d=r["donor_id_demuxlet"].strip()
            if a not in ("dex","no-dex") or not d:
                raise ValueError("UNEXPECTED_CONDITION_OR_DONOR")
            label=r["celltypist_chen2025_fine"].strip()
            sc=r["scpred_vanzyl2022"].strip()
            try:v=float(r["celltypist_chen2025_fine_score"])
            except (TypeError,ValueError):v=-1.
            full=(d,a)
            stats[a]["all"]+=1
            by_donor[full]["all"]+=1
            is_beam=label=="BeamA"
            if is_beam:
                stats[a]["BeamA_any"]+=1
                by_donor[full]["BeamA_any"]+=1
                if v>=.5:
                    stats[a]["BeamA_0.5"]+=1
                    by_donor[full]["BeamA_0.5"]+=1
                if v>=.75:
                    stats[a]["BeamA_0.75"]+=1
                    by_donor[full]["BeamA_0.75"]+=1
            if label=="JCT":
                stats[a]["JCT_any"]+=1
            if sc=="TM_Fibro":
                stats[a]["scpred_TM_Fibro_any"]+=1
            if is_beam and sc=="TM_Fibro":
                cross[(a,"BeamA_and_TM_Fibro")]+=1
            if v>=0:
                conf[(a,"all")].append(v)
                if is_beam:conf[(a,"BeamA")].append(v)
            pool=r["library_id"].split("pool-")[-1].split("_")[0].replace("dex","")
            if len(pool)==1 and pool in "ABCDEFGHI":
                pools[(pool,a)]["all"]+=1
                if is_beam and v>=.75:pools[(pool,a)]["BeamA_0.75"]+=1
    donors=sorted(set(k[0] for k in by_donor if k[1]=="dex")&
                  set(k[0] for k in by_donor if k[1]=="no-dex"))
    if n!=169398 or len(donors)!=59:raise ValueError("DATASET_SINGLET_OR_DONOR_COUNTS_CHANGED")
    arm={}
    for a in ("no-dex","dex"):
        t=stats[a]
        arm[a]={"total_sng":t["all"],"beamA_any":t["BeamA_any"],
                "beamA_score_ge_0.5":t["BeamA_0.5"],
                "beamA_score_ge_0.75":t["BeamA_0.75"],
                "jct_any":t["JCT_any"],
                "scpred_TM_Fibro":t["scpred_TM_Fibro_any"],
                "beamA_fraction_ge_0.75":round(t["BeamA_0.75"]/t["all"],6),
                "beamA_score_ge_0.75_fraction_of_beamA":round(t["BeamA_0.75"]/max(1,t["BeamA_any"]),6),
                "celltypist_score_all":metric(conf[(a,"all")]),
                "celltypist_score_beamA":metric(conf[(a,"BeamA")]),
                "beamA_and_scpred_TM_Fibro":cross[(a,"BeamA_and_TM_Fibro")]}
    diffs={}
    for kind in ("BeamA_any","BeamA_0.5","BeamA_0.75"):
        ds=[(by_donor[d,"dex"][kind]/by_donor[d,"dex"]["all"])-
            (by_donor[d,"no-dex"][kind]/by_donor[d,"no-dex"]["all"])
            for d in donors if by_donor[d,"dex"]["all"] and by_donor[d,"no-dex"]["all"]]
        diffs[kind]={"paired_fraction_delta_dex_minus_vehicle":metric(ds),
                     "paired_sign_test":sign_test(ds)}
    bypool=[]
    for p in "ABCDEFGHI":
        item={"pool":p}
        for a in ("no-dex","dex"):
            t=pools[p,a]
            item[a]={"total_sng":t["all"],
                     "beamA_frac_ge_0_75":round(t["BeamA_0.75"]/max(1,t["all"]),6)}
        bypool.append(item)
    return {"dataset":"GSE313736","run_utc":datetime.now(timezone.utc).isoformat(),
       "status":"PASS_CELL_STATE_SELECTION_DIAGNOSTIC",
       "biological_donors_with_both_arms":len(donors),
       "cells":n,"arms":arm,"donor_paired_differences":diffs,
       "per_pool_aggregate":bypool,
       "interpretation":["DEX-associated changes in annotation frequency may represent biology AND/OR classifier re-labeling",
           "Conditioning on high-confidence BeamA can create selection bias; these diagnostics do not establish a treatment-causal effect",
           "Counts represent cells, not independent donor observations for inference",
           "Only pooled summaries and paired sign-test counts, no donor identifiers saved"]}
def main():
    a=argparse.ArgumentParser()
    a.add_argument("--output",required=True)
    opts=a.parse_args()
    evidence={"status":"HOLD"}
    with tempfile.TemporaryDirectory() as td:
        try:
            p=Path(td)/"meta.tsv.gz"
            download(BASE+"series/GSE313nnn/GSE313736/suppl/GSE313736_HTMC_scRNA_cell_metadata.tsv.gz",
                     p,12*1024*1024)
            evidence=evaluate(p)
        except Exception as e:
            evidence["error"]=type(e).__name__+": "+str(e)[:200]
    Path(opts.output).write_text(json.dumps(evidence,indent=2,ensure_ascii=False)+"\n")
    print("P5E_SELECTION",evidence["status"],"arm_fraction",
          [(x,y.get("beamA_fraction_ge_0.75")) for x,y in evidence.get("arms",{}).items()],
          "error",evidence.get("error"))
    if evidence["status"]!="PASS_CELL_STATE_SELECTION_DIAGNOSTIC":raise SystemExit(2)
if __name__=="__main__":main()
