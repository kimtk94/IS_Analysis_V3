"""Seven-day whole-cohort score adapter. Separate Kallisto sensitivity, not DEE2 replay."""
import argparse
import csv
import hashlib
import json
import math
import statistics
import subprocess
from pathlib import Path
from ocular_fixed_panel_score_preflight import read_panel, evaluate, exact_sign, holm, require

PAIR_BLOB = "20b51fb9cf482d98003160c7875b07fa0666b1db"
MAP_SHA = "cb7305d9d6c0d3f336c388554a801513c4004d163cc7035ea448a3476fdc111e"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def median_interval(values):
    # Exact distribution-free order statistic interval under independent continuous observations.
    # If ties occur, coverage is conservative; if n is too small, report unbounded.
    n=len(values); ordered=sorted(values)
    candidates=[(k,1-2*sum(math.comb(n,i) for i in range(k))/2**n) for k in range(1,(n+1)//2+1)]
    valid=[(k,c) for k,c in candidates if c>=.95]
    if not valid: return {"lower":None,"upper":None,"coverage":1.,"note":"unbounded at >=95%"}
    k,c=valid[-1]
    return {"lower":ordered[k-1],"upper":ordered[-k],"coverage":c}

def fixed_pairs(path):
    require(subprocess.check_output(["git","hash-object",str(path)],text=True).strip()==PAIR_BLOB,"pair archive mismatch")
    p=json.loads(path.read_text())["pairs"]
    runs=[x[k] for x in p for k in ("ETH","DEX")]
    require(len(p)==8 and len(set(runs))==16,"complete eight distinct pairs required")
    require(sum(x["group"]=="GC-Responder" for x in p)==4 and sum(x["group"]=="GC-nonResponder" for x in p)==4,"group inventory")
    return p

def load_sample(root,run):
    d=root/("ocular-seven-day-gene-"+run)
    q=json.loads((d/"qc.json").read_text())
    s=json.loads((d/"summary.json").read_text())
    require(q["run"]==run and q["status"]=="TECHNICAL_COMPLETE_QC_REVIEW_REQUIRED" and q["gene_summary"]==s,"technical QC provenance")
    require(s["tx2gene_sha256"]==MAP_SHA and s["genes"]==62700,"reference identity")
    require("uniformly excluded" in q["orphan_policy"],"orphan policy")
    table=d/"gene_tpm.tsv"
    require(sha(table)==s["gene_table_sha256"],"gene output checksum")
    with table.open(newline="") as f:
        rd=csv.DictReader(f,delimiter="\t")
        require(rd.fieldnames==["gene_id","tpm","estimated_counts"],"gene TPM schema")
        values={}
        for row in rd:
            g=row["gene_id"]; v=float(row["tpm"])
            require(g not in values and math.isfinite(v) and v>=0,"invalid gene value")
            values[g]=v
    require(len(values)==62700 and abs(math.fsum(values.values())-1e6)<=100,"TPM universe or sum")
    return values

def log_tpm_sensitivity(panel,pairs):
    universe=set(pairs[0][0])
    require(all(set(x)==universe for p in pairs for x in p),"different universe")
    result={}
    for name,ids in (("all241",list(panel)),("ECM16",[g for g,(_,e) in panel.items() if e])):
        mapped=[g for g in ids if g in universe]
        require(len(mapped)>=math.ceil(.75*len(ids)),"overall coverage")
        for sign in (-1,1):
            require(sum(panel[g][0]==sign for g in mapped)>=math.ceil(.7*sum(panel[g][0]==sign for g in ids)),"sign coverage")
        scores=[statistics.mean(panel[g][0]*(math.log2(b[g]+1)-math.log2(a[g]+1)) for g in mapped) for a,b in pairs]
        result[name]={"mapped":len(mapped),"frozen":len(ids),"n_pairs":len(scores),"positive":sum(x>0 for x in scores),"negative":sum(x<0 for x in scores),"ties":sum(x==0 for x in scores),"median_shift":statistics.median(scores),"two_sided_sign_p":exact_sign(scores),"median_interval":median_interval(scores)}
    for k,p in holm({k:v["two_sided_sign_p"] for k,v in result.items()}).items():result[k]["holm_two_panel_p"]=p
    return result

def run(panel_path,pair_path,root,inventory_path,run_id):
    panel=read_panel(panel_path);pairs=fixed_pairs(pair_path)
    inventory=json.loads(inventory_path.read_text())
    require(inventory["scope"]=="full" and inventory["expected"]==inventory["technical_complete"]==inventory["cohort_denominator"]==16,"full cohort only")
    expected={x[k] for x in pairs for k in ("ETH","DEX")}
    require(len(inventory["results"])==16 and {x["run"] for x in inventory["results"]}==expected,"inventory identity")
    require(all(x["status"]=="TECHNICAL_COMPLETE_QC_REVIEW_REQUIRED" for x in inventory["results"]),"incomplete technical inventory")
    require({p.name for p in root.glob("ocular-seven-day-gene-*")}=={"ocular-seven-day-gene-"+r for r in expected},"unexpected/missing run directory")
    for row in inventory["results"]:
        require(row["gene_table_sha256"]==sha(root/("ocular-seven-day-gene-"+row["run"])/"gene_tpm.tsv"),"inventory hash mismatch")
    data=[(load_sample(root,x["ETH"]),load_sample(root,x["DEX"])) for x in pairs]
    rank=evaluate(panel,data)
    for v in rank.values():
        v["median_interval"]=median_interval(v.pop("pair_scores"))
    return {"status":"FULL_COHORT_SENSITIVITY_NOT_INDEPENDENT","scientific_validation":False,"source_run":run_id,"panel_sha256":sha(panel_path),"pair_blob":PAIR_BLOB,"inventory_sha256":sha(inventory_path),"n_pairs":8,"rank_sensitivity":rank,"log2_TPM_plus1_sensitivity":log_tpm_sensitivity(panel,data),"note":"Same historical cohort; neither scale is original DEE2 log2 CPM. Separate two-panel Holm families are descriptive sensitivity analyses, not additional confirmatory tests. No responder interaction or clinical inference."}

def selftest(panel_path,pair_path):
    panel=read_panel(panel_path); fixed_pairs(pair_path)
    a={g:1. if s==1 else 3. for g,(s,e) in panel.items()}
    b={g:3. if s==1 else 1. for g,(s,e) in panel.items()}
    for pairs,expected in [([(a,b)]*8,.015625), ([(a,b)]*7+[(b,a)],.140625), ([(a,a)]*8,1.)]:
        result=log_tpm_sensitivity(panel,pairs)
        require(all(abs(v["holm_two_panel_p"]-expected)<1e-12 for v in result.values()),"Holm parity")
    require(median_interval(list(range(8)))=={"lower":0,"upper":7,"coverage":.9921875},"exact interval parity")
    require(median_interval([1]*5)["lower"] is None,"small n interval")
    try:log_tpm_sensitivity(panel,[(a,{**b,"unexpected":0})])
    except ValueError:pass
    else:raise AssertionError("mismatched universe accepted")
    print(json.dumps({"status":"SYNTHETIC_PASS","scientific_validation":False,"checks":"pair identity, fixed panel, 8/8, 7/8, ties, exact median interval, mismatched universe"}))

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--panel",type=Path,required=True);p.add_argument("--pairs",type=Path,required=True)
    p.add_argument("--selftest",action="store_true");p.add_argument("--root",type=Path)
    p.add_argument("--inventory",type=Path);p.add_argument("--source-run");p.add_argument("--out",type=Path)
    a=p.parse_args()
    if a.selftest:selftest(a.panel,a.pairs)
    else:
        require(all((a.root,a.inventory,a.source_run,a.out)),"full inputs required")
        result=run(a.panel,a.pairs,a.root,a.inventory,a.source_run)
        a.out.write_text(json.dumps(result,indent=2)+"\n");print(json.dumps(result,indent=2))
