#!/usr/bin/env python3
"""Fixed-panel paired percentile score. Proposed validation protocol; no real data run in selftest."""
import argparse
import csv
import hashlib
import json
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path

PIN = "d4b266d833853d66f9990e4b00c1d96c0cd0d9db166202cfdf6a4ffd59379d79"

def require(ok, message):
    if not ok:
        raise ValueError(message)

def read_panel(path):
    require(hashlib.sha256(path.read_bytes()).hexdigest() == PIN, "frozen panel hash mismatch")
    with path.open(newline="", encoding="utf8") as f:
        rd = csv.DictReader(f, delimiter="\t")
        require(rd.fieldnames == ["gene_id", "symbol", "direction", "ecm_member"], "panel schema")
        rows = list(rd)
    panel = {r["gene_id"]: (int(r["direction"]), int(r["ecm_member"])) for r in rows}
    require(len(panel) == len(rows) == 241, "panel cardinality")
    require(all(re.fullmatch(r"ENSG\d+", g) and s in (-1,1) and e in (0,1) for g,(s,e) in panel.items()), "panel identity")
    require(sum(s==1 for s,e in panel.values()) == 117 and sum(e for s,e in panel.values()) == 16 and sum(s==1 and e==1 for s,e in panel.values()) == 7, "panel counts")
    return panel

def rank(v):
    order = sorted(v, key=v.__getitem__)
    out = {}
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and v[order[j]] == v[order[i]]:
            j += 1
        for k in order[i:j]:
            out[k] = (i+j)/(2*len(order))
        i = j
    return out

def exact_sign(values):
    p = sum(v>0 for v in values)
    n = sum(v<0 for v in values)
    return min(1., 2*sum(math.comb(p+n,k) for k in range(min(p,n)+1))/2**(p+n)) if p+n else 1.

def holm(p):
    out = {}
    running = 0.
    for i,(k,v) in enumerate(sorted(p.items(), key=lambda kv:kv[1])):
        running = max(running, min(1., (len(p)-i)*v))
        out[k] = running
    return out

def evaluate(panel, pairs):
    require(bool(pairs), "no pairs")
    first = pairs[0][0]
    universe = set(first)
    require(bool(universe), "empty gene universe")
    for a,b in pairs:
        for sample in (a,b):
            require(set(sample)==universe, "different gene universe; do not fill absent rows")
            require(all(math.isfinite(v) and v>=0 for v in sample.values()), "invalid expression")
    # Identical reference universe for all arms; retain zero-valued genes and average tied ranks.
    ranks = [(rank(a),rank(b)) for a,b in pairs]
    out = {}
    for name, ids in (("all241",list(panel)),("ECM16",[g for g,(s,e) in panel.items() if e])):
        mapped = [g for g in ids if g in universe]
        require(len(mapped)>=math.ceil(.75*len(ids)), name+" overall coverage")
        for sign in (-1,1):
            require(sum(panel[g][0]==sign for g in mapped)>=math.ceil(.7*sum(panel[g][0]==sign for g in ids)), name+" direction coverage")
        scores = [statistics.mean((b[g]-a[g])*panel[g][0] for g in mapped) for a,b in ranks]
        out[name] = {"mapped":len(mapped),"frozen":len(ids),"n_pairs":len(scores),"positive":sum(x>0 for x in scores),"negative":sum(x<0 for x in scores),"ties":sum(x==0 for x in scores),"median_shift":statistics.median(scores),"mean_shift":statistics.mean(scores),"two_sided_sign_p":exact_sign(scores),"pair_scores":scores}
    for name,p in holm({k:v["two_sided_sign_p"] for k,v in out.items()}).items():
        out[name]["holm_two_panel_sensitivity_p"]=p
    return out

def selftest(panel):
    require(rank({"a":4,"b":4,"c":1,"d":9})=={"a":.5,"b":.5,"c":.125,"d":.875},"midranks")
    require(exact_sign([1]*5)==.0625 and exact_sign([0]*5)==1.,"sign test")
    require(holm({"a":.0625,"b":.0625})=={"a":.125,"b":.125},"Holm")
    control={g:1. if s==1 else 3. for g,(s,e) in panel.items()}
    dex={g:3. if s==1 else 1. for g,(s,e) in panel.items()}
    result=evaluate(panel,[(control,dex)]*5)
    require(all(v["positive"]==5 and v["two_sided_sign_p"]==.0625 for v in result.values()),"synthetic positive")
    result=evaluate(panel,[(dex,control)]*5)
    require(all(v["negative"]==5 for v in result.values()),"synthetic reverse")
    result=evaluate(panel,[(control,control)]*5)
    require(all(v["ties"]==5 and v["two_sided_sign_p"]==1 for v in result.values()),"synthetic ties")
    for a,b in [({**control,"extra":1.},dex),({**control,next(iter(control)):float("nan")},dex),({g:v for g,v in control.items() if panel[g][0]==1},{g:v for g,v in dex.items() if panel[g][0]==1})]:
        try:
            evaluate(panel,[(a,b)])
        except ValueError:
            pass
        else:
            raise AssertionError("expected fail closed")
    print(json.dumps({"status":"SYNTHETIC_TEST_PASS","scientific_validation":False,"panel_sha256":PIN,"tests":9,"note":"No expression inputs or biological outcomes; proposed rank universe and two-panel sensitivity require protocol review."}))

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--panel",type=Path,required=True)
    p.add_argument("--selftest",action="store_true",required=True)
    args=p.parse_args()
    selftest(read_panel(args.panel))
