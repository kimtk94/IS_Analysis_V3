#!/usr/bin/env python3
"""Quantify effect of expanding narrow 1KG-EAS LD panels, preserving old QC."""
import argparse,csv,json
from pathlib import Path
from collections import defaultdict
BASE=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1")
def tsv(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def compare(base):
    old=tsv(base/"IS_GWS_GWAS_REFERENCE_OVERLAP.tsv")
    v2=base/"expanded_reference_v2"
    new=tsv(v2/"IS_GWS_GWAS_REFERENCE_OVERLAP.tsv")
    keys=lambda rows:{(r["dataset"],r["group_id"]):r for r in rows}
    a=keys(old);b=keys(new)
    if len(a)!=42 or set(a)!=set(b):
        raise ValueError("Study/group dimensions were altered")
    expanded={"IS_XDATA_G0015","IS_XDATA_G0022"}
    out=[]
    for key in sorted(a):
        prev=a[key];curr=b[key]
        nold=int(prev["matched_oriented_nonpal"])
        nnew=int(curr["matched_oriented_nonpal"])
        if curr["group_id"] in expanded and nnew<nold:
            raise ValueError(f"Expanded genotype reference lost matched SNPs at {key}")
        if curr["group_id"] not in expanded and nnew!=nold:
            raise ValueError(f"Unexpected change in unchanged reference at {key}")
        if int(prev["window_variants"])!=int(curr["window_variants"]):
            raise ValueError(f"GWAS source window count changed at {key}")
        out.append({"dataset":key[0],"group_id":key[1],
            "source_window_variants":int(prev["window_variants"]),
            "old_nonpal_oriented":nold,"new_nonpal_oriented":nnew,
            "newly_oriented":nnew-nold,
            "old_position_missing":int(prev["no_reference_position"]),
            "new_position_missing":int(curr["no_reference_position"]),
            "position_missing_reduced":int(prev["no_reference_position"])-int(curr["no_reference_position"]),
            "expanded_panel":int(curr["group_id"] in expanded)})
    path=v2/"IS_EXPANDED_LD_BEFORE_AFTER.tsv"
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter="\t")
        w.writeheader();w.writerows(out)
    bygroup=defaultdict(lambda:{"old_nonpal":0,"new_nonpal":0,"old_missing":0,"new_missing":0})
    for r in out:
        q=bygroup[r["group_id"]]
        for k,target in [("old_nonpal_oriented","old_nonpal"),("new_nonpal_oriented","new_nonpal"),
             ("old_position_missing","old_missing"),("new_position_missing","new_missing")]:
            q[target]+=r[k]
    report={"n_study_region_pairs":len(out),
       "expanded_groups":sorted(expanded),
       "by_group":dict(sorted(bygroup.items())),
       "newly_oriented_total":sum(x["newly_oriented"] for x in out),
       "newly_oriented_expanded":sum(x["newly_oriented"] for x in out if x["expanded_panel"]),
       "old_oriented_total":sum(x["old_nonpal_oriented"] for x in out),
       "new_oriented_total":sum(x["new_nonpal_oriented"] for x in out),
       "old_missing_total":sum(x["old_position_missing"] for x in out),
       "new_missing_total":sum(x["new_position_missing"] for x in out),
       "interpretation":"More position/allele matches only, not new GWAS association or credible set."}
    (v2/"IS_EXPANDED_LD_BEFORE_AFTER_SUMMARY.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    return report
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--base",type=Path,default=BASE)
    compare(p.parse_args().base)
