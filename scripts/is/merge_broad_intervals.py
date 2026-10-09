#!/usr/bin/env python3
"""Unify overlapping provisional intervals across studies without asserting LD independence."""
import argparse,csv,json
from collections import defaultdict
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument("--root",type=Path,default=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1"))
a=p.parse_args()
with (a.root/"BROAD_IS_PROVISIONAL_REGIONS.tsv").open() as f: rows=list(csv.DictReader(f,delimiter="\t"))
by=defaultdict(list)
for r in rows: by[r["chr"]].append(r)
groups=[]
for chrom,items in by.items():
    items.sort(key=lambda r:(int(r["start"]),int(r["end"])))
    active=[]; last=-1
    for r in items:
        if active and int(r["start"])>last: groups.append((chrom,active));active=[];last=-1
        active.append(r);last=max(last,int(r["end"]))
    if active:groups.append((chrom,active))
def ck(chrom): return (0,int(chrom)) if chrom.isdigit() else (1,chrom)
groups.sort(key=lambda x:(ck(x[0]),min(int(r["start"]) for r in x[1])))
out=[]
for i,(chrom,items) in enumerate(groups,1):
    lead=min(items,key=lambda r:float(r["lead_p"]))
    out.append({"group_id":f"IS_XDATA_G{i:04d}","chr":chrom,
      "start":min(int(r["start"]) for r in items),"end":max(int(r["end"]) for r in items),
      "n_source_regions":len(items),"n_datasets":len(set(r["dataset"] for r in items)),
      "n_phenotypes":len(set(r["phenotype"] for r in items)),
      "datasets":";".join(sorted(set(r["dataset"] for r in items))),
      "phenotypes":";".join(sorted(set(r["phenotype"] for r in items))),
      "lead_p":lead["lead_p"],"lead_variant":lead["lead_variant"],
      "has_gws":int(any(r["region_class"]=="GWS" for r in items)),
      "source_region_ids":";".join(r["region_id"] for r in items),
      "status":"CROSS_DATASET_INTERVAL_COMPONENT_NOT_LD_LOCUS"})
with (a.root/"CROSS_DATASET_INTERVAL_GROUPS.tsv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(out[0]) if out else ["group_id"],delimiter="\t");w.writeheader();w.writerows(out)
summary={"input_source_regions":len(rows),"interval_components":len(out),
 "components_with_gws":sum(x["has_gws"] for x in out),
 "components_multi_dataset":sum(x["n_datasets"]>1 for x in out),
 "status":"PROVISIONAL_INTERVAL_MERGE_ONLY"}
(a.root/"CROSS_DATASET_GROUPS_SUMMARY.json").write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
