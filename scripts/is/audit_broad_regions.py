#!/usr/bin/env python3
"""QC of provisional broad IS regions. Never calls them LD-independent."""
import csv, json, argparse
from pathlib import Path
from collections import defaultdict
def read_tsv(path):
    with path.open() as f: return list(csv.DictReader(f,delimiter="\t"))
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input",type=Path,default=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1"))
    p.add_argument("--out",type=Path)
    a=p.parse_args(); out=a.out or a.input
    regions=read_tsv(a.input/"BROAD_IS_PROVISIONAL_REGIONS.tsv")
    qc=read_tsv(a.input/"DATASET_QC.tsv")
    out.mkdir(parents=True,exist_ok=True)
    seen=set(); errs=[]; by=defaultdict(list)
    for r in regions:
        rid=r["region_id"]
        if rid in seen: errs.append("DUPLICATE_REGION_ID:"+rid)
        seen.add(rid)
        try:
            start,end=int(r["start"]),int(r["end"])
            if start>end or start<1: errs.append("INVALID_INTERVAL:"+rid)
            if int(r["n_p_lt_1e6"])<int(r["n_p_lt_5e8"]): errs.append("GWS_GT_SUGGESTIVE:"+rid)
            if r["region_class"]=="GWS" and int(r["n_p_lt_5e8"])==0: errs.append("GWS_WITH_NO_GWS:"+rid)
        except (ValueError,KeyError): errs.append("BAD_ROW:"+rid);continue
        by[(r["dataset"],r["phenotype"],r["chr"])].append((start,end,rid))
    for key,ivs in by.items():
        ivs.sort()
        for left,right in zip(ivs,ivs[1:]):
            if left[1]>=right[0]: errs.append("OVERLAPPING_REGIONS:"+left[2]+":"+right[2])
    report={"status":"PASS" if not errs else "FAIL","n_sources":len(qc),
      "n_regions":len(regions),"n_gws_regions":sum(r["region_class"]=="GWS" for r in regions),
      "n_suggestive_only":sum(r["region_class"]=="SUGGESTIVE" for r in regions),
      "n_regions_by_dataset":dict(sorted(((q["source"].split("/")[-1],sum(r["dataset"] in q["source"] or q["source"].split("/")[-1].startswith(r["dataset"]) for r in regions)) for q in qc))),
      "issues":errs,"interpretation":"Distance clusters only; not LD-independent or validated causal genes"}
    (out/"BROAD_REGION_QC.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    if errs: raise SystemExit(1)
if __name__=="__main__":main()
