#!/usr/bin/env python3
"""Region-specific EAS LD clumping of significant BBJ/GIGASTROKE study signals.

This is a 1000G-EAS marker correlation screen, NOT true fine mapping and
NOT independent cross-cohort replication (GIGASTROKE studies overlap).
"""
import argparse,csv,gzip,json,subprocess
from collections import defaultdict
from pathlib import Path
from audit_gws_lead_eaf import source_prefix,NEW,OLD
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2")
def read(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def prepare(root):
    expected=read(root/"IS_GWS_STUDY_INPUT_READINESS.tsv")
    if len(expected)!=42:raise ValueError("Exactly 42 source/group pairs needed")
    candidates=defaultdict(dict)
    gws=set()
    with gzip.open(root/"IS_GWS_REFERENCE_MATCHED_VARIANTS.tsv.gz","rt") as f:
        for r in csv.DictReader(f,delimiter="\t"):
            if r["qc_status"] not in ("MATCH_ALT_EFFECT","MATCH_REF_EFFECT"):continue
            try:p=float(r["p"])
            except (ValueError,TypeError):continue
            if p<0 or p>1e-4:continue
            key=(r["dataset"],r["group_id"])
            if p<=5e-8:gws.add(key)
            if p==0:p=1e-300
            record=candidates[key].get(r["variant_id"])
            if record is None or p<record:
                candidates[key][r["variant_id"]]=p
    out=root/"clump_screen_v1";out.mkdir(exist_ok=True,parents=True)
    tasks=[]
    for (ds,g),variant_p in sorted(candidates.items()):
        if (ds,g) not in gws:continue
        tag=ds+"__"+g
        assoc=out/(tag+".assoc")
        with assoc.open("w",newline="") as h:
            w=csv.writer(h,delimiter="\t");w.writerow(["ID","P"])
            for vid,p in sorted(variant_p.items()):
                w.writerow([vid,format(p,".12g")])
        prefix=source_prefix(g,NEW,OLD)
        triplet=all(Path(str(prefix)+e).is_file() for e in (".pgen",".pvar",".psam"))
        tasks.append({"dataset":ds,"group_id":g,
            "significant_study_group":1,"tested_variant_count":len(variant_p),
            "ancestry":"BBJ_JAPANESE" if ds=="BBJ" else "GIGASTROKE_EAS",
            "reference_population":"1000G_EAS_504",
            "reference_prefix":str(prefix),"triplet_complete":int(triplet),
            "assoc_path":str(assoc),"clump_output_prefix":str(out/tag),
            "clump_p1":"5e-8","clump_p2":"1e-4","clump_r2":"0.1","clump_kb":"1000",
            "status":"READY_TO_CLUMP" if triplet else "BLOCKED_NO_REFERENCE"})
    if not tasks:raise RuntimeError("No real significant study/group pairs found")
    path=out/"IS_GWS_CLUMP_EXECUTION_MANIFEST.tsv"
    with path.open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=list(tasks[0]),delimiter="\t")
        w.writeheader();w.writerows(tasks)
    return tasks
def execute(tasks):
    out=[]
    for t in tasks:
        if t["status"]!="READY_TO_CLUMP":raise RuntimeError("Reference missing for "+t["group_id"])
        cmd=["plink2","--pfile",t["reference_prefix"],"--clump",t["assoc_path"],
             "--clump-p1","5e-8","--clump-p2","1e-4","--clump-r2","0.1","--clump-kb","1000",
             "--out",t["clump_output_prefix"]]
        r=subprocess.run(cmd,capture_output=True,text=True,timeout=180)
        if r.returncode!=0:raise RuntimeError("LD clumping failed: "+t["dataset"]+" "+t["group_id"]+" "+r.stdout[-800:])
        clumps=Path(t["clump_output_prefix"]+".clumps")
        with clumps.open() as f:rows=list(csv.DictReader(f,delimiter="\t"))
        out.append({"dataset":t["dataset"],"group_id":t["group_id"],
           "candidate_variants":t["tested_variant_count"],
           "clump_index_count":len(rows),
           "clump_index_ids":";".join(row["ID"] for row in rows),
           "reference":"1000G_EAS_504",
           "scientific_status":"LD_CLUMP_MARKER_SCREEN_NOT_INDEPENDENT_CAUSAL_SIGNALS"})
    target=Path(tasks[0]["assoc_path"]).parent
    with (target/"IS_GWS_LD_CLUMP_RESULTS.tsv").open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=list(out[0]),delimiter="\t")
        w.writeheader();w.writerows(out)
    return out
def main(args):
    tasks=prepare(args.root)
    out={"prepared_study_region_pairs":len(tasks),
         "tested_study_region_pairs":0,
         "status":"PLAN_ONLY",
         "warnings":"Within-study EAS-reference clumps only; no independent locus/causal claims."}
    if args.execute:
        results=execute(tasks)
        out.update({"tested_study_region_pairs":len(results),
           "clump_index_count_sum":sum(x["clump_index_count"] for x in results),
           "pairs_with_multiple_index_clumps":sum(x["clump_index_count"]>1 for x in results),
           "status":"REFERENCE_LD_CLUMP_SCREEN_COMPLETE"})
    folder=args.root/"clump_screen_v1"
    (folder/"IS_GWS_LD_CLUMP_SUMMARY.json").write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))
    return out
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=ROOT)
    p.add_argument("--execute",action="store_true")
    main(p.parse_args())
