#!/usr/bin/env python3
"""Evidence coverage matrix for ALL discovered positional genes.

Missing evidence means NOT_TESTED, never negative. Prior deep-screen gene selection
must not be confused with genome-wide molecular screening.
"""
import csv,json,argparse
from pathlib import Path
from collections import defaultdict
def load(path):
    if not path.exists():return []
    with path.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def numeric(value):
    try:return float(value)
    except (TypeError,ValueError):return None
def main():
    a=argparse.ArgumentParser()
    a.add_argument("--root",type=Path,default=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1"))
    a.add_argument("--functional",type=Path,default=Path("/srv/is-analysis/results/is/stage5_functional"))
    o=a.parse_args()
    genes=load(o.root/"IS_ALL_GENE_WINDOW_UNIVERSE.tsv")
    anchors=load(o.root/"LEGACY_ANCHOR_CANDIDATES.tsv")
    coloc=load(o.functional/"phase9c_convergence/COLOC_ABF_MASTER_ANNOTATED_V2.tsv")
    human=load(o.functional/"phase9f_e_human_r3_4_import/HUMAN_TARGET_GENE_SUMMARY.tsv")
    mouse=load(o.functional/"phase9f_c_target_expression/CELLTYPE_TARGET_EVIDENCE_MASTER.tsv")
    anchor_symbols={x["gene"].upper() for x in anchors}
    names={}
    for r in genes:
        key=r["gene_symbol"].upper()
        if not key:continue
        names.setdefault(key,{"gene_symbol":r["gene_symbol"],"gene_ids":set(),"group_ids":set(),"biotypes":set(),"minimum_distance":None})
        x=names[key];x["gene_ids"].add(r["gene_id"]);x["group_ids"].add(r["group_id"]);x["biotypes"].add(r["biotype"])
        n=int(r["lead_distance_bp"]);x["minimum_distance"]=n if x["minimum_distance"] is None else min(n,x["minimum_distance"])
    for s in anchor_symbols:
        names.setdefault(s,{"gene_symbol":s,"gene_ids":set(),"group_ids":set(),"biotypes":set(),"minimum_distance":None})
    col=defaultdict(list)
    for r in coloc:
        key=(r.get("gene_symbol") or "").upper()
        if key and r.get("status")=="PASS":col[key].append(r)
    hum={r["gene"].upper():r for r in human if r.get("gene")}
    mus=defaultdict(list)
    for r in mouse:
        key=(r.get("human_gene") or "").upper()
        if key:mus[key].append(r)
    output=[]
    for key,item in sorted(names.items()):
        valid=[(numeric(r.get("PP.H4")),r) for r in col.get(key,[])]
        valid=[x for x in valid if x[0] is not None]
        best=max(valid,key=lambda x:x[0]) if valid else None
        h=hum.get(key);m=mus.get(key,[])
        output.append({"gene_symbol":item["gene_symbol"],"gene_ids":";".join(sorted(item["gene_ids"])),
          "group_ids":";".join(sorted(item["group_ids"])),"biotypes":";".join(sorted(item["biotypes"])),
          "min_lead_distance_bp":item["minimum_distance"] if item["minimum_distance"] is not None else "",
          "legacy_anchor":int(key in anchor_symbols),
          "positional_status":"WINDOW_MAPPED" if item["group_ids"] else "LEGACY_ANCHOR_OUTSIDE_CURRENT_WINDOWS",
          "coloc_status":"TESTED" if valid else "NOT_TESTED",
          "coloc_best_pp_h4":best[0] if best else "",
          "coloc_best_dataset":best[1].get("dataset_key","") if best else "",
          "coloc_best_pp_h3":best[1].get("PP.H3","") if best else "",
          "human_atlas_status":"PRESENT" if h and h.get("present","").upper()=="TRUE" else ("NOT_DETECTED_IN_SELECTED_ATLAS" if h else "NOT_TESTED"),
          "human_top_celltype":h.get("top_celltype","") if h else "",
          "human_detection_fraction":h.get("top_detection_fraction","") if h else "",
          "mouse_status":"DESCRIPTIVE_DATA_AVAILABLE" if m else "NOT_TESTED",
          "mouse_top_celltypes":";".join(sorted(set(x.get("top_celltype","") for x in m if x.get("top_celltype")))),
          "causal_status":"NOT_ESTABLISHED"})
    fields=list(output[0])
    with (o.root/"IS_GENE_EVIDENCE_COVERAGE_MATRIX.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,delimiter="\t",fieldnames=fields);w.writeheader();w.writerows(output)
    summary={"candidate_genes":len(output),"legacy_anchor_rows":sum(int(x["legacy_anchor"]) for x in output),
      "window_mapped":sum(x["positional_status"]=="WINDOW_MAPPED" for x in output),
      "coloc_tested":sum(x["coloc_status"]=="TESTED" for x in output),
      "human_atlas_assessed":sum(x["human_atlas_status"]!="NOT_TESTED" for x in output),
      "mouse_assessed":sum(x["mouse_status"]!="NOT_TESTED" for x in output),
      "warning":"Current coloc/cell-type sources are targeted legacy screens; coverage flags do not constitute evidence of no effect"}
    (o.root/"IS_GENE_EVIDENCE_COVERAGE_SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
if __name__=="__main__":main()
