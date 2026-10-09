#!/usr/bin/env python3
"""Contrast original IS anchor genes with unbiased positional candidates.
No selection cutoff, no causal/ranking claim.
"""
import argparse,csv,json
from collections import defaultdict
from pathlib import Path
def read(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",type=Path,default=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1"))
    a=ap.parse_args()
    genes=read(a.root/"IS_ALL_GENE_WINDOW_UNIVERSE.tsv")
    anchors=read(a.root/"LEGACY_ANCHOR_CANDIDATES.tsv")
    byname=defaultdict(list)
    for r in genes:byname[r["gene_symbol"].upper()].append(r)
    result=[]
    for arow in anchors:
        symbol=arow["gene"].upper()
        matches=byname.get(symbol,[])
        result.append({"anchor_gene":arow["gene"],"legacy_locus":arow["locus"],
          "found_in_expanded_windows":int(bool(matches)),
          "group_ids":";".join(sorted(set(x["group_id"] for x in matches))),
          "minimum_lead_distance_bp":min((int(x["lead_distance_bp"]) for x in matches),default=""),
          "mapping_status":"WINDOW_OVERLAP_ONLY" if matches else "NOT_IN_CURRENT_GWAS_WINDOWS_NOT_BIOLOGICAL_NEGATIVE"})
    with (a.root/"IS_LEGACY_ANCHOR_COVERAGE.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,delimiter="\t",fieldnames=list(result[0]));w.writeheader();w.writerows(result)
    counts=defaultdict(int)
    for g in genes:counts[g["group_id"]]+=1
    summary={"legacy_anchor_count":len(result),
        "anchors_overlapping_expanded_windows":sum(x["found_in_expanded_windows"] for x in result),
        "unique_positional_genes":len(set(r["gene_id"] for r in genes)),
        "gene_locus_pairs":len(genes),
        "regions_with_genes":len(counts),
        "gene_count_per_region":dict(sorted(counts.items())),
        "caution":"Positional overlap is hypothesis generation; no causal gene, LD or molecular QTL verdict"}
    (a.root/"IS_CANDIDATE_COVERAGE_SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
if __name__=="__main__":main()
