#!/usr/bin/env python3
"""Read-only inventory of LD reference coverage for cross-study IS regions.

This never infers LD-independent loci from physical distance.
Existing BBJ-L001..L004 regional pfiles are classified as local-only references.
"""
import argparse,csv,json
from collections import defaultdict
from pathlib import Path

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",type=Path,default=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1"))
    ap.add_argument("--reference",type=Path,default=Path("/srv/is-analysis/data/is/ld_reference/1kg_eas_grch37/pgen_gt"))
    args=ap.parse_args()
    with (args.root/"CROSS_DATASET_INTERVAL_GROUPS.tsv").open() as f: groups=list(csv.DictReader(f,delimiter="\t"))
    refs=[]
    for f in sorted(args.reference.glob("*.pvar")):
        prefix=f.with_suffix("")
        if not all(Path(str(prefix)+ext).exists() for ext in (".pgen",".psam")): continue
        chroms=defaultdict(lambda:[None,None,0])
        with f.open() as fh:
            for line in fh:
                if line.startswith("#"): continue
                a=line.rstrip("\n").split("\t")
                if len(a)<2:continue
                try: pos=int(a[1])
                except ValueError:continue
                chr_=a[0].removeprefix("chr")
                cur=chroms[chr_]
                cur[0]=pos if cur[0] is None else min(cur[0],pos)
                cur[1]=pos if cur[1] is None else max(cur[1],pos)
                cur[2]+=1
        for chrom,(lo,hi,n) in chroms.items():
            refs.append({"prefix":str(prefix),"chr":chrom,"start":lo,"end":hi,"n_variants":n})
    result=[]
    for g in groups:
        chrom=g["chr"];lo=int(g["start"]);hi=int(g["end"])
        matches=[r for r in refs if r["chr"]==chrom and r["start"]<=hi and r["end"]>=lo]
        full=[r for r in matches if r["start"]<=lo and r["end"]>=hi]
        status="NO_LOCAL_REFERENCE"
        if full:status="BOUNDARY_COVERED_NEEDS_VARIANT_AND_LD_QC"
        elif matches:status="PARTIAL_INTERVAL_REFERENCE"
        result.append({"group_id":g["group_id"],"chr":chrom,"start":lo,"end":hi,
           "n_reference_overlaps":len(matches),"n_boundary_covered":len(full),
           "reference_prefixes":";".join(r["prefix"] for r in matches),
           "ld_readiness":status,
           "independent_locus_claim":"NOT_EVALUATED"})
    fields=list(result[0]) if result else ["group_id"]
    with (args.root/"CROSS_DATASET_LD_COVERAGE_AUDIT.tsv").open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t");w.writeheader();w.writerows(result)
    summary={"groups":len(result),"references":len(refs),"no_reference":sum(x["ld_readiness"]=="NO_LOCAL_REFERENCE" for x in result),
      "partial":sum(x["ld_readiness"]=="PARTIAL_INTERVAL_REFERENCE" for x in result),
      "boundary_covered":sum(x["ld_readiness"]=="BOUNDARY_COVERED_NEEDS_VARIANT_AND_LD_QC" for x in result),
      "caveat":"Regional pfiles from historic BBJ loci cannot establish LD coverage for unrelated genome-wide regions. Boundary overlap alone is not variant overlap."}
    (args.root/"LD_REFERENCE_COVERAGE_SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
if __name__=="__main__":main()
