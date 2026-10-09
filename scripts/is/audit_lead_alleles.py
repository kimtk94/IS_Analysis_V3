#!/usr/bin/env python3
"""Allele-aware lead validation against local regional EAS PVARs; no LD causal claims."""
import argparse,csv,json
from collections import defaultdict
from pathlib import Path

def run(root, reference):
    with (root/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv").open() as h:
        entries=list(csv.DictReader(h,delimiter="\t"))
    index=defaultdict(list)
    for p in sorted(reference.glob("*.pvar")):
        stem=p.with_suffix("")
        if not all(Path(str(stem)+ext).exists() for ext in (".psam",".pgen")):continue
        with p.open() as h:
            for line in h:
                if line.startswith("#"):continue
                t=line.rstrip("\n").split("\t")
                if len(t)<5:continue
                try:position=int(t[1])
                except ValueError:continue
                index[(t[0].removeprefix("chr"),position)].append((t[3].upper(),t[4].upper().split(","),p.name))
    result=[]
    for row in entries:
        chrom,pos,ref,alt=row["lead_variant"].split(":")[:4]
        match=index.get((chrom,int(pos)),[])
        same=[name for r,alts,name in match if r==ref.upper() and alt.upper() in alts]
        reverse=[name for r,alts,name in match if r==alt.upper() and ref.upper() in alts]
        if same:status="EXACT_REF_ALT_PRESENT"
        elif reverse:status="SWAPPED_REF_ALT_PRESENT_REQUIRES_ORIENTATION"
        elif match:status="POSITION_ONLY_ALLELES_DIFFER"
        else:status="NO_POSITION_MATCH"
        result.append({"group_id":row["group_id"],"lead_variant":row["lead_variant"],
                       "status":status,"exact_references":";".join(same),
                       "swapped_references":";".join(reverse),
                       "position_references":";".join(sorted(set(x[2] for x in match))),
                       "fine_mapping_ready":"NO",
                       "reason":"LD sample, variant overlap, LD-matrix and effect-allele orientation remain unverified"})
    with (root/"LEAD_ALLELE_REFERENCE_AUDIT.tsv").open("w",newline="") as h:
        w=csv.DictWriter(h,delimiter="\t",fieldnames=list(result[0]));w.writeheader();w.writerows(result)
    counts={x:sum(r["status"]==x for r in result) for x in sorted(set(r["status"] for r in result))}
    summary={"total":len(result),"status_counts":counts,
             "interpretation":"Allele match is not sufficient for fine-mapping readiness"}
    (root/"LEAD_ALLELE_REFERENCE_SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1"))
    p.add_argument("--reference",type=Path,default=Path("/srv/is-analysis/data/is/ld_reference/1kg_eas_grch37/pgen_gt"))
    a=p.parse_args()
    run(a.root,a.reference)
