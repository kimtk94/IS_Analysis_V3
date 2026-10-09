#!/usr/bin/env python3
"""Map all GENCODE v19 genes overlapping expanded GRCh37 IS locus windows.

Inclusive overlap; not causal assignment. Preserve noncoding genes.
"""
import argparse,csv,gzip,json,re
from collections import defaultdict
from pathlib import Path
ATTR=re.compile(r'([A-Za-z_][A-Za-z0-9_]*) "([^"]*)"')
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",type=Path,default=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1"))
    ap.add_argument("--gtf",type=Path,default=Path("/srv/is-analysis/data/is/reference/gencode_v19/gencode.v19.annotation.gtf.gz"))
    args=ap.parse_args()
    with (args.root/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv").open() as f:
        regions=list(csv.DictReader(f,delimiter="\t"))
    bychr=defaultdict(list)
    for r in regions:
        bychr[r["chr"]].append(r)
    hits=[]; n_genes=0
    with gzip.open(args.gtf,"rt") as f:
        for line in f:
            if line.startswith("#"):continue
            c=line.rstrip("\n").split("\t")
            if len(c)<9 or c[2]!="gene":continue
            chrom=c[0].removeprefix("chr")
            if chrom not in bychr:continue
            start,end=int(c[3]),int(c[4])
            a=dict(ATTR.findall(c[8]))
            geneid=a.get("gene_id","")
            symbol=a.get("gene_name","")
            biotype=a.get("gene_type",a.get("gene_biotype","UNDEFINED"))
            n_genes+=1
            for r in bychr[chrom]:
                lo,hi=int(r["window_start"]),int(r["window_end"])
                if end<lo or start>hi:continue
                lead=int(r["lead_variant"].split(":")[1])
                dist=max(start-lead,lead-end,0)
                hits.append({"group_id":r["group_id"],"chr":chrom,
                 "gene_start":start,"gene_end":end,"gene_id":geneid,
                 "gene_symbol":symbol,"biotype":biotype,"strand":c[6],
                 "lead_distance_bp":dist,"region_start":r["region_start"],
                 "region_end":r["region_end"],"annotation_version":"GENCODE_v19_GRCh37",
                 "evidence":"GENOMIC_WINDOW_OVERLAP_ONLY"})
    fields=["group_id","chr","gene_start","gene_end","gene_id","gene_symbol","biotype",
            "strand","lead_distance_bp","region_start","region_end",
            "annotation_version","evidence"]
    with (args.root/"IS_ALL_GENE_WINDOW_UNIVERSE.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t");w.writeheader();w.writerows(hits)
    group_gene=defaultdict(set)
    for x in hits:group_gene[x["group_id"]].add(x["gene_id"])
    report={"regions":len(regions),"annotation_gene_records_scanned":n_genes,
     "region_gene_associations":len(hits),"unique_gene_ids":len(set(x["gene_id"] for x in hits)),
     "coding_gene_ids":len(set(x["gene_id"] for x in hits if x["biotype"]=="protein_coding")),
     "noncoding_or_other_gene_ids":len(set(x["gene_id"] for x in hits if x["biotype"]!="protein_coding")),
     "regions_without_any_annotated_gene":sorted(set(r["group_id"] for r in regions)-set(group_gene)),
     "interpretation":"Overlapping genes only. Regulatory distal genes and QTL mapping not yet incorporated."}
    (args.root/"IS_ALL_GENE_UNIVERSE_SUMMARY.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
if __name__=="__main__":main()
