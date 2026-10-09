#!/usr/bin/env python3
"""Non-destructive broad IS discovery: BBJ and available GIGASTROKE EAS.
Outputs are exploratory, not LD-independent loci or causal genes."""
import csv, gzip, json, math
from pathlib import Path
from collections import defaultdict

BASE=Path("/srv/is-analysis")
OUT=BASE/"results/is/stage5_functional/broad_discovery_v1"
OUT.mkdir(parents=True, exist_ok=True)
inputs=[BASE/"data/is/processed/japan/bbj/BBJ_IS_GRCh37.canonical.tsv.gz"]
inputs+=sorted((BASE/"data/is/processed/gigastroke/eas").glob("*.canonical.tsv.gz"))
thresholds=[("gws",5e-8),("suggestive",1e-6)]
hits=defaultdict(list)
qc=[]
for path in inputs:
    counts={"total":0,"invalid":0,"gws":0,"suggestive":0}
    with gzip.open(path,"rt") as handle:
        reader=csv.DictReader(handle,delimiter="\t")
        required={"dataset","phenotype","chr","pos","p","variant_id","effect_allele","other_allele"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Missing required columns in {path}")
        for r in reader:
            counts["total"]+=1
            try:
                p=float(r["p"]); pos=int(r["pos"])
                af=float(r["eaf"]) if r.get("eaf") not in ("",None,"NA") else None
                info=float(r["info"]) if r.get("info") not in ("",None,"NA") else None
                if not (0<p<=1 and pos>0 and (af is None or 0<af<1)):
                    raise ValueError("invalid p/position/frequency")
            except (ValueError,TypeError):
                counts["invalid"]+=1; continue
            if info is not None and info<0.7: continue
            if af is not None and min(af,1-af)<0.01: continue
            if p>1e-6: continue
            tier="gws" if p<=5e-8 else "suggestive"
            counts[tier]+=1
            chrom=r["chr"].removeprefix("chr")
            if chrom not in [str(i) for i in range(1,23)]: continue
            hits[(r["dataset"],r["phenotype"],chrom)].append({
                "dataset":r["dataset"],"phenotype":r["phenotype"],
                "ancestry":r.get("ancestry",""),"chr":chrom,"pos":pos,"p":p,
                "variant_id":r["variant_id"],"beta":r.get("beta",""),
                "effect_allele":r["effect_allele"],"other_allele":r["other_allele"],
                "tier":tier
            })
    qc.append({"source":str(path),**counts})
# Distance-only regional clusters are clearly labeled provisional, not LD-independent.
regions=[]
for (ds,phen,chrom), arr in sorted(hits.items()):
    arr.sort(key=lambda r:r["pos"])
    blocks=[]; block=[]
    for r in arr:
        if block and r["pos"]-block[-1]["pos"]>1_000_000:
            blocks.append(block);block=[]
        block.append(r)
    if block:blocks.append(block)
    for ix,block in enumerate(blocks,1):
        lead=min(block,key=lambda r:r["p"])
        regions.append({"region_id":f"{ds}_{phen}_CHR{chrom}_{ix:04d}",
           "dataset":ds,"phenotype":phen,"ancestry":lead["ancestry"],"chr":chrom,
           "start":min(r["pos"] for r in block),"end":max(r["pos"] for r in block),
           "lead_variant":lead["variant_id"],"lead_p":lead["p"],
           "n_p_lt_1e6":len(block),"n_p_lt_5e8":sum(x["p"]<=5e-8 for x in block),
           "region_class":"GWS" if lead["p"]<=5e-8 else "SUGGESTIVE",
           "status":"PROVISIONAL_DISTANCE_CLUSTER_REQUIRES_LD"})
def save(filename, rows, fields):
    with (OUT/filename).open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",extrasaction="ignore")
        w.writeheader();w.writerows(rows)
fields=list(regions[0]) if regions else ["region_id","dataset","phenotype","chr","start","end"]
save("BROAD_IS_PROVISIONAL_REGIONS.tsv",regions,fields)
save("DATASET_QC.tsv",qc,["source","total","invalid","gws","suggestive"])
# Preserve all genes across original 4 BBJ windows, including noncoding.
universe=BASE/"results/is/stage5_functional/phase9/GRCH37_LOCUS_GENE_UNIVERSE.tsv"
genes=[]
if universe.exists():
    with universe.open() as fh:
        genes=list(csv.DictReader(fh,delimiter="\t"))
if genes:
    save("LEGACY_BBJ_FULL_GENE_UNIVERSE.tsv",genes,list(genes[0]))
legacy=BASE/"results/is/stage5_functional/phase9d_literature_benchmark/GENE_EVIDENCE_LITERATURE_MASTER.tsv"
anchor=[]
if legacy.exists():
    with legacy.open() as fh: anchor=list(csv.DictReader(fh,delimiter="\t"))
if anchor:
    save("LEGACY_ANCHOR_CANDIDATES.tsv",anchor,list(anchor[0]))
summary={"sources":len(inputs),"qc":qc,"provisional_regions":len(regions),
 "gws_regions":sum(r["region_class"]=="GWS" for r in regions),
 "suggestive_only_regions":sum(r["region_class"]=="SUGGESTIVE" for r in regions),
 "legacy_gene_rows":len(genes),"legacy_anchor_genes":len(anchor),
 "warnings":["Not LD-clumped, not independent loci","Not genome-wide discovery until GWAS sources are complete",
 "No nearest-five-gene filtering; gene mapping for expanded regions remains pending",
 "EAS GIGASTROKE only; source coverage is not pan-ancestry"]}
(OUT/"BROAD_DISCOVERY_STATUS.json").write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
