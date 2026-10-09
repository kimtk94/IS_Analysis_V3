#!/usr/bin/env python3
"""EUR AIS genome-wide interval discovery without confounding with EAS reference.

Input: independently MD5-verified source, normalized allele-pair GWAS.
GWS p<=5e-8, suggestive p<=1e-6, MAF>=.01.
Distance groupings are PROVISIONAL and cannot be treated as LD-independent loci.
Retain the original 30 EAS/Japanese components and add EUR candidates separately.
"""
import argparse,csv,gzip,json,math
from pathlib import Path
from collections import defaultdict,Counter

BASE=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1")
EUR=Path("/srv/is-analysis/data/is/processed/gigastroke/broad_v1/GCST90104540_AIS_EUR_GRCh37.allele_pairs.tsv.gz")
OUTPUT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v2_ancestry")
CHRS={str(i) for i in range(1,23)}
FIELDS=["group_id","chr","region_start","region_end","window_start","window_end",
    "window_bp","lead_variant","lead_p","phenotypes","has_gws","ancestry",
    "source","provisional_status","linked_eas_components","fine_mapping_status"]
def read(path):
    with path.open() as h:return list(csv.DictReader(h,delimiter="\t"))
def parse(path, distance=1000000):
    sig=defaultdict(list)
    qc=Counter()
    with gzip.open(path,"rt") as f:
        for r in csv.DictReader(f,delimiter="\t"):
            qc["total"]+=1
            if r["chr"] not in CHRS:qc["other_chrom"]+=1;continue
            try:
                pos=int(r["pos"]);p=float(r["p"]);af=float(r["eaf"])
            except (ValueError,TypeError):
                qc["malformed"]+=1;continue
            if not 0<=p<=1 or not 0<af<1 or not math.isfinite(p) or not math.isfinite(af):
                qc["invalid"]+=1;continue
            if min(af,1-af)<.01:
                qc["maf_lt_001"]+=1;continue
            if p>1e-6:continue
            qc["suggestive_or_better"]+=1
            if p<=5e-8:qc["gws"]+=1
            v={"pos":pos,"p":p,"variant_pair_id":r["variant_pair_id"],
               "is_gws":p<=5e-8}
            sig[r["chr"]].append(v)
    groups=[]
    for chrom,rows in sig.items():
        rows.sort(key=lambda x:x["pos"])
        cluster=[];last=None
        for v in rows:
            if cluster and v["pos"]-last>distance:
                groups.append((chrom,cluster));cluster=[]
            cluster.append(v);last=v["pos"]
        if cluster:groups.append((chrom,cluster))
    order=lambda x:(int(x[0]),min(v["pos"] for v in x[1]))
    groups.sort(key=order)
    return groups,dict(qc)
def build(base,eur,out,distance):
    if not eur.is_file():raise FileNotFoundError("MD5-verified EUR normalized file required: "+str(eur))
    original=read(base/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv")
    if len(original)!=30:raise ValueError("Original EAS/Japanese 30-group universe changed")
    original_groups={x["group_id"] for x in original}
    if len(original_groups)!=30:raise ValueError("Duplicate original group ID")
    groups,qc=parse(eur,distance)
    eurs=[]
    for n,(chrom,variants) in enumerate(groups,1):
        lead=min(variants,key=lambda x:(x["p"],x["pos"]))
        lo=min(x["pos"] for x in variants);hi=max(x["pos"] for x in variants)
        links=[r["group_id"] for r in original if r["chr"]==chrom and
               int(r["region_start"])<=hi and int(r["region_end"])>=lo]
        eurs.append({"group_id":f"IS_EUR_AIS_G{n:04d}","chr":chrom,
          "region_start":lo,"region_end":hi,"window_start":max(1,lo-500000),
          "window_end":hi+500000,"window_bp":hi-lo+1000001,
          "lead_variant":lead["variant_pair_id"],"lead_p":lead["p"],
          "phenotypes":"AIS","has_gws":int(any(v["is_gws"] for v in variants)),
          "ancestry":"EUR","source":"GCST90104540",
          "provisional_status":"DISTANCE_CLUSTER_NOT_LD_INDEPENDENT",
          "linked_eas_components":";".join(links),
          "fine_mapping_status":"BLOCKED_PENDING_EUR_LD_AND_REF_ALT_HARMONIZATION"})
    rows=[]
    for r in original:
        item={x:r.get(x,"") for x in FIELDS}
        item.update({"ancestry":"EAS_AND_JAPANESE","source":"BBJ_GIGASTROKE_EAS",
                     "provisional_status":"COORDINATE_MERGE_NOT_LD_INDEPENDENT",
                     "fine_mapping_status":"BLOCKED_PENDING_LD_SIGNAL_QC"})
        rows.append(item)
    rows+=eurs
    out.mkdir(parents=True,exist_ok=True)
    with (out/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,delimiter="\t");w.writeheader();w.writerows(rows)
    with (out/"IS_EUR_AIS_PROVISIONAL_REGIONS.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,delimiter="\t");w.writeheader();w.writerows(eurs)
    report={"status":"COMPLETE_EUR_POSITIONAL_DISCOVERY",
        "input_variant_qc":qc,"original_eas_japanese_components":len(original),
        "eur_distance_clusters":len(eurs),
        "eur_gws_clusters":sum(x["has_gws"] for x in eurs),
        "eur_suggestive_only_clusters":sum(not x["has_gws"] for x in eurs),
        "eur_overlapping_original_components":sum(bool(x["linked_eas_components"]) for x in eurs),
        "combined_rows":len(rows),"distance_threshold_bp":distance,
        "warning":"EUR versus EAS sample overlap possible; intervals are not independent loci; no credible sets"}
    (out/"IS_ANCESTRY_EXPANSION_DISCOVERY_SUMMARY.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    return report
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--base",type=Path,default=BASE)
    p.add_argument("--eur",type=Path,default=EUR)
    p.add_argument("--out",type=Path,default=OUTPUT)
    p.add_argument("--distance-bp",type=int,default=1000000)
    a=p.parse_args()
    build(a.base,a.eur,a.out,a.distance_bp)
