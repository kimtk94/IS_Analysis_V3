#!/usr/bin/env python3
"""Read back the observed GSE313736 metadata schema, returning aggregate design/QC only."""
import argparse,collections,csv,gzip,hashlib,json,re
from pathlib import Path
from ocular_gse313736_design_gate import fetch,BASE,FILES

PIN="41e4a6bf5792fe58e023b92f0129119dcf2b49877f62dd580b9dab5de48f5192"
REQUIRED=("library_id","DEX Treatment","barcode","droplet_type","donor_id_demuxlet","sampleid",
"scpred_vanzyl2022","celltypist_chen2025_coarse","celltypist_chen2025_fine")
def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    path=a.root/FILES[0];source=fetch(BASE+"/suppl/"+FILES[0],path,10000000)
    if source["sha256"]!=PIN or source["bytes"]!=6592040:raise ValueError("Pinned source changed")
    arms=collections.Counter();droplet=collections.Counter();cross=collections.Counter()
    by=collections.defaultdict(collections.Counter);libraries=collections.defaultdict(set)
    labels=collections.defaultdict(collections.Counter);seen=set();n=0;dups=0
    # Exact pseudonymous sample IDs used only transiently. No IDs or individual records emitted.
    with gzip.open(path,"rt",encoding="utf-8-sig",newline="") as f:
        r=csv.DictReader(f,delimiter="\t")
        if not set(REQUIRED)<=set(r.fieldnames or []):raise ValueError("Source schema drift")
        for x in r:
            n+=1
            if None in x or any(x[k] is None for k in REQUIRED):raise ValueError("Malformed row")
            if not re.fullmatch(r"HTMC[-_]?(\d+)",x["sampleid"],re.I):raise ValueError("Unexpected sample token")
            k=(x["library_id"],x["barcode"])
            if k in seen:dups+=1
            seen.add(k)
            # Exact original categorical values for source review.
            arm=x["DEX Treatment"].strip()
            dt=x["droplet_type"].strip()
            if len(arm)>60 or len(dt)>60:raise ValueError("Unexpected category")
            arms[arm]+=1;droplet[dt]+=1;cross[arm,dt]+=1
            by[x["sampleid"]][arm,dt]+=1
            libraries[arm].add(x["library_id"])
            for col in ("scpred_vanzyl2022","celltypist_chen2025_coarse","celltypist_chen2025_fine"):
                labels[col][x[col]]+=1
    if n!=169398 or dups:raise ValueError("Row count/unique library-barcode gate")
    # Report exact two source levels only; do not assign DEX orientation before inspecting labels.
    levels=sorted(arms);paired=[]
    if len(levels)==2:
        for mode in ["all_rows"]+sorted(droplet):
            vals={d:{arm:sum(v for (aa,dd),v in q.items() if aa==arm and (mode=="all_rows" or dd==mode)) for arm in levels} for d,q in by.items()}
            row={"subset":mode,"levels":levels,"donors_any":sum(any(q.values()) for q in vals.values()),
                 "donors_by_level":{arm:sum(q[arm]>0 for q in vals.values()) for arm in levels}}
            for cut in (1,50,100,200):row["paired_min_"+str(cut)]=sum(min(q.values())>=cut for q in vals.values())
            paired.append(row)
    out={"status":"EXACT_METADATA_ROW_AND_PAIR_COVERAGE_NOT_EXPRESSION","source":source,
         "rows":n,"unique_library_barcode_keys":len(seen),"unique_sample_ids":len(by),
         "source_treatment_levels":dict(arms),"source_droplet_levels":dict(droplet),
         "treatment_by_droplet":[{"level":aa,"droplet":dd,"rows":v} for (aa,dd),v in sorted(cross.items())],
         "library_counts_by_level":{k:len(v) for k,v in libraries.items()},"pair_coverage":paired,
         "annotation_category_counts":{k:dict(v) for k,v in labels.items()},
         "limits":["Only aggregate publicly released metadata; no clinical covariates, identifiers, raw genotype or expression retained.",
                   "Source sample labels are not independently verified donor uniqueness; paired counts remain subject to cell QC and barcode/matrix alignment.",
                   "169398 metadata rows, not the publication's full 316799 singlets. No inference that all reported donors or cells retained.",
                   "No expression analysis or P values. Need exact treatment, singlet and cell-type labels before scoring."]}
    path.unlink(missing_ok=True);a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n");print(json.dumps(out,indent=2,allow_nan=False))
if __name__=="__main__":main()
