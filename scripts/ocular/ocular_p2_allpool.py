#!/usr/bin/env python3
"""P2: inspect all 18 GSE313736 cell barcodes/demux links, aggregate-only output."""
import argparse,csv,gzip,json,tempfile
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone
from ocular_p1_integrity import BASE,download,bar,lines,demux_barcodes

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--output",required=True);args=ap.parse_args()
    a={"utc":datetime.now(timezone.utc).isoformat(),"dataset":"GSE313736",
       "scope":"All nine donor pools, 18 public barcode/demux datasets; matrices C only verified by P1",
       "libraries":[],"pool_pair_summaries":[],"status":"HOLD"}
    with tempfile.TemporaryDirectory(prefix="ocular-p2-") as t:
        tmp=Path(t)
        meta=tmp/"cell_metadata.tsv.gz"
        meta_url=BASE+"series/GSE313nnn/GSE313736/suppl/GSE313736_HTMC_scRNA_cell_metadata.tsv.gz"
        a["metadata_transfer"]=download(meta_url,meta,12*1024*1024)
        m={}; droplet=Counter(); all_donors=set(); cells=0
        with gzip.open(meta,"rt",encoding="utf8",errors="replace") as f:
            rd=csv.DictReader(f,delimiter="\t")
            for r in rd:
                cells+=1
                group=(r["library_id"].strip(),r["DEX Treatment"].strip().lower())
                m.setdefault(group,{})[bar(r["barcode"])]=r["donor_id_demuxlet"].strip()
                if "droplet_type" in r:droplet[r["droplet_type"].strip()]+=1
                if r["donor_id_demuxlet"].strip():all_donors.add(r["donor_id_demuxlet"].strip())
        a["metadata_cells"]=cells;a["metadata_distinct_donors"]=len(all_donors)
        a["metadata_droplet_type_counts"]=dict(droplet)
        paired_all=set(); covered=set(); fail=0
        for i in range(9):
            pool=chr(ord("A")+i);ds={}
            for condition,gsm_id,prefix in (
                ("vehicle",9373276+i,"pool_"+pool),
                ("dex",9373285+i,"pool_dex"+pool)):
                gsm="GSM"+str(gsm_id)
                obj={"pool":pool,"condition":condition,"gsm":gsm,"status":"HOLD"}
                try:
                    barcodes=tmp/(gsm+"_barcodes.tsv.gz");dmx=tmp/(gsm+"_dmx.best.gz")
                    obj["barcodes_transfer"]=download(BASE+"samples/"+gsm[:-3]+"nnn/"+gsm+"/suppl/"+gsm+"_"+prefix+"_barcodes.tsv.gz",barcodes,10*1024*1024)
                    obj["dmx_transfer"]=download(BASE+"samples/"+gsm[:-3]+"nnn/"+gsm+"/suppl/"+gsm+"_"+prefix+"_dmx.best.gz",dmx,10*1024*1024)
                    bs={bar(v) for v in lines(barcodes)}
                    demux=demux_barcodes(dmx)
                    obj["barcode_count"]=len(bs);obj["demux_coverage"]=round(len(bs&demux)/max(1,len(bs)),4)
                    candidates=[]
                    for (lib,treated),barcode_to_donor in m.items():
                        if (condition=="dex") != (treated=="dex"):continue
                        n=len(bs&set(barcode_to_donor))
                        candidates.append((n,lib,treated,len(barcode_to_donor)))
                    if not candidates:raise ValueError("NO_META_TREATMENT")
                    candidates.sort(reverse=True)
                    n,lib,treated,nmeta=candidates[0]
                    obj["matched_library"]=lib;obj["metadata_barcodes"]=nmeta
                    obj["barcode_overlap"]=n;obj["metadata_join_rate"]=round(n/max(1,nmeta),4)
                    donor={d for bc,d in m[(lib,treated)].items()
                           if bc in bs and d and d.lower() not in ("na","doublet","unassigned")}
                    obj["matched_donors"]=len(donor)
                    if nmeta<1 or obj["metadata_join_rate"]<.90:
                        raise ValueError("META_JOIN_BELOW_90_PERCENT")
                    if obj["demux_coverage"]<.50:
                        obj["demux_warning"]="DEMUX_COVERAGE_LT_50_PERCENT; not used as primary donor assignment"
                    elif obj["demux_coverage"]<.90:
                        obj["demux_warning"]="PARTIAL_DEMUX_COVERAGE_LT_90_PERCENT"
                    obj["status"]="BARCODE_METADATA_PASS"
                    ds[condition]=donor
                    covered|=donor
                except Exception as e:
                    fail+=1;ds[condition]=set()
                    obj["reason"]=type(e).__name__+": "+str(e)[:180]
                a["libraries"].append(obj)
            both=ds.get("vehicle",set())&ds.get("dex",set())
            paired_all|=both
            a["pool_pair_summaries"].append({
                "pool":pool,"vehicle_donors":len(ds.get("vehicle",set())),
                "dex_donors":len(ds.get("dex",set())),
                "paired_donors":len(both),
                "both_conditions_valid":all(any(x["pool"]==pool and x["condition"]==c and x["status"]=="BARCODE_METADATA_PASS" for x in a["libraries"]) for c in ("vehicle","dex"))})
        a["aggregates"]={"libraries_checked":len(a["libraries"]),
          "libraries_pass":sum(x["status"]=="BARCODE_METADATA_PASS" for x in a["libraries"]),
          "donors_with_any_barcode_mapping":len(covered),
          "donors_with_paired_conditions":len(paired_all),
          "metadata_total_distinct_donors":len(all_donors),
          "unsuccessful_libraries":fail}
        a["status"]="GO_ALLPOOL_METADATA_ONLY" if fail==0 and len(paired_all)>=50 else "HOLD_ALLPOOL_METADATA"
    Path(args.output).write_text(json.dumps(a,indent=2,ensure_ascii=False)+"\n")
    print("FINAL",a["status"],a["aggregates"])
    for x in a["pool_pair_summaries"]:print("POOL",x)
    for x in a["libraries"]:
        if x.get("demux_warning") or x.get("reason"):
            print("QC_WARNING",x["pool"],x["condition"],x.get("demux_warning",x.get("reason")))
    if a["status"]!="GO_ALLPOOL_METADATA_ONLY":raise SystemExit(2)
if __name__=="__main__":main()
