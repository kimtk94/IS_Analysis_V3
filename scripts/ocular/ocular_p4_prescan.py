#!/usr/bin/env python3
"""P4 inexpensive preflight: BeamA score>=.75 donor pairing vs 18 10x barcodes.

No matrices or individual donor identities persist. Gate before ~6.5 GB reload.
"""
import argparse,csv,gzip,json,tempfile
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone
from ocular_p1_integrity import BASE,bar,download,file_url,lines

def check(cutoff=.75,minimum=30):
    result={"dataset":"GSE313736","utc":datetime.now(timezone.utc).isoformat(),
            "status":"HOLD","cell_type":"BeamA","confidence":cutoff,
            "min_cells_per_arm":minimum,"libraries":[]}
    with tempfile.TemporaryDirectory() as td:
        root=Path(td);meta=root/"cell_meta.gz"
        download(BASE+"series/GSE313nnn/GSE313736/suppl/GSE313736_HTMC_scRNA_cell_metadata.tsv.gz",
                 meta,12*1024*1024)
        group={};before={};eligible=Counter();all_metadata=Counter()
        with gzip.open(meta,"rt",encoding="utf8",errors="replace") as f:
            rd=csv.DictReader(f,delimiter="\t")
            for r in rd:
                lib=r["library_id"].strip()
                arm=r["DEX Treatment"].strip().lower()
                if arm not in ("dex","no-dex"):continue
                all_metadata[(lib,arm)]+=1
                if r["droplet_type"]!="SNG" or r["celltypist_chen2025_fine"]!="BeamA":continue
                try:confidence=float(r["celltypist_chen2025_fine_score"])
                except (ValueError,TypeError):continue
                if confidence<cutoff:continue
                donor=r["donor_id_demuxlet"].strip()
                if not donor or donor.lower() in ("na","doublet","unassigned"):continue
                group.setdefault((lib,arm),{})[bar(r["barcode"])]=donor
                eligible[(donor,arm)]+=1
        donor_ids=sorted(set(d for d,arm in eligible if arm=="dex")&
                         set(d for d,arm in eligible if arm=="no-dex"))
        prior=[d for d in donor_ids if eligible[d,"dex"]>=minimum and
               eligible[d,"no-dex"]>=minimum]
        result["prior_eligible_paired_donors"]=len(prior)
        selected=set(prior)
        joined=Counter()
        for pi in range(9):
            p=chr(65+pi)
            for arm,gsm,prefix,lib in [
              ("no-dex",f"GSM{9373276+pi}",f"pool_{p}",f"HTMC_pool-{p}_GRCh38"),
              ("dex",f"GSM{9373285+pi}",f"pool_dex{p}",f"HTMC_pool-dex{p}_GRCh38")]:
                fp=root/(gsm+".gz")
                download(file_url(gsm,prefix,"barcodes.tsv.gz"),fp,12*1024*1024)
                barcodes=set(map(bar,lines(fp)))
                mp=group.get((lib,arm),{})
                candidates={k:v for k,v in mp.items() if v in selected}
                hit=barcodes&set(candidates)
                for bc in hit:joined[(candidates[bc],arm)]+=1
                recall=len(hit)/max(1,len(candidates))
                result["libraries"].append({
                    "gsm":gsm,"pool":p,"arm":arm,
                    "selected_meta":len(candidates),
                    "selected_barcodes_present":len(hit),
                    "recall":round(recall,6),
                    "warning":"RECALL_UNDER_95" if recall<.95 else None})
        post=[d for d in prior if joined[d,"no-dex"]>=minimum and joined[d,"dex"]>=minimum]
        result["post_eligible_paired_donors"]=len(post)
        result["minimum_recall"]=min(x["recall"] for x in result["libraries"])
        result["below_90pct"]=[x["gsm"] for x in result["libraries"] if x["recall"]<.9]
        result["below_95pct"]=[x["gsm"] for x in result["libraries"] if x["recall"]<.95]
        result["status"]="GO_BEAMA_PREFLIGHT" if not result["below_90pct"] and len(post)>=50 else "HOLD_BEAMA_PREFLIGHT"
    return result
def main():
    p=argparse.ArgumentParser();p.add_argument("--output",required=True)
    args=p.parse_args()
    result={"status":"HOLD"}
    try:result=check()
    except Exception as e:result["error"]=type(e).__name__+": "+str(e)[:200]
    Path(args.output).write_text(json.dumps(result,indent=2)+"\n")
    print("P4_PREFLIGHT",result["status"],"prior",result.get("prior_eligible_paired_donors"),
          "post",result.get("post_eligible_paired_donors"),"below90",result.get("below_90pct"),
          "below95",result.get("below_95pct"),"error",result.get("error"))
    if result["status"]!="GO_BEAMA_PREFLIGHT":raise SystemExit(2)
if __name__=="__main__":main()
