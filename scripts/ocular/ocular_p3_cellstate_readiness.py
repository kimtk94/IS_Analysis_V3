#!/usr/bin/env python3
"""Aggregate-only cell-state pseudobulk readiness audit for GSE313736."""
import argparse,csv,gzip,json,tempfile,collections
from pathlib import Path
from datetime import datetime,timezone
from ocular_p1_integrity import BASE,download

ANNOTATIONS=[
 ("celltypist_coarse","celltypist_chen2025_coarse","celltypist_chen2025_coarse_score"),
 ("celltypist_fine","celltypist_chen2025_fine","celltypist_chen2025_fine_score"),
 ("scpred","scpred_vanzyl2022","scpred_vanzyl2022_score"),
]
def score(x):
    try:return float(x)
    except (ValueError,TypeError):return -1.
def audit(path):
    allrows=collections.Counter();groups=collections.defaultdict(lambda:collections.Counter())
    labels=collections.defaultdict(lambda:collections.Counter());weak=collections.Counter()
    donors=set();conditions=set()
    with gzip.open(path,"rt",encoding="utf8",errors="replace") as f:
        reader=csv.DictReader(f,delimiter="\t")
        required={"donor_id_demuxlet","DEX Treatment","barcode","droplet_type"}
        for _,lbl,conf in ANNOTATIONS: required.update((lbl,conf))
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("MISSING_ANNOTATION_COLUMNS: "+str(sorted(required-set(reader.fieldnames or []))))
        for r in reader:
            donor=r["donor_id_demuxlet"].strip()
            condition=r["DEX Treatment"].strip().lower()
            if condition not in ("dex","no-dex"):raise ValueError("UNKNOWN_TREATMENT")
            if r["droplet_type"]!="SNG":continue
            conditions.add(condition);donors.add(donor);allrows[condition]+=1
            for kind,lbl,conf in ANNOTATIONS:
                label=r[lbl].strip()
                value=score(r[conf])
                labels[(kind,"all")][label]+=1
                if value<.5:
                    weak[kind]+=1
                    continue
                for threshold in (.5,.75):
                    if value>=threshold:
                        groups[(kind,threshold,label,donor,condition)]["cells"]+=1
    summary=[]
    for kind,lbl,conf in ANNOTATIONS:
        top=labels[(kind,"all")].most_common(40)
        for label,ncells in top:
            if not label or label.lower() in ("na","nan","unknown"):continue
            x={"method":kind,"cell_state":label,"annotated_cells_total":ncells}
            for threshold in (.5,.75):
                donors_by_cond={}
                ncell=0
                for d in donors:
                    c=groups[(kind,threshold,label,d,"no-dex")]["cells"]
                    e=groups[(kind,threshold,label,d,"dex")]["cells"]
                    ncell+=c+e
                    for minimum in (20,30,50):
                        if c>=minimum and e>=minimum:
                            donors_by_cond[minimum]=donors_by_cond.get(minimum,0)+1
                x[f"score_ge_{threshold}_cells"]=ncell
                x[f"score_ge_{threshold}_paired_donors_min20"]=donors_by_cond.get(20,0)
                x[f"score_ge_{threshold}_paired_donors_min30"]=donors_by_cond.get(30,0)
                x[f"score_ge_{threshold}_paired_donors_min50"]=donors_by_cond.get(50,0)
            x["primary_model_candidate"]=x["score_ge_0.5_paired_donors_min30"]>=30
            x["exploratory_model_candidate"]=x["score_ge_0.5_paired_donors_min20"]>=15
            summary.append(x)
    return {"dataset":"GSE313736","run_utc":datetime.now(timezone.utc).isoformat(),
        "status":"PASS_METADATA_CELL_STATE_COUNTS" if len(donors)==59 and conditions=={"dex","no-dex"} else "HOLD",
        "conditions":dict(allrows),"unique_donor_count":len(donors),
        "cell_label_count":len(summary),"unconfident_annotation_rows_by_method":dict(weak),
        "cell_states":summary,
        "note":"Paired donor counts from public metadata are necessary but not sufficient for pseudobulk raw-count QC. Cell labels are model predictions and require biological validation."}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True);p.add_argument("--metadata",default="")
    args=p.parse_args()
    if args.metadata:
        data=audit(args.metadata)
    else:
        with tempfile.TemporaryDirectory(prefix="ocular-cellstate-") as t:
            loc=Path(t)/"metadata.tsv.gz"
            download(BASE+"series/GSE313nnn/GSE313736/suppl/GSE313736_HTMC_scRNA_cell_metadata.tsv.gz",loc,12*1024*1024)
            data=audit(loc)
    dst=Path(args.output);dst.parent.mkdir(parents=True,exist_ok=True)
    dst.write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n")
    print("QC",data["status"],"donors",data["unique_donor_count"],"states",data["cell_label_count"])
    for a in data["cell_states"]:
        if a["primary_model_candidate"] or a["exploratory_model_candidate"]:
            print("MODEL_CANDIDATE",a["method"],a["cell_state"],
                "primary",a["score_ge_0.5_paired_donors_min30"],
                "exploratory",a["score_ge_0.5_paired_donors_min20"])
    if data["status"]!="PASS_METADATA_CELL_STATE_COUNTS":raise SystemExit(2)
if __name__=="__main__":main()
