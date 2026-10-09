#!/usr/bin/env python3
"""Calculate lead-centered reference r² for NEW GWS regions only.
Only 1000G EAS, not SuSiE or independence proof.
"""
import argparse,csv,json,subprocess
from pathlib import Path
from build_anchor_ld_proxies import calc
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1")
REF=Path("/srv/is-analysis/data/is/ld_reference/broad_v1/1kg_eas")
NEW=["IS_XDATA_G0004","IS_XDATA_G0008","IS_XDATA_G0021"]
def run(root,ref,group):
    with (root/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv").open() as f:
        rows={r["group_id"]:r for r in csv.DictReader(f,delimiter="\t")}
    if group not in NEW:raise ValueError("Only three novel GWS groups are permitted")
    g=rows[group]
    prefix=ref/(group+".stable")
    triplet=[Path(str(prefix)+x) for x in (".pgen",".pvar",".psam")]
    if not all(x.is_file() for x in triplet):
        raise FileNotFoundError("Complete 504-person regional PGEN is required")
    study=g["lead_variant"].split(":")
    lead_id=None;orientation=None
    with triplet[1].open() as f:
        for line in f:
            if line.startswith("#"):continue
            r=line.rstrip("\n").split("\t")
            if r[0]==study[0] and r[1]==study[1] and set(r[3:5])==set(study[2:4]):
                lead_id=r[2]
                orientation="EXACT" if r[3:5]==study[2:4] else "REVERSED_ALLELE_ORDER"
                break
    if not lead_id or lead_id==".":
        raise ValueError("Lead must have stable chr:pos:REF:ALT ID: "+group)
    out=root/"new_gws_ld_proxy"
    out.mkdir(parents=True,exist_ok=True)
    export=out/(group+"_STABLE_DOSAGE")
    traw=Path(str(export)+".traw")
    if not traw.is_file():
        cmd=["plink2","--pfile",str(prefix),"--export","Av","--out",str(export)]
        p=subprocess.run(cmd,capture_output=True,text=True,timeout=360)
        if p.returncode:raise RuntimeError(p.stderr+"\n"+p.stdout[-1500:])
    result=calc(traw,lead_id,out/(group+"_LEAD_R2_PROXIES.tsv"))
    result.update({"group_id":group,"lead_allele_order":orientation,
        "claim":"Reference r² only, not LD-independent or causal"})
    (out/(group+"_R2_QC.json")).write_text(json.dumps(result,indent=2))
    return result
if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--root",type=Path,default=ROOT)
    parser.add_argument("--ref",type=Path,default=REF)
    parser.add_argument("--group",choices=NEW,required=True)
    args=parser.parse_args()
    print(json.dumps(run(args.root,args.ref,args.group),indent=2))
