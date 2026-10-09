#!/usr/bin/env python3
"""EAS anchor LD proxy audit using PLINK2 variant-major dosage export.

Reference-only r²; not independent signals, statistical fine-mapping, or proof of causality.
"""
import argparse,csv,json,subprocess
from pathlib import Path
import numpy as np
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1")
REF=Path("/srv/is-analysis/data/is/ld_reference/1kg_eas_grch37/pgen_gt")
LOCUS={"IS_XDATA_G0007":"L001","IS_XDATA_G0015":"L002",
       "IS_XDATA_G0022":"L003","IS_XDATA_G0023":"L004"}
def read(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def correlation(x,y):
    mask=np.isfinite(x)&np.isfinite(y)
    if mask.sum()<50:return None,int(mask.sum())
    a=x[mask];b=y[mask];da=a-a.mean();db=b-b.mean()
    den=np.sqrt(np.dot(da,da)*np.dot(db,db))
    if den==0:return None,int(mask.sum())
    r=float(np.dot(da,db)/den)
    return max(0.,min(1.,r*r)),int(mask.sum())
def genotypes(path):
    with path.open() as f:
        h=next(f).rstrip("\n").split("\t")
        n=len(h)-6
        for line in f:
            p=line.rstrip("\n").split("\t",6)
            if len(p)!=7:continue
            values=np.fromstring(p[6].replace("NA","nan"),sep="\t")
            if len(values)!=n:raise ValueError("Unexpected dosage sample count at "+p[1])
            yield p[:6],values
def calc(traw,lead_id,output):
    lead=None
    for fields,dose in genotypes(traw):
        if fields[1]==lead_id:
            lead=dose;break
    if lead is None:raise ValueError("Lead variant not exported: "+lead_id)
    records=[]
    for fields,dose in genotypes(traw):
        if fields[1]==lead_id:continue
        r2,n=correlation(lead,dose)
        if r2 is not None and r2>=0.2:
            records.append({"lead_variant_reference":lead_id,"proxy_variant":fields[1],
                            "chrom":fields[0],"pos":fields[3],"ref_counted_allele":fields[4],
                            "other_allele":fields[5],"r2":round(r2,8),"n_pairwise":n})
    with output.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["lead_variant_reference","proxy_variant","chrom",
            "pos","ref_counted_allele","other_allele","r2","n_pairwise"],delimiter="\t")
        w.writeheader();w.writerows(sorted(records,key=lambda x:-x["r2"]))
    return {"lead":lead_id,"n_samples":len(lead),"lead_nonmissing":int(np.isfinite(lead).sum()),
            "n_proxy_r2_ge_0.2":len(records),
            "n_proxy_r2_ge_0.8":sum(x["r2"]>=0.8 for x in records),
            "max_proxy_r2":max((x["r2"] for x in records),default=None)}
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",type=Path,default=ROOT)
    ap.add_argument("--reference",type=Path,default=REF)
    ap.add_argument("--plink2",default="plink2")
    ap.add_argument("--window-kb",type=int,default=500)
    a=ap.parse_args()
    audit=read(a.root/"LEAD_ALLELE_REFERENCE_AUDIT.tsv")
    out=a.root/"anchor_ld_probe"
    out.mkdir(exist_ok=True,parents=True)
    results=[]
    for r in audit:
        group=r["group_id"]
        if group not in LOCUS:continue
        tag=LOCUS[group]
        prefix=a.reference/f"BBJ_IS_{tag}.1KG_EAS.GRCh37"
        chrom,pos,allele1,allele2=r["lead_variant"].split(":")
        lead_id=None
        with Path(str(prefix)+".pvar").open() as f:
            for line in f:
                if line.startswith("#"):continue
                row=line.rstrip("\n").split("\t")
                if row[0]==chrom and row[1]==pos and set([row[3],row[4]])==set([allele1,allele2]):
                    lead_id=row[2]
                    break
        if lead_id is None:
            raise RuntimeError("Lead allele pair absent from ref PVAR: "+r["lead_variant"])
        export=out/f"{group}_{tag}_AV"
        traw=Path(str(export)+".traw")
        cmd=[a.plink2,"--pfile",str(prefix),"--chr",chrom,
             "--from-bp",str(max(1,int(pos)-1000*a.window_kb)),
             "--to-bp",str(int(pos)+1000*a.window_kb),
             "--export","Av","--out",str(export)]
        if not traw.is_file():
            proc=subprocess.run(cmd,capture_output=True,text=True,timeout=180)
            if proc.returncode!=0:raise RuntimeError(proc.stderr+"\n"+proc.stdout[-2000:])
        report=calc(traw,lead_id,out/f"{group}_LEAD_R2_PROXIES.tsv")
        report.update({"group_id":group,"legacy_locus":tag,
                       "orientation":"REVERSED" if r["status"].startswith("SWAPPED") else "EXACT",
                       "window_kb":a.window_kb,
                       "status":"REFERENCE_R2_ONLY_NOT_FINE_MAPPING"})
        results.append(report)
        print(json.dumps(report),flush=True)
    summary={"completed_anchor_groups":len(results),"n_expected":len(LOCUS),
        "reference_population":"1KG EAS; genotype n=504",
        "per_group":results,
        "limitation":"r2 proxies estimated against 1KG EAS only; no genome-wide LD clumping; other 26 regions lack reference"}
    (out/"ANCHOR_LD_PROXY_SUMMARY.json").write_text(json.dumps(summary,indent=2))
    return summary
if __name__=="__main__":main()
