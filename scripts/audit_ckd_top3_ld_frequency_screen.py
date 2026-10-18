#!/usr/bin/env python3
"""Read-only CKD reference frequency compatibility *screen*, not proof of LD suitability."""
import argparse,csv,gzip,json,math
from pathlib import Path
GENES=("SDCCAG8","INHBC","UMOD")
def read(p):
    op=gzip.open if p.suffix==".gz" else open
    with op(p,"rt",newline="",encoding="utf-8-sig") as f:return list(csv.DictReader(f,delimiter="\t"))
def screen(g,summary,ldroot):
    A={r["snp"]:r for r in read(summary/(g+".tsv.gz"))}
    M=read(ldroot/(g+".ld.meta.tsv"))
    diffs=[];unavailable=0;bad=0
    for m in M:
        r=A.get(m["snp"])
        if not r:bad+=1;continue
        try:
            maf=float(r["maf_pqtl"])
            alt=float(m["reference_alt_freq"])
            if not all(map(math.isfinite,(maf,alt))) or not(0<=maf<=.5 and 0<=alt<=1):raise ValueError()
            refmaf=min(alt,1-alt)
            diffs.append(abs(maf-refmaf))
        except (KeyError,ValueError,TypeError):unavailable+=1
    diffs.sort()
    def quant(q):return diffs[min(len(diffs)-1,int(q*(len(diffs)-1)))] if diffs else None
    return {"gene":g,"ld_variants":len(M),"maf_compared":len(diffs),"maf_missing":unavailable,
            "unmatched":bad,"abs_maf_delta_median":quant(.5),"abs_maf_delta_p95":quant(.95),
            "abs_maf_delta_max":quant(1),"fraction_maf_delta_gt_0_1":sum(x>.1 for x in diffs)/len(diffs) if diffs else None,
            "screen_only":True,"ancestry_cohort_ld_verified":False,
            "raw_consortrium_effect_alleles_verified":False,
            "status":"FREQUENCY_SCREEN_AVAILABLE" if diffs and bad==0 else "INSUFFICIENT_INPUT"}
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--stage2c",type=Path,required=True);p.add_argument("--lddir",type=Path,required=True);a=p.parse_args()
    results=[screen(g,a.stage2c,a.lddir) for g in GENES]
    print(json.dumps({"kind":"CKD_REFERENCE_MAF_COMPATIBILITY_SCREEN","results":results},indent=2))
    raise SystemExit(0 if all(x["status"]=="FREQUENCY_SCREEN_AVAILABLE" for x in results) else 2)
