#!/usr/bin/env python3
"""Read-only CKD top3 source-level effect allele and reference LD consistency.

Checks Stage2B harmonized GWAS/pQTL rows, Stage2C subset, metadata
effect/reference allele and sign derivation. Does not inspect raw original GWAS
nor infer whether cohort-matched LD would agree. No sensitive records printed.
"""
import argparse,csv,gzip,json,math
from pathlib import Path
GENES=("SDCCAG8","INHBC","UMOD")
def read(p):
    opener=gzip.open if p.suffix==".gz" else open
    with opener(p,"rt",newline="",encoding="utf-8-sig") as f:
        return list(csv.DictReader(f,delimiter="\t"))
def audit(g,b2,b3,ld):
    a,b,m=(read(b2/(g+".tsv.gz")),read(b3/(g+".tsv.gz")),read(ld/(g+".ld.meta.tsv")))
    if not a or not b or not m:raise ValueError(g+" empty inputs")
    a={r["snp"]:r for r in a};b={r["snp"]:r for r in b}
    mismatch={"position":0,"effect_allele":0,"reference_pair":0,"reference_sign":0,"beta_pqtl":0,"beta_outcome_aligned":0}
    checked=0;missing_in_b2=0
    for row in m:
        k=row["snp"];r=a.get(k);s=b.get(k)
        if r is None or s is None:missing_in_b2+=1;continue
        checked+=1
        if r["pos37"]!=row["pos37"] or s["pos37"]!=row["pos37"]:mismatch["position"]+=1
        ea=r["allele1_pqtl"].upper();other=r["allele0_pqtl"].upper()
        if ea!=row["effect_allele"].upper() or other!=row["other_allele"].upper():mismatch["effect_allele"]+=1
        ref=row["reference_ref"].upper();alt=row["reference_alt"].upper()
        if {ea,other}!={ref,alt}:mismatch["reference_pair"]+=1
        altfreq=float(row["reference_alt_freq"])
        if not 0<=altfreq<=1:raise ValueError(g+" bad reference frequency")
        major=alt if altfreq>0.5 else ref
        desired=1 if ea==major else -1
        if int(float(row["effect_vs_ld_major_sign"]))!=desired:mismatch["reference_sign"]+=1
        for label in ("beta_pqtl","beta_outcome_aligned"):
            v=float(r[label]);w=float(s[label])
            if not (math.isfinite(v) and math.isfinite(w) and abs(v-w)<=1e-10*max(1,abs(v),abs(w))):
                mismatch[label]+=1
    result={"gene":g,"ld_rows":len(m),"matched_rows":checked,"missing_rows":missing_in_b2,
        "mismatch":mismatch,"source_harmonization_consistency_only":True,
        "raw_gwas_effect_allele_independently_verified":False,
        "reference_ld_cohort_compatibility_verified":False}
    result["status"]="SOURCE_CONSISTENCY_PASS_RAW_UNVERIFIED" if checked==len(m) and not any(mismatch.values()) else "SOURCE_CONSISTENCY_FAIL"
    return result
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--stage2b",required=True,type=Path);p.add_argument("--stage2c",required=True,type=Path);p.add_argument("--lddir",required=True,type=Path)
    x=p.parse_args();results=[audit(g,x.stage2b,x.stage2c,x.lddir) for g in GENES]
    print(json.dumps({"audit":"CKD_TOP3_EFFECT_SIGN_READONLY","results":results},indent=2))
    raise SystemExit(0 if all(r["status"]=="SOURCE_CONSISTENCY_PASS_RAW_UNVERIFIED" for r in results) else 2)
