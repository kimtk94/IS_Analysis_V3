#!/usr/bin/env python3
"""Audit only high-MAF-discordance variants; aggregate, existing data, read-only."""
import argparse,csv,gzip,json,math
from pathlib import Path
GENES=("SDCCAG8","INHBC","UMOD")
def read(p):
 op=gzip.open if p.suffix==".gz" else open
 with op(p,"rt",newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def main():
 p=argparse.ArgumentParser();p.add_argument("--input",type=Path,required=True);p.add_argument("--lddir",type=Path,required=True);a=p.parse_args()
 out=[]
 for g in GENES:
  s={r["snp"]:r for r in read(a.input/(g+".tsv.gz"))}
  m=read(a.lddir/(g+".ld.meta.tsv"))
  for r in m:
   q=s.get(r["snp"])
   if not q:continue
   try:
    pm=float(q["maf_pqtl"]);am=float(r["reference_alt_freq"])
    if not all(math.isfinite(z) for z in (pm,am)):continue
   except (KeyError,TypeError,ValueError):continue
   d=abs(pm-min(am,1-am))
   if d<=.1:continue
   ea=r["effect_allele"].upper();other=r["other_allele"].upper();ref=r["reference_ref"].upper();alt=r["reference_alt"].upper()
   out.append({"gene":g,"snp":r["snp"],"pos37":r["pos37"],"effect_allele":ea,"other_allele":other,
     "reference_ref":ref,"reference_alt":alt,"reference_alt_freq":am,"pqtl_maf":pm,
     "absolute_maf_difference":d,"reference_alleles_match":{ea,other}=={ref,alt},
     "raw_consortium_verified":False})
 print(json.dumps({"screen":"MAF_OUTLIERS_GT_0_1","outlier_count":len(out),"outliers":out},indent=2))
 return 0 if len(out)==4 else 2
if __name__=="__main__":raise SystemExit(main())
