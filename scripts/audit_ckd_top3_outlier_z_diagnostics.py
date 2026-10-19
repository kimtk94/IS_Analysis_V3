#!/usr/bin/env python3
"""Read-only SNP-level Z and MAF anomaly diagnostics, not a SuSiE convergence test."""
import csv,gzip,json,math,statistics
from pathlib import Path
B=Path("/srv/is-analysis/results/ckd/stage2c_susie/susie_input")
L=Path("/srv/is-analysis/data/ckd/stage2c_ld/ld")
TARGET={"SDCCAG8":{"rs35646972","rs6663222"},"INHBC":{"rs11834807"},"UMOD":{"rs11860566"}}
output=[]
for gene,snps in TARGET.items():
 with gzip.open(B/(gene+".tsv.gz"),"rt") as f:
  rr=list(csv.DictReader(f,delimiter="\t"))
 with (L/(gene+".ld.meta.tsv")).open() as f:
  mm={r["snp"]:r for r in csv.DictReader(f,delimiter="\t")}
 scores=[]
 badse=0
 for r in rr:
  try:
   bp=float(r["beta_pqtl"]);bg=float(r["beta_outcome_aligned"])
   # Stage2C source may use varied SE header spelling; fail closed if unavailable.
   sp=next(float(r[k]) for k in ("se_pqtl","se_pqtl_raw") if r.get(k))
   sg=next(float(r[k]) for k in ("se_outcome_aligned","se_outcome") if r.get(k))
   if sp<=0 or sg<=0 or not all(map(math.isfinite,(bp,bg,sp,sg))):raise ValueError()
   scores.append((r["snp"],bp/sp,bg/sg,r))
  except (ValueError,KeyError,StopIteration):badse+=1
 aligned={s:(zp,zg,r) for s,zp,zg,r in scores}
 events=[]
 for snp in sorted(snps):
  if snp not in aligned:
   events.append({"snp":snp,"status":"NO_VALID_Z_IN_STAGE2C"})
   continue
  zp,zg,r=aligned[snp];m=mm.get(snp)
  events.append({"snp":snp,"status":"AVAILABLE","z_pqtl":zp,"z_egfr_aligned":zg,
   "abs_z_pqtl_rank_fraction":sum(abs(z)<=abs(zp) for _,z,_,_ in scores)/len(scores),
   "abs_z_egfr_rank_fraction":sum(abs(z)<=abs(zg) for _,_,z,_ in scores)/len(scores),
   "source_maf_pqtl":r.get("maf_pqtl"),
   "ld_reference_alt_freq":m.get("reference_alt_freq") if m else None})
 output.append({"gene":gene,"stage2c_rows":len(rr),"usable_z_pairs":len(scores),"excluded_se_or_z":badse,
  "events":events,"status":"Z_SUMMARY_SCREEN_ONLY"})
print(json.dumps({"audit":"CKD_TOP3_OUTLIER_Z_DIAGNOSTICS","results":output,
 "not_a_ld_z_consistency_test":True,"not_a_susie_causal_diagnosis":True},indent=2))
