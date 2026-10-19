#!/usr/bin/env python3
"""Read-only UMOD signed-LD versus Z neighborhood descriptive diagnostics."""
import csv,gzip,json,numpy as np
from pathlib import Path
base=Path("/srv/is-analysis")
ld=base/"data/ckd/stage2c_ld/ld"
with gzip.open(base/"results/ckd/stage2c_susie/susie_input/UMOD.tsv.gz","rt") as f:
 rows={r["snp"]:r for r in csv.DictReader(f,delimiter="\t")}
with (ld/"UMOD.ld.meta.tsv").open() as f:
 meta=list(csv.DictReader(f,delimiter="\t"))
n=len(meta)
sign=np.array([float(r["effect_vs_ld_major_sign"]) for r in meta])
mat=np.memmap(ld/"UMOD.ld.bin",dtype="<f4",mode="r",shape=(n,n))
z=np.array([float(rows[m["snp"]]["beta_outcome_aligned"])/float(rows[m["snp"]]["se_outcome"]) for m in meta])
zp=np.array([float(rows[m["snp"]]["beta_pqtl"])/float(rows[m["snp"]]["se_pqtl"]) for m in meta])
targets=["rs11860566",meta[int(np.argmax(abs(z)))]["snp"],meta[int(np.argmax(abs(zp)))]["snp"]]
out=[]
for target in dict.fromkeys(targets):
 idx=next(i for i,m in enumerate(meta) if m["snp"]==target)
 r=np.asarray(mat[idx,:],dtype=float)*sign[idx]*sign
 tagged=np.where((abs(r)>=0.8)&(np.arange(n)!=idx))[0]
 diag={"snp":target,"gw_z":float(z[idx]),"pqtl_z":float(zp[idx]),"r_ge_0_8_count":int(len(tagged))}
 for key,v in (("gwas",z),("pqtl",zp)):
  if len(tagged):
   pred=r[tagged]*v[idx]
   residual=v[tagged]-pred
   diag[key+"_tagged_z_residual_rms"]=float(np.sqrt(np.mean(residual**2)))
   diag[key+"_tagged_z_signed_correlation"]=float(np.corrcoef(v[tagged],r[tagged])[0,1]) if np.std(v[tagged])>0 else None
  else:
   diag[key+"_tagged_z_residual_rms"]=None
 out.append(diag)
print(json.dumps({"audit":"UMOD_LD_Z_TAG_NEIGHBOR_SCREEN","ld_n":n,"targets":out,
 "limitation":"descriptive tag neighborhood; not a formal DENTIST-S or SuSiE convergence test"},indent=2))
