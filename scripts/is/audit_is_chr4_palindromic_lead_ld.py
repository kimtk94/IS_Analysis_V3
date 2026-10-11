#!/usr/bin/env python3
"""Genomic ALT-dose correlation at chr4 palindrome and BBJ lead; strand cautious."""
import csv,json,math
from pathlib import Path
import pandas as pd,numpy as np
ROOT=Path("/srv/is-analysis/results/is/audits")
O=ROOT/"is_chr4_palindrome_lead_ld_20261011_v1"
P="/srv/is-analysis/data/is/ld_reference/1kg_eas_grch37/integrated_call_samples_v3.20130502.ALL.panel"
CH="4:81182554:T:C"
rows=[]
panel=pd.read_csv(P,sep="\t")
jpt=set(panel.loc[panel["pop"]=="JPT","sample"])
for anc,root,name in [
 ("EAS504",ROOT/"is_priority4_eur503_eas504_conditional_20261011_v1","BBJ_IS_L001.EAS504.fullcis.raw"),
 ("EUR503",ROOT/"is_priority4_eur503_fullcis_signed_ld_qc_20261011_v1","BBJ_IS_L001.EUR503.raw")]:
 f=pd.read_csv(root/name,sep="\t")
 col=next(c for c in f if c.startswith(CH+"_"))
 lead=f[col].to_numpy(float)
 if col.rsplit("_",1)[-1]=="T":lead=2-lead
 calls={}
 with (O/(anc+"_4_81181072_GT.tsv")).open() as x:
  for line in x:
   sample,gt=line.rstrip().split("\t")
   alleles=gt.replace("|","/").split("/")
   calls[sample]=float(sum(int(a) for a in alleles)) if len(alleles)==2 and all(a in ("0","1") for a in alleles) else np.nan
 alt=np.array([calls[i] for i in f.IID.astype(str)],float)
 if not np.isfinite(alt).all():raise ValueError("Missing palindrome genotype")
 populations={anc:np.ones(len(f),bool)}
 if anc=="EAS504":
  populations={"EAS504":np.ones(len(f),bool),"JPT104":f.IID.isin(jpt).to_numpy(),"otherEAS400":~f.IID.isin(jpt).to_numpy()}
 for pop,mask in populations.items():
  xx=lead[mask];yy=alt[mask]
  rr=float(np.corrcoef(xx,yy)[0,1])
  rows.append(dict(population=pop,N=len(xx),
   BBJ_lead_ALT_C_EAF=float(xx.mean()/2),palindrome_genomic_ALT_T_EAF=float(yy.mean()/2),
   signed_r_ALT_C_vs_ALT_T=rr,r2=rr*rr,
   reference_alleles="4:81181072:A:T",GIGA_variant_ID="4:81181072:T:A",
   strand_provenance_UNRESOLVED=True,
   no_conditional_GWAS_or_causal_inference=True))
with (O/"IS_CHR4_PALINDROME_REFERENCE_LD.tsv").open("w") as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
(O/"IS_CHR4_PALINDROME_REFERENCE_LD.json").write_text(json.dumps(rows,indent=2)+"\n")
for a in rows:print(a)
