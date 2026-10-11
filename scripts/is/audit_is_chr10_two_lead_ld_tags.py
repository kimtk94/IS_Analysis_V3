#!/usr/bin/env python3
"""Descriptive signed LD to two chr10 stroke leads, not fine-mapping credible sets."""
import csv,json,math
from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path("/srv/is-analysis/results/is/audits")
source=ROOT/"is_priority4_fullcis_bbj_giga_snp_exact_20261011_v1/IS_PRIORITY4_FULLCIS_EBI_BBJ_GIGASTROKE_SNP_EVIDENCE.tsv"
raw=ROOT/"is_priority4_eur503_eas504_conditional_20261011_v1/BBJ_IS_L002.EAS504.fullcis.raw"
panel="/srv/is-analysis/data/is/ld_reference/1kg_eas_grch37/integrated_call_samples_v3.20130502.ALL.panel"
out=ROOT/"is_chr10_two_lead_ld_tags_20261011_v1"
leads=("10:105459834:T:C","10:105599770:T:G")
r=pd.read_csv(raw,sep="\t");x=r.iloc[:,6:].to_numpy(float,copy=True);keys=[]
for j,c in enumerate(r.columns[6:]):
 k,a=c.rsplit("_",1);v=k.split(":")
 if a==v[2]:x[:,j]=2-x[:,j]
 elif a!=v[3]:raise ValueError("Allele count mismatch")
 keys.append(k)
idx={k:i for i,k in enumerate(keys)}
if not all(k in idx for k in leads):raise ValueError("Missing leads")
meta=pd.read_csv(panel,sep="\t")
jpt=set(meta.loc[meta["pop"]=="JPT","sample"])
pop={"EAS504":np.ones(len(r),bool),"JPT104":r.IID.isin(jpt).to_numpy(),"otherEAS400":(~r.IID.isin(jpt)).to_numpy()}
records=[]
for group,sel in pop.items():
 z=x[sel];z=z-z.mean(0);ss=np.sqrt((z*z).sum(0))
 for label,k in zip(("BBJ","GIGASTROKE"),leads):
  lead=z[:,idx[k]]
  vals=(z.T@lead)/(ss*ss[idx[k]])
  for snp,v in zip(keys,vals):
   if not np.isfinite(v):continue
   records.append(dict(population=group,lead=label,lead_b37=k,variant_b37=snp,
    signed_r=float(v),r2=float(v*v),r2_ge_0p8=bool(v*v>=.8),r2_ge_0p2=bool(v*v>=.2),
    not_credible_set=True))
# Per genomic SNP: annotate EBI molecular QTL significance without double counting.
seen={}
with source.open() as f:
 for t in csv.DictReader(f,delimiter="\t"):
  if t["locus"]!="BBJ_IS_L002":continue
  k=t["SNP_GRCh37"];gene=t["symbol"]
  if k not in idx:continue
  seen.setdefault(k,{})[gene]=float(t["QTD000161_nominal_p"])
for t in records:
 q=seen.get(t["variant_b37"],{})
 t["CALHM2_QTD_p"]=q.get("CALHM2","")
 t["NEURL_QTD_p"]=q.get("NEURL","")
out.mkdir(parents=True,exist_ok=True)
with (out/"IS_CHR10_EAS_JPT_TWO_LEAD_PER_SNP_LD.tsv").open("w") as f:
 w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter="\t");w.writeheader();w.writerows(records)
summary=[]
for group in pop:
 for label in ("BBJ","GIGASTROKE"):
  a=[t for t in records if t["population"]==group and t["lead"]==label]
  summary.append(dict(population=group,lead=label,evaluable=len(a),
   n_r2_ge_0p8=sum(t["r2_ge_0p8"] for t in a),
   n_r2_ge_0p2=sum(t["r2_ge_0p2"] for t in a),
   n_QTD_CALHM2_nominal_p_lt_1e4_and_r2_ge_0p8=sum(t["r2_ge_0p8"] and isinstance(t["CALHM2_QTD_p"],float) and t["CALHM2_QTD_p"]<1e-4 for t in a),
   not_independent_signal_or_credible_set=True))
with (out/"IS_CHR10_EAS_JPT_TWO_LEAD_LD_QC.tsv").open("w") as f:
 w=csv.DictWriter(f,fieldnames=list(summary[0]),delimiter="\t");w.writeheader();w.writerows(summary)
(out/"IS_CHR10_EAS_JPT_TWO_LEAD_LD_MANIFEST.json").write_text(json.dumps(dict(
 result="EXPLORATORY_REFERENCE_LD_TAGS_NOT_SUSIE_CREDIBLE_SETS",summary=summary,
 source_variants=len(keys),existing_2225_genes_unmodified=True),indent=2)+"\n")
print(json.dumps(summary,indent=2))
