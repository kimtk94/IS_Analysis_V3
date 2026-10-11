#!/usr/bin/env python3
"""Exploratory 1000G EUR LD conditional GWAS tagging: alpha sensitivity only.

Stepwise lead-conditioned z: z_cond(j|lead)=(z_j-r_jL*z_L)/sqrt(1-r_jL²)
where z=beta/se and r from EUR503 ALT dosages. Never interpret as multivariable
GWAS conditioning, SuSiE, colocalization or causal effect. Shrinkage uses
r_alpha=(1-alpha)*r for off-diagonal correlations.
"""
import argparse,csv,json,math,hashlib
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
ALPHAS=(0.,0.01,0.05,0.1,0.2)
LEADS={"BBJ_IS_L001":"4:81182554:T:C","BBJ_IS_L002":"10:105459834:T:C"}
def rows(path):
 with path.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def run(source,genotypes,out):
 audit=rows(source)
 if len(audit)!=11050:raise ValueError("Full source drift")
 by_locus=defaultdict(dict)
 for r in audit:
  k=r["variant_GRCh37_REF_ALT"];loc=r["locus"]
  if k in by_locus[loc]:
   prev=by_locus[loc][k]
   if prev["BBJ_beta_ALT"]!=r["BBJ_beta_ALT"] or prev["BBJ_se"]!=r["BBJ_se"]:
    raise ValueError("GWAS z differs across genes sharing same SNP")
  by_locus[loc][k]=r
 final=[];summ=[]
 for loc,lead in LEADS.items():
  df=pd.read_csv(genotypes/f"{loc}.EUR503.raw",sep="\t")
  if len(df)!=503:raise ValueError("Ancestry reference sample count drift")
  m=df.iloc[:,6:].to_numpy(dtype=float,copy=True)
  keys=[]
  for i,col in enumerate(df.columns[6:]):
   key,allele=col.rsplit("_",1);ch,pos,ref,alt=key.split(":")
   if allele==ref:m[:,i]=2-m[:,i]
   elif allele!=alt:raise ValueError("Counted allele does not match SNP REF/ALT")
   keys.append(key)
  if lead not in keys:raise ValueError("Lead SNP absent")
  idx={k:i for i,k in enumerate(keys)}
  y=m[:,idx[lead]]
  if not np.isfinite(m).all():raise ValueError("Nonmissing reference required")
  Y=y-y.mean()
  sy=np.sqrt(np.sum(Y*Y))
  if sy<=0:raise ValueError("Monomorphic lead")
  X=m-m.mean(axis=0)
  norm=np.sqrt(np.sum(X*X,axis=0))
  rr=np.full(len(keys),np.nan)
  ok=norm>0
  rr[ok]=(Y@X[:,ok])/(sy*norm[ok])
  rr=np.clip(rr,-1,1)
  targets=by_locus[loc]
  selected=[(k,rr[idx[k]],r) for k,r in targets.items() if k in idx]
  if len(selected)<2500:raise ValueError("Incomplete GWAS-to-reference key selection")
  snpc=0
  for k,rho,x in selected:
   if not np.isfinite(rho):continue
   z=float(x["BBJ_beta_ALT"])/float(x["BBJ_se"])
   zl=float(targets[lead]["BBJ_beta_ALT"])/float(targets[lead]["BBJ_se"])
   if not np.isfinite(z) or not np.isfinite(zl):raise ValueError("Nonfinite GWAS z")
   ch,pos,*_=k.split(":");block=int(pos)//250000
   for a in ALPHAS:
    ra=(1-a)*rho
    den=1-ra*ra
    if k==lead or den<1e-9:
     zc=None;reason="LEAD_OR_PERFECT_LD"
    else:
     zc=(z-ra*zl)/math.sqrt(den)
     reason="REF_EUR_LEAD_CONDITIONED_DIAGNOSTIC"
    final.append(dict(locus=loc,variant_b37=k,block_250kb=block,BBJ_z=z,
       GWAS_lead_b37=lead,GWAS_lead_z=zl,EUR503_signed_r=rho,
       shrinkage_alpha=a,shrunken_signed_r=ra,
       lead_conditional_z=zc if zc is not None else "",
       abs_z_gte_5p45=(abs(zc)>=5.45 if zc is not None else ""),
       valid=zc is not None,status=reason))
   snpc+=1
  summ.append(dict(locus=loc,n_reference_samples=503,
    n_unique_GWAS_SNPs=len(targets),n_with_EUR503_ref=int(len(selected)),
    n_informative_pairwise_r=int(snpc),
    max_locus_GWAS_lead_abs_z=float(abs(float(targets[lead]["BBJ_beta_ALT"])/float(targets[lead]["BBJ_se"]))),
    scientific_gate="EUR_REFERENCE_LEAD_CONDITIONED_ASSOCIATION_SCREEN_NOT_CAUSAL_NOT_SUSIE"))
 out.mkdir(parents=True,exist_ok=True)
 for name,data in (("IS_PRIORITY4_EUR_REFERENCE_LEAD_CONDITIONAL_Z.tsv",final),
                   ("IS_PRIORITY4_EUR_REFERENCE_LEAD_CONDITIONAL_QC.tsv",summ)):
  with (out/name).open("w",newline="") as f:
   w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter="\t");w.writeheader();w.writerows(data)
 c=defaultdict(list)
 for r in final:
  if r["valid"]:c[(r["locus"],r["shrinkage_alpha"])].append(abs(r["lead_conditional_z"]))
 stats=[]
 for (loc,a),v in sorted(c.items()):
  stats.append(dict(locus=loc,alpha=a,n_evaluable=int(len(v)),n_absZ_gt_5p45=int(sum(z>=5.45 for z in v)),
    max_absZ=max(v),n_absZ_gt_3=int(sum(z>=3 for z in v)),
    no_independence_inference=True))
 with (out/"IS_PRIORITY4_EUR_REF_LEAD_CONDITIONAL_ALPHA_SUMMARY.tsv").open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(stats[0]),delimiter="\t");w.writeheader();w.writerows(stats)
 manifest=dict(science_status="EXPLORATORY_REFERENCE_LD_PAIRWISE_LEAD_CONDITION_ONLY",
  n_genotype_individuals=503,n_loci=2,number_of_alpha_settings=len(ALPHAS),
  qc=summ,alpha_grid=stats,no_GTex_in_study_LD=True,
  no_conditional_GWAS_summary_calibration=True,no_SuSiE_coloc_performed=True,
  candidate_genes_preserved=2225)
 (out/"IS_PRIORITY4_EUR_REF_LEAD_CONDITIONAL_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")
 return manifest
if __name__=="__main__":
 p=argparse.ArgumentParser()
 for name in ("source","genotypes","out"):p.add_argument("--"+name,type=Path,required=True)
 a=p.parse_args();print(json.dumps(run(a.source,a.genotypes,a.out),indent=2))
