#!/usr/bin/env python3
"""JPT104, other-EAS400, EAS504, EUR503: BBJ reference LD sensitivity only."""
import argparse,csv,json
from pathlib import Path
import numpy as np,pandas as pd
LEADS={"BBJ_IS_L001":"4:81182554:T:C","BBJ_IS_L002":"10:105459834:T:C"}
CAND="10:105694487:G:T"
ALPHAS=(0,.01,.05,.1,.2)
def load(path):
 df=pd.read_csv(path,sep="\t");X=df.iloc[:,6:].to_numpy(float,copy=True);keys=[]
 for i,col in enumerate(df.columns[6:]):
  key,a=col.rsplit("_",1);v=key.split(":")
  if len(v)!=4:raise ValueError("Bad variant ID")
  if a==v[2]:X[:,i]=2-X[:,i]
  elif a!=v[3]:raise ValueError("Bad counted allele")
  keys.append(key)
 if not np.isfinite(X).all():raise ValueError("Missing genotypes")
 if len(set(keys))!=len(keys):raise ValueError("Duplicate variants")
 return df.IID.astype(str).to_numpy(),keys,X
def corr(X,li):
 xx=X-X.mean(0);norm=np.sqrt((xx*xx).sum(0));r=np.full(X.shape[1],np.nan)
 ok=norm>0;r[ok]=(xx[:,ok].T@xx[:,li])/(norm[ok]*norm[li]);return np.clip(r,-1,1)
def save(out,name,rows):
 with (out/name).open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
def boot(x,y,z,zlead,rng,n=2000):
 out=np.zeros((n,2))
 for j in range(n):
  k=rng.integers(0,len(x),len(x));a=x[k];b=y[k]
  rr=np.corrcoef(a,b)[0,1]
  if not np.isfinite(rr):raise ValueError("Bootstrap monomorphic pair")
  q=.95*rr;out[j]=[rr,(z-q*zlead)/np.sqrt(1-q*q)]
 return out
def run(source,eas,eur,panel,out):
 pop=pd.read_csv(panel,sep="\t")[["sample","pop","super_pop"]]
 jpt=set(pop.loc[pop["pop"]=="JPT","sample"])
 if len(jpt)!=104:raise ValueError("JPT population size drift")
 G={}
 with source.open() as f:
  for r in csv.DictReader(f,delimiter="\t"):
   key=(r["locus"],r["variant_GRCh37_REF_ALT"])
   if key in G and G[key]["BBJ_beta_ALT"]!=r["BBJ_beta_ALT"]:raise ValueError("Duplicate GWAS mismatch")
   G[key]=r
 if len(G)!=5666:raise ValueError("BBJ unique SNP count drift")
 summary=[];scenarios=[];cand=[];uncertainty=[];rng=np.random.default_rng(20261011)
 for locus,lead in LEADS.items():
  ids,kE,X=load(eas/(locus+".EAS504.fullcis.raw"))
  idsU,kU,U=load(eur/(locus+".EUR503.raw"))
  if len(ids)!=504 or len(idsU)!=503 or set(kE)!=set(kU):raise ValueError("Unmatched reference")
  if not set(ids).issubset(set(pop["sample"])):raise ValueError("Unknown EAS subject")
  mask=np.array([i in jpt for i in ids]);assert sum(mask)==104
  groups={"JPT104":(kE,X[mask]),"otherEAS400":(kE,X[~mask]),
          "EAS504":(kE,X),"EUR503":(kU,U)}
  gwas={k[1]:v for k,v in G.items() if k[0]==locus}
  zlead=float(gwas[lead]["BBJ_beta_ALT"])/float(gwas[lead]["BBJ_se"])
  ref={}
  for name,(keys,A) in groups.items():
   if lead not in keys:raise ValueError("Missing lead")
   idx={key:i for i,key in enumerate(keys)}
   r=corr(A,idx[lead]);ref[name]=(keys,A,r,idx)
   shared=[key for key in gwas if key in idx]
   counts={a:0 for a in ALPHAS};evals={a:0 for a in ALPHAS}
   for key in shared:
    if key==lead:continue
    rr=r[idx[key]]
    if not np.isfinite(rr):continue
    z=float(gwas[key]["BBJ_beta_ALT"])/float(gwas[key]["BBJ_se"])
    for alpha in ALPHAS:
     q=(1-alpha)*rr
     if 1-q*q<=1e-9:continue
     zz=(z-q*zlead)/np.sqrt(1-q*q)
     evals[alpha]+=1;counts[alpha]+=abs(zz)>=5.45
   for alpha in ALPHAS:
    scenarios.append(dict(locus=locus,population=name,alpha=alpha,
      evaluable=evals[alpha],heuristic_abs_z_ge_5p45=counts[alpha],
      inference_gate="REFERENCE_ONLY_NOT_CONDITIONAL_GWAS"))
   summary.append(dict(locus=locus,population=name,n_samples=len(A),
     n_SNP=len(keys),n_lead_informative=int(np.isfinite(r).sum()),
     n_monomorphic=int(np.isnan(r).sum()),lead_ALT_AF=float(A[:,idx[lead]].mean()/2)))
  if locus=="BBJ_IS_L002":
   z=float(gwas[CAND]["BBJ_beta_ALT"])/float(gwas[CAND]["BBJ_se"])
   bootstraps={}
   for name,(keys,A,rr,idx) in ref.items():
    r=float(rr[idx[CAND]]);x=A[:,idx[lead]];y=A[:,idx[CAND]]
    for alpha in ALPHAS:
     q=(1-alpha)*r;zz=(z-q*zlead)/np.sqrt(1-q*q)
     cand.append(dict(population=name,reference_n=len(A),alpha=alpha,
       signed_r=r,r2=r*r,lead_ALT_AF=float(x.mean()/2),
       candidate_ALT_AF=float(y.mean()/2),BBJ_marginal_z=z,
       BBJ_lead_z=zlead,approx_residual_z=zz,abs_z_ge_5p45=abs(zz)>=5.45))
    b=boot(x,y,z,zlead,rng);bootstraps[name]=b
    uncertainty.append(dict(population=name,bootstrap_reps=len(b),
       signed_r_observed=r,bootstrap_r_CI_low=float(np.quantile(b[:,0],.025)),
       bootstrap_r_CI_high=float(np.quantile(b[:,0],.975)),
       bootstrap_r_positive_rate=float(np.mean(b[:,0]>0)),
       heuristic_z_005_CI_low=float(np.quantile(b[:,1],.025)),
       heuristic_z_005_CI_high=float(np.quantile(b[:,1],.975)),
       bootstrap_heuristic_abs_z_ge_5p45_rate=float(np.mean(np.abs(b[:,1])>=5.45))))
   for name in ("otherEAS400","EUR503"):
    delta=bootstraps["JPT104"][:,0]-bootstraps[name][:,0]
    uncertainty.append(dict(population="JPT104_minus_"+name,bootstrap_reps=len(delta),
      signed_r_observed=float(ref["JPT104"][2][ref["JPT104"][3][CAND]]-ref[name][2][ref[name][3][CAND]]),
      bootstrap_r_CI_low=float(np.quantile(delta,.025)),
      bootstrap_r_CI_high=float(np.quantile(delta,.975)),
      bootstrap_r_positive_rate=float(np.mean(delta>0)),
      heuristic_z_005_CI_low="",heuristic_z_005_CI_high="",
      bootstrap_heuristic_abs_z_ge_5p45_rate=""))
 out.mkdir(parents=True,exist_ok=True)
 for file,rows in (("JPT_REFERENCE_LOCUS_QC.tsv",summary),
                  ("JPT_EAS_EUR_ALPHA_SCREEN.tsv",scenarios),
                  ("CHR10_JPT_EAS_EUR_CANDIDATE_ALPHA.tsv",cand),
                  ("CHR10_JPT_EAS_EUR_BOOTSTRAP.tsv",uncertainty)):save(out,file,rows)
 meta=dict(status="EXTERNAL_REFERENCE_ONLY_NOT_FORMAL_CONDITIONAL_GWAS",
   JPT_samples=104,EAS_non_JPT_samples=400,full_EAS_samples=504,EUR_samples=503,
   bootstrap_replicates=2000,seed=20261011,alpha_grid=ALPHAS,
   historical_coloc_646_unchanged=True,new_ABF_120_unchanged=True,
   original_discovery_genes=2225,study_donor_LD_unverified=True)
 (out/"JPT_REFERENCE_BOOTSTRAP_MANIFEST.json").write_text(json.dumps(meta,indent=2)+"\n")
 return meta
if __name__=="__main__":
 p=argparse.ArgumentParser()
 for k in ("source","eas","eur","panel","out"):p.add_argument("--"+k,type=Path,required=True)
 a=p.parse_args();print(json.dumps(run(a.source,a.eas,a.eur,a.panel,a.out),indent=2))
