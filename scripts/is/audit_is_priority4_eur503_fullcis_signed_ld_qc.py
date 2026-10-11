#!/usr/bin/env python3
"""EUR503 1000G external reference ALT dosage QC / lead signed LD; not cohort LD."""
import argparse,csv,json,math,hashlib
from pathlib import Path
import numpy as np
import pandas as pd

LEADS={"BBJ_IS_L001":"4:81182554:T:C","BBJ_IS_L002":"10:105459834:T:C"}
EXPECTED={"BBJ_IS_L001":2834,"BBJ_IS_L002":2828}
def run(indir,out):
 out.mkdir(parents=True,exist_ok=True)
 summaries=[];detail=[]
 for locus,lead in LEADS.items():
  p=indir/f"{locus}.EUR503.raw"
  frame=pd.read_csv(p,sep="\t")
  if len(frame)!=503:raise ValueError("Unexpected EUR population size")
  a=frame.iloc[:,6:].to_numpy(dtype=float,copy=True)
  names=[]
  for i,col in enumerate(frame.columns[6:]):
   key,allele=col.rsplit("_",1)
   chrom,pos,ref,alt=key.split(":")
   if allele==ref:a[:,i]=2-a[:,i]
   elif allele!=alt:raise ValueError("PLINK allele not REF or ALT")
   names.append(key)
  if len(names)!=EXPECTED[locus] or len(set(names))!=len(names) or lead not in names:
   raise ValueError("SNP variant or lead mismatch")
  if np.any(np.isfinite(a)&((a<0)|(a>2))):raise ValueError("Invalid ALT dosage")
  means=np.nanmean(a,axis=0)
  miss=np.mean(~np.isfinite(a),axis=0)
  maf=np.minimum(means/2,1-means/2)
  monomorphic=np.nanvar(a,axis=0)==0
  leadvec=a[:,names.index(lead)]
  r=[]
  for i in range(len(names)):
   x=a[:,i];v=np.isfinite(x)&np.isfinite(leadvec)
   if v.sum()<100 or monomorphic[i] or np.nanvar(leadvec[v])==0:r.append(float("nan"));continue
   r.append(float(np.corrcoef(x[v],leadvec[v])[0,1]))
  for i,key in enumerate(names):
   detail.append(dict(locus=locus,snp_b37=key,EUR503_ALT_EAF=means[i]/2,EUR503_MAF=maf[i],
    EUR503_missingness=miss[i],EUR503_monomorphic=bool(monomorphic[i]),
    signed_r_to_BBJ_GWAS_lead=r[i],r2_to_BBJ_GWAS_lead=(r[i]**2 if math.isfinite(r[i]) else ""),
    reference_only_not_GTex_donor_LD=True))
  informative=np.isfinite(r)
  # Small, fixed SNP slice for PSD QC and reproducible numerical diagnostics.
  keep=np.flatnonzero(np.isfinite(a).all(axis=0)&(~monomorphic))
  if len(keep)<100:raise ValueError("Too few complete-polymorphic LD SNPs")
  pick=keep[np.linspace(0,len(keep)-1, min(120,len(keep)),dtype=int)]
  sub=np.corrcoef(a[:,pick],rowvar=False)
  eig=np.linalg.eigvalsh(sub)
  if np.min(eig)<-1e-8:raise ValueError("Signed LD correlation matrix not PSD")
  summaries.append(dict(locus=locus,n_samples=len(frame),n_source_SNPs=len(names),
    n_missing_any=int(np.sum(miss>0)),max_missingness=float(np.max(miss)),
    n_monomorphic=int(np.sum(monomorphic)),n_informative_lead_r=int(np.sum(informative)),
    lead_variant=lead,lead_ALT_EAF=float(means[names.index(lead)]/2),
    ld_subset_n_variants=len(pick),ld_subset_min_eigenvalue=float(np.min(eig)),
    ld_subset_max_eigenvalue=float(np.max(eig)),
    status="REFERENCE_EUR503_SIGNED_LD_QC_PASS_NOT_GTex_STUDY_LD"))
 def save(path,rows):
  with path.open("w",newline="") as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
 save(out/"IS_PRIORITY4_EUR503_FULLCIS_LD_PER_VARIANT.tsv",detail)
 save(out/"IS_PRIORITY4_EUR503_FULLCIS_LD_QC_SUMMARY.tsv",summaries)
 report=dict(loci=summaries,n_total_variant_locus_rows=len(detail),
  eQTL_donor_cohort_LD_verified=False,full_SuSiE_coloc_performed=False,
  prior_ABF_H4_modified=False,candidate_universe=2225,
  source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(indir.glob("*.raw"))})
 (out/"IS_PRIORITY4_EUR503_FULLCIS_LD_QC_MANIFEST.json").write_text(json.dumps(report,indent=2)+"\n")
 return report
if __name__=="__main__":
 p=argparse.ArgumentParser()
 p.add_argument("--input",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
 a=p.parse_args()
 print(json.dumps(run(a.input,a.out),indent=2))
