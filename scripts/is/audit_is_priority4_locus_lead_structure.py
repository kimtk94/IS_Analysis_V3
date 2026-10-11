#!/usr/bin/env python3
"""BBJ versus GIGASTROKE regional lead-variant structure; reference-LD descriptive."""
import argparse,csv,json,math
from pathlib import Path
import pandas as pd,numpy as np
LOC={"BBJ_IS_L001":("4","4:81182554:T:C"),"BBJ_IS_L002":("10","10:105459834:T:C")}
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def run(source,eas,out):
 rr=read(source)
 if len(rr)!=20689:raise ValueError("Frozen QTD/GWAS source changed")
 D={}
 for r in rr:
  if r["BBJ_exact_Fasta37"]!="True" or r["GIGA_EAS_AIS_exact_Fasta37"]!="True" or r["palindromic_review"]!="False":continue
  k=(r["locus"],r["SNP_GRCh37"])
  if k in D:
   old=D[k]
   if any(old[x]!=r[x] for x in ("BBJ_p","GIGA_EAS_AIS_p","BBJ_beta_ALT","GIGA_EAS_AIS_beta_ALT")):
    raise ValueError("GWAS source differs for same SNP across two gene annotations")
  D[k]=r
 summary=[];points=[]
 for locus,(chrom,prior_bbj_lead) in LOC.items():
  g={s:v for (loc,s),v in D.items() if loc==locus}
  bbj=min(g,key=lambda x:float(g[x]["BBJ_p"]))
  giga=min(g,key=lambda x:float(g[x]["GIGA_EAS_AIS_p"]))
  frame=pd.read_csv(eas/f"{locus}.EAS504.fullcis.raw",sep="\t")
  a=frame.iloc[:,6:].to_numpy(dtype=float,copy=True);keys=[]
  if len(frame)!=504:raise ValueError("Incorrect EAS reference sample N")
  for j,col in enumerate(frame.columns[6:]):
   k,count=col.rsplit("_",1);parts=k.split(":")
   if count==parts[2]:a[:,j]=2-a[:,j]
   elif count!=parts[3]:raise ValueError("Genotype counted allele drift")
   keys.append(k)
  idx={k:i for i,k in enumerate(keys)}
  if len(idx)!=len(keys):raise ValueError("Duplicate variants")
  if bbj not in idx or giga not in idx or prior_bbj_lead not in idx:raise ValueError("Missing GWAS lead in EAS")
  x=a[:,idx[bbj]];y=a[:,idx[giga]]
  if not np.isfinite(x).all() or not np.isfinite(y).all():raise ValueError("Missing lead dosage")
  r=float(np.corrcoef(x,y)[0,1])
  if not math.isfinite(r):raise ValueError("Monomorphic lead")
  identical=bbj==giga
  concordance=int(sum(float(v["BBJ_beta_ALT"])*float(v["GIGA_EAS_AIS_beta_ALT"])>0 for v in g.values()))
  # effective independent SNP count is not number of tags: never label as replication.
  summary.append(dict(locus=locus,common_nonpal_allele_exact_snps=len(g),
     minP_BBJ_variant=bbj,minP_BBJ=float(g[bbj]["BBJ_p"]),
     minP_GIGA_variant=giga,minP_GIGA=float(g[giga]["GIGA_EAS_AIS_p"]),
     same_lead_variant=identical,EAS504_signed_r_between_leads=r,EAS504_r2_between_leads=r*r,
     lead_sign_BBJ_GIGA_agreement=concordance/len(g),
     both_gwas_genomewide_significant_same_snp=sum(float(v["BBJ_p"])<5e-8 and float(v["GIGA_EAS_AIS_p"])<5e-8 for v in g.values()),
     loci_not_independent_of_other_gene_labels=True,
     LD_proxy_not_in_study=True,independent_replication_not_proven=True))
  for snp in (bbj,giga):
   v=g[snp]
   points.append(dict(locus=locus,lead_for="BBJ" if snp==bbj else "GIGASTROKE",
    snp=snp,BBJ_p=v["BBJ_p"],GIGA_EAS_AIS_p=v["GIGA_EAS_AIS_p"],
    BBJ_beta_ALT=v["BBJ_beta_ALT"],GIGA_beta_ALT=v["GIGA_EAS_AIS_beta_ALT"]))
 out.mkdir(parents=True,exist_ok=True)
 for nm,items in (("IS_PRIORITY4_BBJ_GIGA_REGIONAL_LEAD_STRUCTURE.tsv",summary),
                  ("IS_PRIORITY4_BBJ_GIGA_LEAD_CROSS_LOOKUP.tsv",points)):
  with (out/nm).open("w") as f:
   w=csv.DictWriter(f,fieldnames=list(items[0]),delimiter="\t");w.writeheader();w.writerows(items)
 (out/"IS_PRIORITY4_REGIONAL_LEAD_STRUCTURE_MANIFEST.json").write_text(json.dumps({
  "data_scope":"FOUR_GENE_EBI_CIS_QTL_SHARED_TWO_GWAS_NONPAL_SNP_UNIVERSE",
  "not_entire_studywide_GWAS_region":True,"no_causal_replicability_inference":True,
  "original_2225_candidate_universe_unchanged":True,"summary":summary},indent=2)+"\n")
 return summary
if __name__=="__main__":
 p=argparse.ArgumentParser()
 for n in ("source","eas","out"):p.add_argument("--"+n,type=Path,required=True)
 a=p.parse_args()
 print(json.dumps(run(a.source,a.eas,a.out),indent=2))
