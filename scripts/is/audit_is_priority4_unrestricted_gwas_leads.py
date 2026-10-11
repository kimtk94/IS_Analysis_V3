#!/usr/bin/env python3
"""Unrestricted BBJ IS vs GIGASTROKE EAS AIS regional association leads.
GWAS regions NOT filtered by QTL, reference LD, or SNP allele class.
"""
import csv,gzip,json,argparse,math
from pathlib import Path
import pandas as pd,numpy as np
REG={"4":(80087872,82356697),"10":(104113502,106353104)}
def select(path):
 data={c:[] for c in REG};total=0
 with gzip.open(path,"rt") as f:
  for x in csv.DictReader(f,delimiter="\t"):
   total+=1;c=x["chr"]
   if c not in REG:continue
   pos=int(x["pos"])
   if not REG[c][0]<=pos<=REG[c][1]:continue
   try:pv=float(x["p"])
   except ValueError:continue
   if not math.isfinite(pv) or pv<0 or pv>1:continue
   data[c].append(x)
 if total<6000000:raise ValueError("GWAS source truncated")
 return total,data
def do(bbj,giga,reference,out):
 nb,B=select(bbj);ng,G=select(giga)
 summaries=[];cross=[]
 for ch in REG:
  b=min(B[ch],key=lambda r:float(r["p"]))
  g=min(G[ch],key=lambda r:float(r["p"]))
  # study A versus study B variant id; preserve IDs verbatim and flag unavailable cross lookup
  bb={x["variant_id"]:x for x in B[ch]};gg={x["variant_id"]:x for x in G[ch]}
  for tag,source in (("BBJ",b),("GIGA",g)):
   v=source["variant_id"]
   cross.append(dict(chr=ch,lead_source=tag,source_variant=v,
      BBJ_p=bb[v]["p"] if v in bb else "",
      GIGA_p=gg[v]["p"] if v in gg else "",
      BBJ_beta=bb[v]["beta"] if v in bb else "",
      GIGA_beta=gg[v]["beta"] if v in gg else "",
      variant_id_canonical_allele_match_unverified=True))
  loc="BBJ_IS_L001" if ch=="4" else "BBJ_IS_L002"
  df=pd.read_csv(reference/f"{loc}.EAS504.fullcis.raw",sep="\t")
  ids=[x.rsplit("_",1)[0] for x in df.columns[6:]]
  r=""
  if b["variant_id"] in ids and g["variant_id"] in ids:
   i=ids.index(b["variant_id"]);j=ids.index(g["variant_id"])
   xx=df.iloc[:,i+6].to_numpy(float);yy=df.iloc[:,j+6].to_numpy(float)
   if np.std(xx)>0 and np.std(yy)>0:r=float(np.corrcoef(xx,yy)[0,1])
  summaries.append(dict(chr=ch,interval=f"{REG[ch][0]}-{REG[ch][1]}",
   n_BBJ_regional=len(B[ch]),n_GIGA_regional=len(G[ch]),
   BBJ_minp_variant=b["variant_id"],BBJ_minp=b["p"],
   GIGA_minp_variant=g["variant_id"],GIGA_minp=g["p"],
   same_variant=b["variant_id"]==g["variant_id"],
   EAS_reference_r_if_direct_ID_found=r,
   EAS_reference_r2_if_direct_ID_found=r*r if r!="" else "",
   scope="ENTIRE_CANONICAL_REGION_IN_EACH_GWAS_NOT_QTL_RESTRICTED",
   LD_and_variant_alleles_recheck_required=True))
 out.mkdir(parents=True,exist_ok=True)
 for nm,rows in (("IS_UNRESTRICTED_BBJ_GIGA_LEADS.tsv",summaries),("IS_UNRESTRICTED_BBJ_GIGA_CROSS_LOOKUP.tsv",cross)):
  with (out/nm).open("w") as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
 result=dict(BBJ_full_source_rows=nb,GIGA_AIS_EAS_full_source_rows=ng,
  results=summaries,reference_LD_diagnostic_only=True,no_independent_signal_inference=True,
  QTL_not_used_for_lead_selection=True,original_2225_IS_gene_universe_preserved=True)
 (out/"IS_UNRESTRICTED_GWAS_LEADS_MANIFEST.json").write_text(json.dumps(result,indent=2)+"\n")
 return result
if __name__=="__main__":
 p=argparse.ArgumentParser()
 for v in ("bbj","giga","reference","out"):p.add_argument("--"+v,required=True,type=Path)
 a=p.parse_args()
 print(json.dumps(do(a.bbj,a.giga,a.reference,a.out),indent=2))
