#!/usr/bin/env python3
"""Verify genotype-count vs association semantics of archived and reprocessed GTEx.

This is a read-only allele-count and nominal QTL provenance diagnostic.
It MUST NOT convert beta to another allele or authorize mixed-source coloc.
"""
import argparse,csv,hashlib,json,math,statistics
from pathlib import Path
from collections import Counter
G={"ENSG00000138675":("FGF5","BBJ_IS_L001"),
   "ENSG00000197826":("C4orf22","BBJ_IS_L001"),
   "ENSG00000138172":("CALHM2","BBJ_IS_L002"),
   "ENSG00000107954":("NEURL","BBJ_IS_L002")}
TISSUE="GTEx_V8__Brain_Cerebellar_Hemisphere"
def read(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def safe_float(v):
 x=float(v)
 if not math.isfinite(x):raise ValueError("Nonfinite frequency, association or sample measure")
 return x
def classify(a,b,an):
 if a==b:return "AC_SAME"
 if a+b==an:return "AC_COMPLEMENT"
 return "AC_OTHER"
def run(source,qdir,out):
 q={}
 for locus in ("BBJ_IS_L001","BBJ_IS_L002"):
  rows=read(qdir/f"{locus}.QTD000161.region.priority4.tsv")
  for r in rows:
   if r["gene_id"] not in G:continue
   if G[r["gene_id"]][1]!=locus:raise ValueError("Unexpected molecular locus")
   k=(r["gene_id"],r["chromosome"],r["position"],r["ref"].upper(),r["alt"].upper())
   if k in q:
    # EBI may duplicate a variant when there are multiple dbSNP mappings.
    old=q[k]
    fields=("ac","an","maf","beta","se","pvalue")
    if any(old[x]!=r[x] for x in fields):raise ValueError("rsid duplicate has inconsistent QTL statistics")
   else:q[k]=r
 records=[];summaries=[]
 for gid,(symbol,locus) in G.items():
  p=source/f"{locus}__{TISSUE}__{gid}.tsv"
  oldrows=read(p)
  if len(oldrows)!=(1694 if locus=="BBJ_IS_L001" else 1583):raise ValueError("Source count changed")
  stats=Counter();fracdiff=[];p_logdiff=[];beta_same=[]
  for a in oldrows:
   parts=a["variant_id_hg38"].split(":")
   k=(gid,*parts)
   if k not in q:continue
   b=q[k]
   n=safe_float(a["an"]);m=safe_float(b["an"])
   ca=safe_float(a["ac"]);cb=safe_float(b["ac"])
   fa=safe_float(a["maf"]);fb=safe_float(b["maf"])
   if not (0<=ca<=n and 0<=cb<=m and 0<=fa<=.5+1e-7 and 0<=fb<=.5+1e-7):
    raise ValueError("Invalid allele count or MAF")
   # Allow six decimals of summary rounding, not numerical identity of inputs.
   if abs(fa-min(ca/n,1-ca/n))>2e-6:raise ValueError("Archived AC and MAF inconsistent")
   if abs(fb-min(cb/m,1-cb/m))>2e-6:raise ValueError("Catalogue AC and MAF inconsistent")
   match="AN_DIFFERENT" if n!=m else classify(ca,cb,n)
   stats[match]+=1
   if ca<=n/2:stats["ARCHIVE_AC_LE_HALF"]+=1
   if cb<=m/2:stats["CATALOGUE_AC_LE_HALF"]+=1
   if abs(fa-fb)<=2e-6:stats["MAF_MATCH_2E6"]+=1
   p1=safe_float(a["eqtl_p"]);p2=safe_float(b["pvalue"])
   be1=safe_float(a["eqtl_beta"]);be2=safe_float(b["beta"])
   if p1<0 or p1>1 or p2<0 or p2>1:raise ValueError("Invalid association p")
   p_logdiff.append(abs(math.log10(max(p1,1e-300))-math.log10(max(p2,1e-300))))
   fracdiff.append(abs(fa-fb))
   if be1*be2!=0:beta_same.append(be1*be2>0)
   records.append(dict(gene_id=gid,symbol=symbol,locus=locus,
     variant_hg38=":".join(parts),source_an=n,catalogue_an=m,
     source_ac=ca,catalogue_ac=cb,ac_relation=match,
     source_maf=fa,catalogue_maf=fb,MAF_abs_diff=abs(fa-fb),
     source_nominal_p=p1,catalogue_nominal_p=p2,
     source_beta=be1,catalogue_beta=be2,
     beta_sign_matches=(be1*be2>0 if be1*be2!=0 else ""),
     effect_allele_provenance="SOURCE_BETA_ALLELE_DIRECTION_NOT_INDEPENDENTLY_ATTESTED",
     inference_gate="NO_CAUSAL_REINTERPRETATION"))
  nrec=sum(stats[k] for k in ("AN_DIFFERENT","AC_SAME","AC_COMPLEMENT","AC_OTHER"))
  if nrec<=0:raise ValueError("No exact-allele matched associations")
  summaries.append(dict(symbol=symbol,gene_id=gid,locus=locus,exact_variant_rows=nrec,
     an_identical= nrec-stats["AN_DIFFERENT"],
     ac_identical=stats["AC_SAME"],ac_complement=stats["AC_COMPLEMENT"],
     ac_other=stats["AC_OTHER"],an_different=stats["AN_DIFFERENT"],
     archival_ac_at_most_an_half=stats["ARCHIVE_AC_LE_HALF"],
     catalogue_ac_at_most_an_half=stats["CATALOGUE_AC_LE_HALF"],
     maf_equal_2e6=stats["MAF_MATCH_2E6"],
     median_maf_abs_difference=statistics.median(fracdiff),
     median_abs_log10_nominal_p_difference=statistics.median(p_logdiff),
     beta_sign_concordance=sum(beta_same)/len(beta_same),
     beta_effect_direction_automatically_flipped=False,
     original_gtEx_association_model_verified=False,
     scientific_gate="GENOTYPE_MAF_EQUALITY_DOES_NOT_IMPLY_QTL_ESTIMATE_IDENTITY"))
 out.mkdir(parents=True,exist_ok=True)
 for fn,rr in [("IS_PRIORITY4_SNP_AC_MAF_BETA_LINEAGE.tsv",records),
               ("IS_PRIORITY4_AC_MAF_BETA_LINEAGE_SUMMARY.tsv",summaries)]:
  with (out/fn).open("w",newline="") as f:
   w=csv.DictWriter(f,fieldnames=list(rr[0]),delimiter="\t")
   w.writeheader();w.writerows(rr)
 manifest=dict(genes=4,snp_gene_allele_exact_matches=len(records),
   full_cis_qtl_reanalysis_is_not_original_GTEx=True,
   GTEx_donor_in_study_LD_verified=False,
   can_recalculate_archived_coloc_pp_h4=False,
   original_candidate_universe=2225,
   no_automatic_allele_beta_inversion=True,
   summaries=summaries,
   source_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(qdir.glob("BBJ_IS_L00*.QTD000161.region.priority4.tsv"))})
 (out/"IS_PRIORITY4_ALLELE_COUNT_LINEAGE_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")
 return manifest
if __name__=="__main__":
 a=argparse.ArgumentParser()
 for name in ("source","qdir","out"):a.add_argument("--"+name,required=True,type=Path)
 p=a.parse_args();r=run(p.source,p.qdir,p.out)
 for x in r["summaries"]:print("GENE",x)
