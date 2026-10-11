#!/usr/bin/env python3
"""Verify GCTA .ma and PLINK BED same variants and GWAS/EAS frequency compatibility.

Never claim COJO is publishable with 1000 Genomes EAS504 reference (< 4k).
"""
import argparse,csv,json,hashlib,math
from collections import Counter
from pathlib import Path
LOCI=("BBJ_IS_L001","BBJ_IS_L002")
TOL=.2
def source_ma(fp):
 with fp.open(newline="") as f:r=list(csv.DictReader(f,delimiter="\t"))
 if len(r)<2500:raise ValueError("Too few GWAS variants")
 if any(set(x)!={"SNP","A1","A2","freq","b","se","p","N"} for x in r):
  raise ValueError("Required GCTA SNP columns changed")
 return {x["SNP"]:x for x in r}
def run(pre,ref,out):
 if out.exists() and any(out.iterdir()):raise FileExistsError("Do not overwrite output")
 out.mkdir(parents=True,exist_ok=True)
 res=[]
 for locus in LOCI:
  s=source_ma(pre/f"{locus}.candidate.ma")
  p=ref/f"{locus}.EAS504.COJO"
  bim=Path(str(p)+".bim");bed=Path(str(p)+".bed");fam=Path(str(p)+".fam")
  afreq=ref/f"{locus}.EAS504.freq.afreq"
  if not all(x.is_file() for x in (bim,bed,fam,afreq)):
   raise ValueError("Missing PLINK alleles/genotype sample frequencies")
  n_samples=sum(1 for _ in fam.open())
  if n_samples!=504:raise ValueError("Reference sample count not 504")
  genos={}
  for ln in bim.read_text().splitlines():
   v=ln.split()
   if len(v)!=6:raise ValueError("PLINK BIM schema changed")
   chrom,vid,cm,pos,a1,a2=v
   parts=vid.split(":")
   if len(parts)!=4 or chrom!=parts[0] or pos!=parts[1]:
    raise ValueError("BIM variant genomic ID mismatch")
   if set((a1.upper(),a2.upper()))!=set((parts[2].upper(),parts[3].upper())):
    raise ValueError("BIM genotype allele mismatch")
   if vid in genos:raise ValueError("Duplicate typed variant")
   genos[vid]=v
  fr={}
  with afreq.open() as f:
   rd=csv.DictReader(f,delimiter="\t")
   needed={"ID","REF","ALT","ALT_FREQS","OBS_CT"}
   if not needed.issubset(rd.fieldnames or ()):raise ValueError("PLINK2 afreq columns changed")
   for row in rd:
    k=row["ID"]
    if k in fr:raise ValueError("Duplicate frequency row")
    parts=k.split(":")
    if len(parts)!=4 or row["REF"]!=parts[2] or row["ALT"]!=parts[3]:
     raise ValueError("PLINK reference frequency allele not canonical ALT")
    freq=float(row["ALT_FREQS"])
    if not math.isfinite(freq) or freq<0 or freq>1:raise ValueError("Invalid ALT reference AF")
    fr[k]=freq
  if set(fr)!=set(genos):raise ValueError("Reference BIM and frequency genotype sets differ")
  valid={};q=Counter();absdiff=[]
  for k,x in s.items():
   if k not in genos:q["NOT_TYPED_IN_EAS_BED"]+=1;continue
   v=fr[k];freq=float(x["freq"])
   if not 0.01<=v<=0.99:q["EAS_MAF_LT_0P01"]+=1;continue
   if abs(freq-v)>TOL:q["EAS_BBJ_FREQ_DIFFERENCE_GT_0P2"]+=1;continue
   if x["A1"]!=k.split(":")[3] or x["A2"]!=k.split(":")[2]:
    raise ValueError("GCTA A1 frequency not exact ALT")
   valid[k]=x;absdiff.append(abs(freq-v))
   q["GCTA_FILE_ALLELE_FREQ_READY"]+=1
  if len(valid)<2000:raise ValueError("Insufficient reference-matched SNPs")
  with (out/f"{locus}.EAS504.INFO0p8.MAF0p01.diffFreq0p2.ma").open("w",newline="") as f:
   w=csv.DictWriter(f,fieldnames=["SNP","A1","A2","freq","b","se","p","N"],delimiter="\t")
   w.writeheader()
   for k in sorted(valid,key=lambda k:int(k.split(":")[1])):w.writerow(valid[k])
  res.append(dict(locus=locus,reference_EAS_N=n_samples,
   input_GWAS_candidate_snps=len(s),genotyped_EAS_SNPs=len(genos),
   GCTA_format_QC_ready=len(valid),not_typed=q["NOT_TYPED_IN_EAS_BED"],
   excluded_ref_MAF=q["EAS_MAF_LT_0P01"],
   excluded_ancestry_AF_delta_gt_0p2=q["EAS_BBJ_FREQ_DIFFERENCE_GT_0P2"],
   median_abs_EAS_BBJ_AF_delta=sorted(absdiff)[len(absdiff)//2],
   max_abs_EAS_BBJ_AF_delta=max(absdiff),
   GCTA_COJO_publication_gate="FAIL_EAS504_IS_1000G_TOO_SMALL",
   formal_conditional_GWAS_executed=False))
 (out/"IS_PRIORITY4_COJO_EAS504_REFERENCE_MATCH_QC.json").write_text(json.dumps({
  "results":res,"scientific_readiness":"SOURCE_FILE_PREP_ONLY_NOT_ACCEPTABLE_FORMAL_COJO_REFERENCE",
  "GCTA_official_recommended_min_reference_gt_4000":True,
  "no_formal_COJO_run":True,"historical_coloc_unchanged":True,
  "IS_2225_positional_gene_universe_preserved":True},indent=2)+"\n")
 return res
if __name__=="__main__":
 p=argparse.ArgumentParser()
 for v in ("pre","reference","out"):p.add_argument("--"+v,type=Path,required=True)
 a=p.parse_args()
 for row in run(a.pre,a.reference,a.out):print("COJO_READY_DIAGNOSTIC",row,flush=True)
