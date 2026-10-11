#!/usr/bin/env python3
"""Prepare exact-allele BBJ IS GWAS COJO inputs against EAS504 genotypes.

Critical fail-closed gate: GCTA explicitly discourages 1000G LD (N=504);
no formal in-sample COJO inference is licensed by this script.
"""
import argparse,csv,gzip,hashlib,json,math,subprocess
from collections import Counter
from pathlib import Path

LOCI={"BBJ_IS_L001":("4",80087872,82356697,"4:81182554:T:C"),
      "BBJ_IS_L002":("10",104113502,106353104,"10:105459834:T:C")}
MIN_INFO=.8
MIN_FREQ=.01
def fmt(x):
 return format(float(x),".16g")
def hashfile(path):
 h=hashlib.sha256()
 with path.open("rb") as f:
  for part in iter(lambda:f.read(4*1024*1024),b""):h.update(part)
 return h.hexdigest()
def load_ref(path):
 raw=subprocess.run(["bcftools","query","-f","%CHROM\t%POS\t%REF\t%ALT\n",str(path)],
  capture_output=True,text=True,check=True)
 keys=set();multi=0;ambiguous=0
 for ln in raw.stdout.splitlines():
  q=ln.split("\t")
  if len(q)!=4:raise ValueError("Malformed VCF record")
  ch,pos,ref,alts=q
  if "," in alts:multi+=1;continue
  if len(ref)!=1 or len(alts)!=1 or ref not in "ACGT" or alts not in "ACGT":continue
  if {ref,alts} in ({"A","T"},{"C","G"}):ambiguous+=1;continue
  k=f"{ch}:{pos}:{ref}:{alts}"
  if k in keys:raise ValueError("Duplicated exact VCF key")
  keys.add(k)
 return keys,dict(multiallelic_records=multi,palindromic_excluded=ambiguous,total_valid_nonpal=len(keys))
def run(bbj,refdir,out):
 if out.exists() and any(out.iterdir()):raise FileExistsError("Non-overwrite scientific audit")
 references={}
 for locus,(ch,lo,hi,lead) in LOCI.items():
  fp=refdir/f"chr{ch}.EAS504.fullcis.b37.vcf.gz"
  keys,summary=load_ref(fp)
  if lead not in keys:raise ValueError("GWAS lead unavailable in EAS reference")
  references[locus]=(keys,summary,fp)
 staged={loc:{} for loc in LOCI};qc={loc:Counter() for loc in LOCI}
 with gzip.open(bbj,"rt") as f:
  r=csv.DictReader(f,delimiter="\t")
  req={"chr","pos","ref","alt","variant_id","effect_allele","other_allele","beta","se","p","eaf","info","n"}
  if not req.issubset(r.fieldnames or ()):raise ValueError("Changed BBJ GWAS source")
  scanned=0
  for x in r:
   scanned+=1
   chrom=x["chr"]
   if chrom not in ("4","10"):continue
   pos=int(x["pos"])
   locus="BBJ_IS_L001" if chrom=="4" else "BBJ_IS_L002"
   _,lo,hi,_=LOCI[locus]
   if not lo<=pos<=hi:continue
   q=qc[locus];q["BBJ_region_rows"]+=1
   ref=x["ref"].upper();alt=x["alt"].upper();ea=x["effect_allele"].upper();oa=x["other_allele"].upper()
   if len(ref)!=1 or len(alt)!=1 or ref not in "ACGT" or alt not in "ACGT":
    q["not_SNP"]+=1;continue
   k=f"{chrom}:{pos}:{ref}:{alt}"
   if x["variant_id"]!=k:raise ValueError("BBJ canonical genomic key mismatched")
   if k not in references[locus][0]:
    q["no_exact_nonpal_EAS_VCF_SNP"]+=1;continue
   if (ea,oa) not in ((ref,alt),(alt,ref)):raise ValueError("BBJ alleles not canonical")
   try:
    beta=float(x["beta"]);se=float(x["se"]);pv=float(x["p"]);eaf=float(x["eaf"])
    info=float(x["info"]);n=float(x["n"])
   except ValueError:
    q["missing_measure"]+=1;continue
   if not all(math.isfinite(z) for z in (beta,se,pv,eaf,info,n)) or se<=0 or not 0<pv<=1 or not 0<eaf<1 or n<=0:
    q["invalid_summary"]+=1;continue
   if info<MIN_INFO:q["low_INFO"]+=1;continue
   if min(eaf,1-eaf)<MIN_FREQ:q["low_GWAS_MAF"]+=1;continue
   beta_alt=beta if ea==alt else -beta
   freq_alt=eaf if ea==alt else 1-eaf
   ma=dict(SNP=k,A1=alt,A2=ref,freq=fmt(freq_alt),b=fmt(beta_alt),
           se=fmt(se),p=fmt(pv),N=fmt(n))
   if k in staged[locus]:
    if staged[locus][k]!=ma:raise ValueError("Conflicting duplicated SNP")
   else:staged[locus][k]=ma
 if scanned<13_000_000:raise ValueError("BBJ GWAS not scanned in full")
 for locus,(_,lo,hi,lead) in LOCI.items():
  if lead not in staged[locus]:raise ValueError("Required BBJ lead excluded")
  if len(staged[locus])<2500:raise ValueError("Unexpectedly low BBJ-EAS overlap")
 out.mkdir(parents=True)
 results=[]
 for locus,(chrom,lo,hi,lead) in LOCI.items():
  sample=staged[locus]
  order=sorted(sample,key=lambda k:int(k.split(":")[1]))
  ma=out/f"{locus}.candidate.ma"
  with ma.open("w") as f:
   w=csv.DictWriter(f,fieldnames=["SNP","A1","A2","freq","b","se","p","N"],delimiter="\t")
   w.writeheader();w.writerows(sample[k] for k in order)
  (out/f"{locus}.extract.snplist").write_text("\n".join(order)+"\n")
  (out/f"{locus}.lead.snplist").write_text(lead+"\n")
  nvals=Counter(float(x["N"]) for x in sample.values())
  q=qc[locus]
  results.append(dict(locus=locus,source_region=f"{chrom}:{lo}-{hi}",lead=lead,
    BBJ_regional_rows=q["BBJ_region_rows"],BBJ_exact_nonpal_INFO0p8_MAF0p01_N_valid=len(sample),
    region_not_SNP=q["not_SNP"],not_in_EAS_ref=q["no_exact_nonpal_EAS_VCF_SNP"],
    low_INFO=q["low_INFO"],low_BBJ_MAF=q["low_GWAS_MAF"],failed_numeric=q["missing_measure"]+q["invalid_summary"],
    gwas_N_unique_values=len(nvals),most_common_gwas_N=nvals.most_common(1)[0][0],
    reference_1000G_EAS_samples=504,GCTA_min_recommended_reference=4000,
    formal_COJO_reference_gate="FAIL_1000G_EAS504_NOT_RECOMMENDED",
    source_VCF=str(references[locus][2]),source_QC=references[locus][1]))
 manifest=dict(scanned_BBJ_rows=scanned,regions=results,COJO_candidate_inputs_prepared=True,
    full_COJO_scientific_gate="BLOCKED_UNTIL_LARGE_ANCESTRY_MATCHED_LD_PANEL",
    official_GCTA_reference_guidance="https://yanglab.westlake.edu.cn/software/gcta/",
    no_formal_COJO_executed=True,no_susie_coloc_executed=True,
    BBJ_IS_N_not_guessed=True,summary_INFO_threshold=MIN_INFO,
    BBJ_MAF_threshold=MIN_FREQ,no_palindromic_or_multiallelic_variant_ambiguity=True,
    historical_646_coloc_preserved=True,independent_ABF_120_preserved=True,
    original_discovery_universe_2225_preserved=True,sha256_source_BBJ=hashfile(bbj))
 (out/"IS_PRIORITY4_BBJ_EAS_COJO_PRELIMINARY_INPUT_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")
 return manifest
if __name__=="__main__":
 p=argparse.ArgumentParser()
 p.add_argument("--bbj",type=Path,required=True)
 p.add_argument("--reference",type=Path,required=True)
 p.add_argument("--out",type=Path,required=True)
 a=p.parse_args()
 d=run(a.bbj,a.reference,a.out)
 for r in d["regions"]:print("LOCUS",r["locus"],r["BBJ_exact_nonpal_INFO0p8_MAF0p01_N_valid"],r["most_common_gwas_N"],flush=True)
 print("NO_FORMAL_COJO_RAN",flush=True)
