#!/usr/bin/env python3
"""Reciprocal BBJ/GIGA leads via chr:pos unordered alleles; audit FASTA orientation."""
import csv,gzip,json,math,subprocess,argparse
from pathlib import Path
PAIRS={"4":(81182554,81181072),"10":(105459834,105599770)}
def gather(path):
 out={}
 with gzip.open(path,"rt") as f:
  for x in csv.DictReader(f,delimiter="\t"):
   ch=x["chr"]
   if ch not in PAIRS or int(x["pos"]) not in PAIRS[ch]:continue
   k=(ch,int(x["pos"]))
   out.setdefault(k,[]).append(x)
 return out
def reference(fasta,ch,pos):
 p=subprocess.run(["samtools","faidx",str(fasta),f"{ch}:{pos}-{pos}"],text=True,capture_output=True,check=True)
 return p.stdout.splitlines()[-1].upper()
def parse_alleles(x,is_bbj):
 if is_bbj:
  return (x["ref"].upper(),x["alt"].upper())
 parts=x["variant_id"].split(":")
 if len(parts)!=4:raise ValueError("invalid GIGA ID")
 return (parts[2].upper(),parts[3].upper())
def run(bbj,giga,fasta,out):
 B=gather(bbj);G=gather(giga);r=[]
 for ch,positions in PAIRS.items():
  for pos in positions:
   base=reference(fasta,ch,pos);key=(ch,pos)
   ba=B.get(key,[]);ga=G.get(key,[])
   if len(ba)!=1 or len(ga)!=1:raise ValueError(f"Missing/duplicated lead site {key}: {len(ba)} {len(ga)}")
   b=ba[0];g=ga[0];br,bt=parse_alleles(b,True);gr,gt=parse_alleles(g,False)
   if set((br,bt))!=set((gr,gt)):raise ValueError("Locus allele pair mismatch")
   if base not in (br,bt) or base not in (gr,gt):raise ValueError("FASTA reference base incompatible")
   e1=b["effect_allele"].upper();e2=g["effect_allele"].upper()
   if e1 not in (br,bt) or e2 not in (gr,gt):raise ValueError("Effect allele unsupported")
   beta_b=float(b["beta"])*(1 if e1==base else -1)
   beta_g=float(g["beta"])*(1 if e2==base else -1)
   freq_b=float(b["eaf"]) if e1==base else 1-float(b["eaf"])
   freq_g=float(g["eaf"]) if e2==base else 1-float(g["eaf"])
   pal=set((br,bt)) in ({"A","T"},{"C","G"})
   r.append(dict(chrom=ch,pos=pos,REF_fasta37=base,BBJ_ID=b["variant_id"],
    GIGA_ID=g["variant_id"],BBJ_p=b["p"],GIGA_p=g["p"],
    BBJ_beta_fastaREF=beta_b,GIGA_beta_fastaREF=beta_g,
    BBJ_eaf_fastaREF=freq_b,GIGA_eaf_fastaREF=freq_g,
    same_FastaREF_effect_direction=beta_b*beta_g>0,
    palindrome=pal,strand_provenance_verified=not pal,
    frequency_abs_difference=abs(freq_b-freq_g),
    caution="PALINDROMIC_STRAND_UNRESOLVED" if pal else "FASTA_REF_BASE_ALIGNED_NOT_INDEPENDENT_REPLICATION"))
 out.mkdir(parents=True,exist_ok=True)
 with (out/"IS_PRIORITY4_UNRESTRICTED_RECIPROCAL_LEAD_SITE_ALLELES.tsv").open("w") as f:
  w=csv.DictWriter(f,fieldnames=list(r[0]),delimiter="\t");w.writeheader();w.writerows(r)
 (out/"IS_PRIORITY4_UNRESTRICTED_LEAD_ALLELE_MANIFEST.json").write_text(json.dumps({
  "sites":r,"no_independent_replication_inference":True,"palindrome_requires_strand_resolution":True,
  "full_source_GWAS_site_scan":True,"original_2225_gene_universe_unchanged":True},indent=2)+"\n")
 return r
if __name__=="__main__":
 p=argparse.ArgumentParser()
 for v in ("bbj","giga","fasta","out"):p.add_argument("--"+v,type=Path,required=True)
 a=p.parse_args()
 for v in run(a.bbj,a.giga,a.fasta,a.out):print(v,flush=True)
