#!/usr/bin/env python3
"""Read-only EUR503 regional VCF to original BBJ IS source GRCh37 allele overlap."""
import argparse,csv,json,hashlib,subprocess
from pathlib import Path
def read_tsv(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def vcf_records(p):
 cmd=["bcftools","query","-f","%CHROM\t%POS\t%REF\t%ALT\n",str(p)]
 out=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True)
 seen=set()
 for l in out.stdout.splitlines():
  c,pos,ref,alts=l.split("\t")
  for alt in alts.split(","):seen.add((c,pos,ref.upper(),alt.upper()))
 return seen
def run(replay,source,vcfdir,out):
 r=[x for x in read_tsv(replay) if x["status"]=="PASS"]
 if len(r)!=30:raise ValueError("Original replay PASS count changed")
 panels={}
 for c in ("4","10"):
  p=vcfdir/f"IS_BBJ_LD_chr{c}.EUR503.b37.vcf.gz"
  if not p.is_file() or not Path(str(p)+".tbi").is_file():raise ValueError("Missing indexed EUR VCF "+str(p))
  samples=subprocess.check_output(["bcftools","query","-l",str(p)],text=True).splitlines()
  if len(samples)!=503 or len(set(samples))!=503:raise ValueError("EUR donor panel count failed")
  panels[c]=(vcf_records(p),p)
 records=[]
 for x in r:
  if x["gene_base"] not in ("ENSG00000138675","ENSG00000197826","ENSG00000138172","ENSG00000107954"):continue
  source_rows=read_tsv(source/x["filename"])
  chrom=source_rows[0]["chromosome"]
  if chrom not in panels:raise ValueError("Unexpected target chromosome")
  target=set()
  for v in source_rows:
   parts=v["variant_id"].split(":")
   if len(parts)!=4 or parts[0]!=chrom:raise ValueError("Original GRCh37 SNP variant ID invalid")
   target.add((parts[0],parts[1],parts[2].upper(),parts[3].upper()))
  panel,vcf=panels[chrom]
  direct=len(target&panel)
  flipped=sum((c,pos,alt,ref) in panel for c,pos,ref,alt in target)
  records.append({"gene_id":x["gene_base"],"locus":x["locus"],"dataset_key":x["dataset_key"],
   "original_snp_rows":len(source_rows),"unique_b37_snp_ids":len(target),
   "eur503_exact_b37_chrpos_ref_alt_overlap":direct,
   "eur503_ref_alt_swapped_variants":flipped,
   "coverage_fraction":direct/len(target),
   "regional_eur_vcf":str(vcf),
   "GTEx_cohort_LD_available":False,"coloc_validated":False})
 if not records:raise ValueError("No priority-four assays")
 out.mkdir(parents=True,exist_ok=True)
 with (out/"IS_PRIORITY4_EUR503_REFERENCE_OVERLAP.tsv").open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter="\t");w.writeheader();w.writerows(records)
 j={"audited_priority_four_gene_assays":len(records),
  "unique_stable_genes":len({x["gene_id"] for x in records}),
  "eur503_chr4_region_variants":len(panels["4"][0]),"eur503_chr10_region_variants":len(panels["10"][0]),
  "assays_90pct_allele_matched":sum(x["coverage_fraction"]>=.9 for x in records),
  "assays_50pct_allele_matched":sum(x["coverage_fraction"]>=.5 for x in records),
  "max_coverage":max(x["coverage_fraction"] for x in records),
  "source_vcf_sha256":{chrom:hashlib.sha256(p.read_bytes()).hexdigest() for chrom,(v,p) in panels.items()},
  "interpretation":"1000G EUR external reference source overlap, not matched GTEx cohort LD or validated conditional coloc"}
 (out/"IS_PRIORITY4_EUR503_REFERENCE_OVERLAP.json").write_text(json.dumps(j,indent=2)+"\n")
 return j
if __name__=="__main__":
 ap=argparse.ArgumentParser()
 for n in ("replay","source","vcfdir","out"):ap.add_argument("--"+n,type=Path,required=True)
 a=ap.parse_args();print(json.dumps(run(a.replay,a.source,a.vcfdir,a.out),indent=2))
