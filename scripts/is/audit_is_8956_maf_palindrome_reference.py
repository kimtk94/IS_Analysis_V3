#!/usr/bin/env python3
"""Read-only audit: frozen 8956 IS SNPs, BBJ vs annotated 1KG-EAS frequency.

This compares allele-aligned frequencies on the same REF/ALT. PVAR EAS_AF is
an annotation, not the PLINK EAS genotype-frequency test or GTEx genotype AF.
Never auto-flip palindromes or alter GWAS / QTL samples or causal conclusions.
"""
import argparse,csv,hashlib,json,math,re
from collections import Counter,defaultdict
from pathlib import Path
LOCI=("BBJ_IS_L001","BBJ_IS_L002","BBJ_IS_L003","BBJ_IS_L004")
def read(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for x in iter(lambda:f.read(4*1024*1024),b""):h.update(x)
 return h.hexdigest()
def save(p,rows):
 with p.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
  w.writeheader();w.writerows(rows)
def pvar(path,keys):
 matches={}
 with path.open() as f:
  rr=csv.DictReader((line for line in f if not line.startswith("##")),delimiter="\t")
  if not {"#CHROM","POS","ID","REF","ALT"}.issubset(rr.fieldnames or []):
   raise ValueError("Unexpected pvar schema")
  for r in rr:
   name=r["ID"]
   if name not in keys:continue
   if name in matches:raise ValueError("Duplicate PVAR variant")
   chrom,pos,ref,alt=name.split(":")
   if (r["#CHROM"],r["POS"],r["REF"],r["ALT"])!=(chrom,pos,ref,alt):
    raise ValueError("PVAR ID/alleles inconsistent")
   hit=re.search(r"(?:^|;)EAS_AF=([0-9.]+)(?:;|$)",r.get("INFO",""))
   if hit:
    af=float(hit.group(1))
    if not 0<=af<=1:raise ValueError("Bad EAS_AF")
    matches[name]=af
   else:matches[name]=None
 return matches
def audit(rows,refdir,qcdir,stage3dir):
 if len(rows)!=8956 or set(x["locus"] for x in rows)!=set(LOCI):
  raise ValueError("Expected 8956 SNPs at four loci")
 out=[];loc=[];provenance={}
 for locus in LOCI:
  rs=[x for x in rows if x["locus"]==locus]
  variants={x["variant_id_GRCh37"] for x in rs}
  if len(variants)!=len(rs):raise ValueError("Duplicate locus SNP")
  orig=refdir/f"{locus}.1KG_EAS.GRCh37.pvar"
  qc=qcdir/f"{locus}.QC.pvar"
  stage3=stage3dir/f"{locus}.tsv"
  for p in (orig,qc,stage3):
   if not p.is_file():raise FileNotFoundError(p)
   provenance[str(p)]=sha(p)
  original=pvar(orig,variants); filtered=pvar(qc,variants)
  st=read(stage3); stage={r["variant_id"]:r for r in st}
  if len(st)!=len(stage):raise ValueError("Duplicate Stage3 ID")
  if variants-original.keys() or variants-filtered.keys() or variants-stage.keys():
   raise ValueError("Archived coloc SNP absent from upstream source/QC/Stage3")
  for r in rs:
   if r["status"]!="PASS_NATIVE_EFFECT_BETA_AF_AND_CHAIN" or r["native_effect_vs_GRCh37_alt"]!="ALT":
    raise ValueError("Original BBJ effect direction not verified")
   id=r["variant_id_GRCh37"];ref,alt=r["variant_id_GRCh38"].split(":")[2:]
   af=float(r["native_AF_Allele2"]);maf=min(af,1-af)
   if abs(maf-float(stage[id]["maf"]))>1e-9:raise ValueError("Stage3 MAF does not equal BBJ MAF")
   ea=original[id]
   if ea is None:raise ValueError("Original EAS population AF missing")
   pal={ref,alt} in ({"A","T"},{"C","G"})
   near=pal and (.4<=af<=.6 or .4<=ea<=.6)
   delta=abs(af-ea)
   out.append({
    "locus":locus,"variant_GRCh37":id,"variant_GRCh38":r["variant_id_GRCh38"],
    "BBJ_effect_ALT_AF":af,"BBJ_MAF":maf,
    "original_1KG_EAS_ALT_AF_INFO":ea,"original_1KG_EAS_MAF_INFO":min(ea,1-ea),
    "absolute_BBJ_vs_EAS_ALT_AF_difference":delta,
    "palindromic_AT_or_CG":pal,"palindromic_near_half_either_population":near,
    "BBJ_EAS_AF_abs_delta_ge_0p1":delta>=.1,
    "BBJ_EAS_AF_abs_delta_ge_0p2":delta>=.2,
    "in_original_1KG_EAS_PVAR":True,"in_post_MAF_1pct_PVAR":True,
    "in_Stage3_summary":True,
    "risk_label":"PALINDROMIC_NEAR_HALF_REVIEW" if near else (
       "BBJ_EAS_AF_DELTA_GE_0p2_REVIEW" if delta>=.2 else (
       "PALINDROMIC_OTHER" if pal else "NON_PALINDROMIC")),
    "science_gate":"ALLELE_FREQUENCY_DESCRIPTIVE_NOT_MOLECULAR_LD_NOT_CAUSAL"
   })
  x=[t for t in out if t["locus"]==locus]
  loc.append({
   "locus":locus,"Stage3_regional_SNPs":len(st),"archived_coloc_SNPs":len(rs),
   "ratio_archived_over_stage3":len(rs)/len(st),
   "BBJ_MAF_lt_0p01":sum(y["BBJ_MAF"]<.01 for y in x),
   "BBJ_MAF_lt_0p05":sum(y["BBJ_MAF"]<.05 for y in x),
   "EAS_INFO_MAF_lt_0p01":sum(y["original_1KG_EAS_MAF_INFO"]<.01 for y in x),
   "palindromic_SNPs":sum(y["palindromic_AT_or_CG"] for y in x),
   "palindromic_near_half_either_population":sum(y["palindromic_near_half_either_population"] for y in x),
   "BBJ_EAS_AF_delta_ge_0p1":sum(y["BBJ_EAS_AF_abs_delta_ge_0p1"] for y in x),
   "BBJ_EAS_AF_delta_ge_0p2":sum(y["BBJ_EAS_AF_abs_delta_ge_0p2"] for y in x),
   "site_qc_status":"RETAINED_SNP_UNIVERSE_ONLY"
  })
 return out,loc,provenance
def main():
 p=argparse.ArgumentParser()
 for opt in ["native-audit","ref-dir","qc-dir","stage3-dir","out-dir"]:
  p.add_argument("--"+opt,required=True,type=Path)
 a=p.parse_args()
 for path in [a.native_audit]:
  if not path.is_file():raise FileNotFoundError(path)
 if a.out_dir.exists() and any(a.out_dir.iterdir()):raise FileExistsError("No output overwrite")
 out,loc,provenance=audit(read(a.native_audit),a.ref_dir,a.qc_dir,a.stage3_dir)
 a.out_dir.mkdir(parents=True,exist_ok=True)
 save(a.out_dir/"IS_8956_BBJ_EAS_PALINDROME_MAF_AUDIT.tsv",out)
 save(a.out_dir/"IS_FOUR_BBJ_LOCUS_SNP_SELECTION_FUNNEL.tsv",loc)
 summary={
  "status":"8956_BBJ_SOURCE_EFFECT_ALLELE_AND_EAS_PVAR_POPULATION_FREQUENCY_DESCRIPTIVE",
  "BBJ_archived_SNPs":len(out),"Stage3_SNPs":sum(x["Stage3_regional_SNPs"] for x in loc),
  "palindromic_SNPs":sum(x["palindromic_SNPs"] for x in loc),
  "palindromic_near_half_either_population":sum(x["palindromic_near_half_either_population"] for x in loc),
  "BBJ_EAS_alt_AF_delta_ge_0p1":sum(x["BBJ_EAS_AF_delta_ge_0p1"] for x in loc),
  "BBJ_EAS_alt_AF_delta_ge_0p2":sum(x["BBJ_EAS_AF_delta_ge_0p2"] for x in loc),
  "original_GWAS_MAF_below_0p01_in_archived_8956":sum(x["BBJ_MAF_lt_0p01"] for x in loc),
  "retained_subset_not_full_unfiltered_GWAS":True,
  "population_ALT_AF_annotation_not_actual_molecular_QTL_LD":True,
  "no_recomputed_coloc_H4_or_causal_gene_promotion":True,
  "input_SHA256":{str(a.native_audit):sha(a.native_audit),**provenance}
 }
 (a.out_dir/"IS_8956_MAF_PALINDROME_AUDIT_MANIFEST.json").write_text(json.dumps(summary,indent=2)+"\n")
 print("IS_8956_MAF_PALINDROME_PASS",json.dumps({k:v for k,v in summary.items() if isinstance(v,int)}))
if __name__=="__main__":main()
