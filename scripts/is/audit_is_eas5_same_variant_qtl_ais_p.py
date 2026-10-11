#!/usr/bin/env python3
"""Same-SNP QTL and East Asian AIS association intersections, AF-stratified.

This is a DESCRIPTIVE prior to colocalization: two marginal P values on the
same SNP can reflect one LD block, pleiotropy or cohort/sample structure.
No shared causal variant, no conditional analysis, no causal gene promotion.
"""
import argparse,csv,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
def hash_file(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for b in iter(lambda:f.read(4194304),b""):h.update(b)
 return h.hexdigest()
def write(p,rows):
 if not rows:raise ValueError("Empty")
 with p.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
def scan(rows):
 by=defaultdict(Counter);regions=defaultdict(Counter);rg_uniq=defaultdict(set)
 seen=set()
 for r in rows:
  key=(r["region_id"],r["QTL_dataset"],r["gene_id"],r["mapped_variant_GRCh37"])
  if key in seen:raise ValueError("Unremoved repeated QTL gene-SNP effect mapping")
  seen.add(key)
  x=by[key[:3]]; y=regions[r["region_id"]]
  x["n_QTL_source_variant_any_GWAS"]+=1;y["gene_SNP_any_GWAS_rows"]+=1
  b=r["BBJ_exact_SNP"]=="True"
  a=r["EAS_AIS_exact_SNP"]=="True"
  pal=r["palindromic_SNP"]=="True"
  if not (a and b) or pal:continue
  x["n_SNP_both_GWAS_nonpal"]+=1
  qa=float(r["QTL_p"]);wp=float(r["EAS_AIS_p"]);af=float(r["abs_QTL_vs_EAS_AIS_alt_eaf_delta"])
  if not 0<=qa<=1 or not 0<=wp<=1 or not 0<=af<=1 or not all(map(math.isfinite,(qa,wp,af))):
   raise ValueError("Invalid P-value or ALT frequency delta")
  y["gene_SNP_both_nonpal"]+=1
  low=af<=.1;verylow=af<=.05
  if low:x["n_both_nonpal_AF_delta_le_0p1"]+=1
  if verylow:x["n_both_nonpal_AF_delta_le_0p05"]+=1
  if qa<=1e-4:x["n_QTL_p_le_1e4"]+=1
  if wp<5e-8:x["n_AIS_p_lt_5e8"]+=1
  if qa<=1e-4 and wp<5e-8:
   x["n_same_SNP_QTL_1e4_AIS_5e8"]+=1;y["gene_SNP_same_variant_joint"]+=1
   rg_uniq[r["region_id"]].add(r["mapped_variant_GRCh37"])
   if low:
    x["n_joint_AF_delta_le_0p1"]+=1
   if verylow:
    x["n_joint_AF_delta_le_0p05"]+=1
   if not x.get("lead_joint_variant") or wp<float(x["lead_joint_AIS_p"]):
    x["lead_joint_variant"]=r["mapped_variant_GRCh37"]
    x["lead_joint_AIS_p"]=wp
    x["lead_joint_QTL_p"]=qa
  if "min_QTL_p_of_matched" not in x or qa<x["min_QTL_p_of_matched"]:
   x["min_QTL_p_of_matched"]=qa
  if "min_EAS_AIS_p_of_matched" not in x or wp<x["min_EAS_AIS_p_of_matched"]:
   x["min_EAS_AIS_p_of_matched"]=wp
 if not seen:raise ValueError("No allele-harmonized data")
 return by,regions,rg_uniq
def analyze(alleles,reference,output):
 if not alleles.is_file() or not reference.is_file():raise FileNotFoundError("missing source")
 with alleles.open(newline="") as f: s,regions,uniq=scan(csv.DictReader(f,delimiter="\t"))
 with reference.open(newline="") as f:
  gene_meta={(r["region_id"],r["QTL_dataset"],r["gene_id"]):r for r in csv.DictReader(f,delimiter="\t")}
 if len(gene_meta)!=393:raise ValueError("Original 393 source-tested gene-tissue records not preserved")
 out=[]
 for key,meta in sorted(gene_meta.items()):
  counts=s[key]
  if counts["n_SNP_both_GWAS_nonpal"]!=int(meta["exact_NONPALINDROMIC_SNP_in_BOTH_BBJ_AIS"]):
   raise ValueError("Same-gene SNP alignment count drift")
  out.append({
   "region_id":key[0],"QTL_dataset":key[1],
   "gene_id":key[2],"gene_symbol":meta["gene_symbol"],
   "archived_646_coloc_assays":meta["original_archived_ABF_assays"],
   "both_GWAS_nonpal_SNPs":counts["n_SNP_both_GWAS_nonpal"],
   "both_GWAS_nonpal_deltaAF_le_0p1":counts["n_both_nonpal_AF_delta_le_0p1"],
   "both_GWAS_nonpal_deltaAF_le_0p05":counts["n_both_nonpal_AF_delta_le_0p05"],
   "QTL_p_le_1e4_on_both_GWAS_SNPs":counts["n_QTL_p_le_1e4"],
   "AIS_p_lt_5e8_on_both_GWAS_SNPs":counts["n_AIS_p_lt_5e8"],
   "same_SNP_QTL_le_1e4_and_AIS_lt_5e8":counts["n_same_SNP_QTL_1e4_AIS_5e8"],
   "same_SNP_joint_and_deltaAF_le_0p1":counts["n_joint_AF_delta_le_0p1"],
   "same_SNP_joint_and_deltaAF_le_0p05":counts["n_joint_AF_delta_le_0p05"],
   "lead_joint_variant_GRCh37":counts.get("lead_joint_variant",""),
   "lead_joint_AIS_p":counts.get("lead_joint_AIS_p",""),
   "lead_joint_QTL_p":counts.get("lead_joint_QTL_p",""),
   "QTL_effect_allele_officially_ALT":True,
   "QTL_study_ancestry_independently_verified":False,
   "genotype_LD_matched_across_QTL_and_GWAS":False,
   "colocalization_or_causal_gene_proven":False,
   "science_gate":"SAME_VARIANT_MARGINAL_SIGNAL_ONLY"
  })
 if len(out)!=393:raise ValueError("Gene-tissue quantity drift")
 rg=[]
 for rid in sorted({x["region_id"] for x in out}):
  o=[x for x in out if x["region_id"]==rid]
  rg.append({
   "region_id":rid,
   "n_gene_tissue_source_tests":len(o),
   "n_gene_tissue_with_same_variant_QTL_and_AIS_threshold":sum(r["same_SNP_QTL_le_1e4_and_AIS_lt_5e8"]>0 for r in o),
   "n_distinct_genes_with_same_variant_QTL_and_AIS_threshold":len({r["gene_id"] for r in o if r["same_SNP_QTL_le_1e4_and_AIS_lt_5e8"]>0}),
   "n_gene_tissue_with_joint_and_AF_delta_le_0p1":sum(r["same_SNP_joint_and_deltaAF_le_0p1"]>0 for r in o),
   "n_distinct_genes_joint_and_AF_delta_le_0p1":len({r["gene_id"] for r in o if r["same_SNP_joint_and_deltaAF_le_0p1"]>0}),
   "n_unique_GRCh37_variants_meeting_same_SNP_threshold_in_any_QTL":len(uniq[rid]),
   "analysis_gate":"SAME_SNP_ASSOCIATIONS_NOT_COLOCALIZATION"
  })
 output.mkdir(parents=True,exist_ok=True)
 write(output/"IS_EAS5_SAME_SNP_ASSOCIATION_GENE_TISSUE_QC.tsv",out)
 write(output/"IS_EAS5_SAME_SNP_ASSOCIATION_BY_REGION_QC.tsv",rg)
 m={
  "status":"OFFICIAL_ALT_ALIGNED_MARGINAL_SAME_SNP_QTL_GWAS_THRESHOLD_AUDIT",
  "n_gene_tissue_tests":len(out),
  "n_gene_tissue_with_same_SNP_joint":sum(x["same_SNP_QTL_le_1e4_and_AIS_lt_5e8"]>0 for x in out),
  "n_gene_tissue_joint_with_AF_delta_le_0p1":sum(x["same_SNP_joint_and_deltaAF_le_0p1"]>0 for x in out),
  "n_distinct_genes_joint":len({x["gene_id"] for x in out if x["same_SNP_QTL_le_1e4_and_AIS_lt_5e8"]>0}),
  "n_distinct_new_untested_ABF_genes_joint":len({x["gene_id"] for x in out if x["same_SNP_QTL_le_1e4_and_AIS_lt_5e8"]>0 and int(x["archived_646_coloc_assays"])==0}),
  "n_distinct_new_genes_joint_with_AF_delta_le_0p1":len({x["gene_id"] for x in out if x["same_SNP_joint_and_deltaAF_le_0p1"]>0 and int(x["archived_646_coloc_assays"])==0}),
  "5e8_and_1e4_thresholds_descriptive_not_formal_coloc":True,
  "same_SNP_and_same_sign_not_shared_causality":True,
  "source_SHA256":{str(alleles):hash_file(alleles),str(reference):hash_file(reference)}
 }
 (output/"IS_EAS5_SAME_SNP_MARGINAL_QTL_GWAS_AUDIT.json").write_text(json.dumps(m,indent=2)+"\n")
 print("IS_EAS5_EXACT_VARIANT_SAME_SNP_QTL_AIS_MARGINAL_P_AUDIT_COMPLETE",json.dumps({
  k:m[k] for k in ("n_gene_tissue_tests","n_gene_tissue_with_same_SNP_joint",
   "n_distinct_genes_joint","n_distinct_new_untested_ABF_genes_joint",
   "n_distinct_new_genes_joint_with_AF_delta_le_0p1")}),flush=True)
 return m
def main():
 p=argparse.ArgumentParser()
 for k in ("alleles","gene-readiness","out"):
  p.add_argument("--"+k,type=Path,required=True)
 a=p.parse_args()
 if a.out.exists() and any(a.out.iterdir()):raise FileExistsError("No overwrite of prior QC output")
 analyze(a.alleles,a.gene_readiness,a.out)
if __name__=="__main__":main()
