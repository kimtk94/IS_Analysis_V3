#!/usr/bin/env python3
"""Track 240 EAS5 positional genes with independent public cis-eQTL source/FDR,
older narrow-window GWAS SNP match and explicit causal science gates.

Priority is planning/technical readiness, NOT evidence ranking or causal
stroke gene list. Study QTL samples aren't independent GWAS replications.
"""
import argparse,csv,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
REGION_TIER={
 "IS_XDATA_G0015":"T1_BBJ_AND_EAS_AIS_REGION_ASSOCIATED_NOT_INDEPENDENT_REPLICATION",
 "IS_XDATA_G0023":"T1_BBJ_AND_EAS_AIS_REGION_ASSOCIATED_NOT_INDEPENDENT_REPLICATION",
 "IS_XDATA_G0007":"T1_BBJ_AND_EAS_AIS_REGION_ASSOCIATED_NOT_INDEPENDENT_REPLICATION",
 "IS_XDATA_G0025":"T4_SUGGESTIVE_REGION_NO_5E8_SIGNAL_IN_PRIMARY_OUTCOMES",
 "IS_XDATA_G0018":"T4_SUGGESTIVE_REGION_NO_5E8_SIGNAL_IN_PRIMARY_OUTCOMES"
}
def rd(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for b in iter(lambda:f.read(1048576),b""):h.update(b)
 return h.hexdigest()
def write(p,rs):
 with p.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter="\t")
  w.writeheader();w.writerows(rs)
def build(gene,allele,regional,case):
 if len(gene)!=240 or len(allele)!=393 or len(regional)!=30:
  raise ValueError("Expect complete 240 genes and original 393 narrow-window QTL+GWAS intersection")
 if len({x["gene_id"] for x in gene})!=240:raise ValueError("Duplicated positional gene ID")
 rby={x["region_id"]:x for x in regional}
 if set(REGION_TIER)-set(rby):raise ValueError("Region GWAS source tier unavailable")
 for reg,expected in REGION_TIER.items():
  if rby[reg]["primary_GWAS_EAS_BBJ_AIS_tier"]!=expected:
   raise ValueError("Original T1/T4 regional GWAS tier drift")
 linked=defaultdict(list)
 for r in allele:linked[(r["region_id"],r["gene_id"])].append(r)
 out=[]
 for g in gene:
  rid,gid=g["region_id"],g["gene_id"]
  if rid not in REGION_TIER:raise ValueError("Outside pilot source region")
  prior=int(g["number_archived_BBJ_GTEx_ABF_assays"])
  observed=linked.get((rid,gid),[])
  cis=int(g["number_QTL_studies_with_official_fullcis_test"])
  qbh=int(g["number_QTL_studies_with_genomewide_BH05_eQTL"])
  if cis>3 or qbh>cis:raise ValueError("Invalid official permutation study-test count")
  isT1=rid in ("IS_XDATA_G0015","IS_XDATA_G0023","IS_XDATA_G0007")
  if qbh>0:category="PUBLIC_REPROCESSED_EGENE_BH05_SOURCE_ONLY"
  elif cis>0:category="PUBLIC_CIS_TESTED_BH05_NOT_PASS"
  else:category="PUBLIC_PERMUTATION_GENE_TEST_NOT_REPORTED"
  pipeline=("T1_QTL_EGENE_GWAS_RECHECK_REQUIRE_LD_AND_EXTENDED_MATCH" if isT1 and qbh>0 else
            "T1_QTL_NOT_SIGNIFICANT_OR_PERMUTED_NOT_REPORTED" if isT1 else
            "T4_QTL_EGENE_MOLECULAR_ONLY_AIS_ASSOCIATION_SUGGESTIVE" if qbh>0 else
            "T4_QTL_NOT_SIGNIFICANT_OR_PERMUTED_NOT_REPORTED")
  qtl_narrow=sum(int(x["exact_NONPALINDROMIC_SNP_in_BOTH_BBJ_AIS"]) for x in observed)
  if qbh and not observed:raise ValueError("FDR source gene without any prior source overlap records")
  out.append({
   "region_id":rid,"gene_id_stable":gid,"gene_symbol":g["gene_symbol"],
   "original_region_GWAS_phenotype_support":rby[rid]["primary_GWAS_EAS_BBJ_AIS_tier"],
   "original_gene_archived_646_ABF_assays":prior,
   "source_QTL_permuted_cis_studies":cis,
   "source_QTL_eGene_studywide_BH_q_le_0p05":qbh,
   "source_QTL_eGene_datasets":g["dataset_eGene_FDR05"],
   "source_QTL_min_BH_q":g["lowest_BH_q_across_three_QTL_studies"],
   "number_of_QTL_datasets_with_OLD_window_BBJ_AIS_source_intersections":len(observed),
   "original_2MB_window_both_GWAS_exact_nonpalindromic_geneSNP_records":qtl_narrow,
   "extended_full_cis_SNP_variant_GWAS_matching_rerun":False,
   "new_H4_colocalization_not_run":True,
   "official_QTL_genotype_LD_not_stroke_ancestry_matched":True,
   "stroke_GWAS_independent_replication_not_established":True,
   "technical_readiness_category":category,
   "workflow_priority_NOT_biological_gene_rank":pipeline,
   "gene_causal_status":"NOT_ESTABLISHED"
  })
 if len(out)!=240:raise ValueError("Original 240 gene list changed")
 if sum(int(r["source_QTL_eGene_studywide_BH_q_le_0p05"])>0 for r in out)!=66:
  raise ValueError("Published 66 QTL source eGenes changed")
 if sum(int(r["source_QTL_eGene_studywide_BH_q_le_0p05"])>0 and int(r["original_gene_archived_646_ABF_assays"])==0 for r in out)!=56:
  raise ValueError("56 unarchived QTL source eGene counts changed")
 regrows=[]
 for rid in REGION_TIER:
  items=[r for r in out if r["region_id"]==rid]
  regrows.append({
   "region_id":rid,"original_region_GWAS_support":REGION_TIER[rid],
   "n_positional_genes":len(items),
   "n_source_eGenes_BH05":sum(int(r["source_QTL_eGene_studywide_BH_q_le_0p05"])>0 for r in items),
   "n_source_eGenes_BH05_without_historical_ABF":sum(
     int(r["source_QTL_eGene_studywide_BH_q_le_0p05"])>0 and int(r["original_gene_archived_646_ABF_assays"])==0 for r in items),
   "n_source_cis_tested_genes":sum(int(r["source_QTL_permuted_cis_studies"])>0 for r in items),
   "n_fully_GWAS_harmonized_NEW_extended_cis_gene_sets":0,
   "n_new_positive_coloc_or_causal_genes":0,
   "status":"QTL_EGENE_SOURCE_ONLY_NOT_IS_CAUSALITY"
  })
 return out,regrows
def main():
 p=argparse.ArgumentParser()
 for k in ("egene240","old-gwas-overlap","region-queue","case-fraction","out-dir"):
  p.add_argument("--"+k,type=Path,required=True)
 a=p.parse_args()
 paths=(a.egene240,a.old_gwas_overlap,a.region_queue,a.case_fraction)
 for fp in paths:
  if not fp.is_file():raise FileNotFoundError(fp)
 if a.out_dir.exists() and any(a.out_dir.iterdir()):raise FileExistsError("No overwrite")
 case=[r for r in rd(a.case_fraction) if r["trait"]=="AIS"]
 if len(case)!=1 or case[0]["secondary_fraction_cases"]=="":
  raise ValueError("Case fraction provenance unavailable")
 genes,regions=build(rd(a.egene240),rd(a.old_gwas_overlap),rd(a.region_queue),case[0])
 a.out_dir.mkdir(parents=True,exist_ok=True)
 write(a.out_dir/"IS_EAS5_240_GENE_QTL_EGENE_GWAS_READINESS.tsv",genes)
 write(a.out_dir/"IS_EAS5_5_REGION_EGENE_STATUS.tsv",regions)
 summary={
  "status":"PUBLIC_STUDYWIDE_FDR_CIS_EGENE_AND_REGIONAL_GWAS_TECHNICAL_QUEUE_NOT_CAUSAL",
  "EAS5_original_positional_gene_ids":240,
  "QTL_source_eGenes_studywide_BH05":66,
  "new_QTL_source_eGenes_not_archived_ABF":56,
  "T1_candidate_eGenes_source_only":sum(int(r["n_source_eGenes_BH05"]) for r in regions if r["original_region_GWAS_support"].startswith("T1")),
  "T1_new_source_eGenes_not_archived_ABF":sum(int(r["n_source_eGenes_BH05_without_historical_ABF"]) for r in regions if r["original_region_GWAS_support"].startswith("T1")),
  "new_causal_genes":0,
  "new_coloc_results":0,
  "source_SHA256":{str(x):sha(x) for x in paths}
 }
 (a.out_dir/"IS_EAS5_240_GENE_EGENE_GWAS_QUEUE_MANIFEST.json").write_text(json.dumps(summary,indent=2)+"\n")
 print("IS_EAS5_240_STUDYWIDE_BH_EGENE_SOURCE_ONLY_TRIAGE_PASS",json.dumps({
  k:summary[k] for k in ("T1_candidate_eGenes_source_only","T1_new_source_eGenes_not_archived_ABF",
  "new_QTL_source_eGenes_not_archived_ABF","new_causal_genes")}),flush=True)
if __name__=="__main__":main()
