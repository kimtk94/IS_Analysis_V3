#!/usr/bin/env python3
"""Merge 720 EAS5 gene-tissue cis-tested-completeness records before/after extension.

Do NOT silently convert 'n_variants equals count' to confirmed identical
test variant identities, EAS QTL ancestry, full GWAS overlap or positive coloc.
Five provisional discovery regions, 240 original positional genes, three public
eQTL Catalogue studies, no controlled-access studies.
"""
import argparse,csv,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
REGIONS={
 "IS_XDATA_G0015":71,"IS_XDATA_G0023":14,"IS_XDATA_G0007":10,
 "IS_XDATA_G0025":92,"IS_XDATA_G0018":53
}
FULL="MATCHES_PERMUTED_FULL_CIS_VARIANT_COUNT"
UNTESTED="NO_PERMUTED_GENE_TEST_ON_OFFICIAL_SOURCE"
def read(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def write(p,rows):
 with p.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
  w.writeheader();w.writerows(rows)
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for part in iter(lambda:f.read(4194304),b""):h.update(part)
 return h.hexdigest()
def audit(baseline,extended,case_fraction):
 if len(baseline)!=720 or len(extended)!=720:
  raise ValueError("Original 720 gene-tissue cases and extended audit required")
 before={(x["region_id"],x["QTL_dataset"],x["gene_id"]):x for x in baseline}
 after={(x["region_id"],x["QTL_dataset"],x["gene_id"]):x for x in extended}
 if len(before)!=720 or len(after)!=720 or set(before)!=set(after):
  raise ValueError("Incomplete or duplicated frozen gene-tissue grid")
 if not case_fraction or case_fraction["case_fraction_status"]!="SECONDARY_PUBLISHED_RECONCILED_NOT_ORIGINAL_COHORT_LEVEL":
  raise ValueError("AIS source sample status needs verified metadata provenance")
 if int(case_fraction["secondary_published_cases"])+int(case_fraction["secondary_published_controls"])!=256274:
  raise ValueError("AIS case-control source arithmetic changed")
 rows=[]
 byregion=defaultdict(Counter)
 for key,b in sorted(before.items()):
  a=after[key]
  old_n=b["ebi_permuted_full_cis_n_variants"]
  new_n=a["full_cis_official_ebi_permuted_n_variants"]
  if old_n!=new_n:raise ValueError("Official EBI permutation variant universe changed")
  if b["gene_symbol"]!=a["gene_symbol"] or (int(b["archived_ABF_assays"])>0)!=(a["historical_ABF_tested"]=="True"):
   raise ValueError("Original gene symbol or ABF lineage drift")
  observed_before=int(b["nominal_regional_distinct_SNP_INDEL_variants"])
  observed_after=int(a["distinct_extended_cis_variant_ids"])
  if observed_after<observed_before:
   raise ValueError("A wider cis genomic QTL window lost source SNPs")
  prior=b["ebi_source_test_status"];current=a["cis_completeness_status"]
  if prior=="FULL_CIS_EQUIVALENCE_VS_EBI_NVARIANTS_COUNT_ONLY" and current!=FULL:
   raise ValueError("Previously count-complete cis tested universe was lost")
  if prior=="NO_EBI_PERMUTED_GENE_RECORD_CANNOT_ASSESS" and current!=UNTESTED:
   raise ValueError("Original permutation-tested target universe changed")
  if current==UNTESTED and observed_after!=0:
   raise ValueError("Unreported source permutation test still has QTL association rows")
  if current==FULL and (int(new_n)<=0 or observed_after!=int(new_n)):
   raise ValueError("Cis complete source count does not equal official permutation variant count")
  if current==FULL and a.get("official_permutation_lead_in_merged_nominal_QTL")!="True":
   raise ValueError("Official full-cis lead variant is not independently present in the extended nominal source")
  if prior=="REGIONAL_SLICE_INCOMPLETE_CIS_NOT_SAFE_FOR_FULL_LOCUS_COLOC" and observed_before>=int(old_n):
   raise ValueError("Baseline partial-cis label inconsistent with count")
  rid=key[0];stat=byregion[rid]
  stat["all_gene_dataset_pairs"]+=1
  stat["old_count_complete"]+=int(prior=="FULL_CIS_EQUIVALENCE_VS_EBI_NVARIANTS_COUNT_ONLY")
  stat["new_count_complete"]+=int(current==FULL)
  stat["still_incomplete"]+=int(current=="EXTENDED_CIS_STILL_TRUNCATED")
  stat["perm_not_reported"]+=int(current==UNTESTED)
  stat["newly_recovered_complete"]+=int(prior!="FULL_CIS_EQUIVALENCE_VS_EBI_NVARIANTS_COUNT_ONLY" and current==FULL)
  rows.append({
   "region_id":rid,"QTL_dataset":key[1],"gene_id":key[2],
   "gene_symbol":b["gene_symbol"],
   "archived_BBJ_GTEx_ABF_gene_tissue_assays":b["archived_ABF_assays"],
   "original_2MB_slice_unique_cis_variants":observed_before,
   "extended_3to4MB_source_unique_cis_variants":observed_after,
   "official_permuted_n_variants":new_n,
   "ebi_permuted_p_beta_source":b.get("ebi_permuted_gene_p_beta",""),
   "permuted_official_lead_variant_in_extended_QTL":a.get("official_permutation_lead_in_merged_nominal_QTL",""),
   "official_permuted_gene_test_exists":current!=UNTESTED,
   "original_cis_completeness":prior,
   "extended_cis_completeness":current,
   "restored_full_cis_n_variants_count":current==FULL,
   "original_QTL_slice_is_full_multiset_subset_of_wide":True,
   "QTL_tested_cis_variant_identity_not_independently_proven":True,
   "extended_region_GWAS_nonpalindromic_variant_harmonization_done":False,
   "EAS_AIS_secondary_case_fraction":case_fraction["secondary_fraction_cases"],
   "EAS_AIS_original_cohort_case_ratio_attested":False,
   "molecular_QTL_ancestry_matched_GWAS_and_LD":False,
   "coloc_recomputed":False,
   "is_causal_gene_established":False
  })
 counts=[]
 for rid,n in REGIONS.items():
  c=byregion[rid]
  if c["all_gene_dataset_pairs"]!=n*3:raise ValueError("Region gene-tissue universe drift")
  counts.append({
   "region_id":rid,"distinct_positional_genes":n,
   "gene_QTL_dataset_pairs":c["all_gene_dataset_pairs"],
   "original_full_cis_count_matched":c["old_count_complete"],
   "new_full_cis_count_matched":c["new_count_complete"],
   "newly_restored_full_cis_count_matched":c["newly_recovered_complete"],
   "remaining_partial":c["still_incomplete"],
   "no_permutation_record":c["perm_not_reported"],
   "new_coloc_causing_variant_claim":"NONE"
  })
 if sum(x["original_full_cis_count_matched"] for x in counts)!=10:
  raise ValueError("Baseline 10 count-equal cis tests changed")
 if sum(x["no_permutation_record"] for x in counts)!=327:
  raise ValueError("327 source permutation nonreported statuses changed")
 return rows,counts
def main():
 p=argparse.ArgumentParser()
 p.add_argument("--baseline",required=True,type=Path)
 p.add_argument("--case-table",required=True,type=Path)
 p.add_argument("--extended-root",required=True,type=Path)
 p.add_argument("--out-dir",required=True,type=Path)
 a=p.parse_args()
 if a.out_dir.exists() and any(a.out_dir.iterdir()):
  raise FileExistsError("Refuse output overwrite")
 for z in (a.baseline,a.case_table):
  if not z.is_file():raise FileNotFoundError(z)
 studies=[]
 paths=[a.baseline,a.case_table]
 for rid in REGIONS:
  f=a.extended_root/rid/f"{rid}_EXTENDED_CIS_GENE_QTL_UNIVERSE.tsv"
  if not f.is_file():raise FileNotFoundError(f)
  studies+=read(f);paths.append(f)
 cases=[x for x in read(a.case_table) if x["trait"]=="AIS"]
 if len(cases)!=1:raise ValueError("EAS AIS metadata duplicate or missing")
 rows,reg=audit(read(a.baseline),studies,cases[0])
 a.out_dir.mkdir(parents=True,exist_ok=True)
 write(a.out_dir/"IS_EAS5_720_QTL_FULLCIS_BEFORE_AFTER_READINESS.tsv",rows)
 write(a.out_dir/"IS_EAS5_5REGION_QTL_FULLCIS_RECOVERY_SUMMARY.tsv",reg)
 m={
  "status":"ALL_5_PUBLIC_NOMINAL_QTL_REGIONS_PERMUTED_CIS_COUNT_RESTORATION",
  "original_positional_gene_universe":2225,
  "5_regions_candidate_gene_ids":sum(REGIONS.values()),
  "QTL_genes_times_tissues":len(rows),
  "cis_complete_before":sum(x["original_full_cis_count_matched"] for x in reg),
  "cis_complete_after":sum(x["new_full_cis_count_matched"] for x in reg),
  "newly_recovered_count_equal":sum(x["newly_restored_full_cis_count_matched"] for x in reg),
  "still_partial":sum(x["remaining_partial"] for x in reg),
  "not_reported_in_public_permutation_table":sum(x["no_permutation_record"] for x in reg),
  "is_complete_variant_identity_verified_independently":False,
  "is_molecular_QTL_ancestry_matched":False,
  "is_original_EAS_AIS_case_control_cohort_table_attested":False,
  "is_ext_GWAS_allele_harmonization_done":False,
  "can_claim_shared_causal_variant":False,
  "SHA256_sources":{str(z):sha(z) for z in paths}
 }
 (a.out_dir/"IS_EAS5_COMPLETE_CIS_SCOPE_STATUS.json").write_text(json.dumps(m,indent=2)+"\n")
 print("IS_EAS5_CIS_RECOVERY_INTEGRATION_PASS",json.dumps({
  k:m[k] for k in ("QTL_genes_times_tissues","cis_complete_before","cis_complete_after",
  "newly_recovered_count_equal","still_partial","not_reported_in_public_permutation_table")}),flush=True)
if __name__=="__main__":main()
