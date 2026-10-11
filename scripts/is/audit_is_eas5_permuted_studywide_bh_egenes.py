#!/usr/bin/env python3
"""Read-only genome-wide study-specific BH adjustment of eQTL Catalogue p_beta.

Uses ALL molecular traits in each original QTL study's permuted.tsv.gz, not
only stroke-region candidate genes. q<=0.05 means a published-source
REPROCESSED cis-eQTL evidence signal, NOT GWAS coloc/shared causal allele.
GTEx Aorta/Cortex donor samples are related; do not count as independence.
"""
import argparse,csv,gzip,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
DATASETS={"QTD000131":"GTEx Aorta","QTD000171":"GTEx Cortex",
          "QTD000434":"ROSMAP DLPFC"}
REGIONS=("IS_XDATA_G0015","IS_XDATA_G0023","IS_XDATA_G0007",
         "IS_XDATA_G0025","IS_XDATA_G0018")
def rd(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def write(p,rows):
 with p.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for buf in iter(lambda:f.read(4194304),b""):h.update(buf)
 return h.hexdigest()
def bh_all(ps):
 """Return Benjamini-Hochberg adjusted p-values in original entry order."""
 if not ps:return []
 vals=[float(x) for x in ps]
 if any(not math.isfinite(x) or not 0<=x<=1 for x in vals):
  raise ValueError("Invalid permutation p_beta passed to BH correction")
 order=sorted(range(len(vals)),key=lambda i:vals[i])
 q=[None]*len(vals);high=1.
 for rankidx in range(len(vals)-1,-1,-1):
  i=order[rankidx]
  high=min(high, vals[i]*len(vals)/(rankidx+1),1.)
  q[i]=high
 return q
def read_perm(path):
 rows=[]
 with gzip.open(path,"rt") as f:
  reader=csv.DictReader(f,delimiter="\t")
  need={"molecular_trait_object_id","n_variants","p_beta","variant"}
  if not need.issubset(reader.fieldnames or ()):raise ValueError("EBI permuted schema changed")
  for row in reader:rows.append(row)
 gene_ids=[x["molecular_trait_object_id"].split(".")[0] for x in rows]
 if len(gene_ids)!=len(set(gene_ids)):raise ValueError("Duplicated official gene permutation traits")
 good=[i for i,r in enumerate(rows) if r["p_beta"] not in ("","NA",".")]
 p=[float(rows[i]["p_beta"]) for i in good]
 q=bh_all(p)
 d={}
 for r in rows:
  gene=r["molecular_trait_object_id"].split(".")[0]
  d[gene]={"p_beta":None,"BH_q":None,"n_variants":int(r["n_variants"]),
           "lead_variant":r["variant"]}
 for i,adjusted in zip(good,q):
  gene=gene_ids[i]
  d[gene]["p_beta"]=float(rows[i]["p_beta"])
  d[gene]["BH_q"]=adjusted
 return d,len(rows),len(p)
def integrate(merged,original,permd):
 if len(merged)!=720 or len(original)!=2225:raise ValueError("Original 720/2225 read-only input required")
 keys={(x["region_id"],x["QTL_dataset"],x["gene_id"]):x for x in merged}
 if len(keys)!=720:raise ValueError("Duplicated gene-region-dataset")
 by_gene={x["gene_id_stable"]:x for x in original}
 selected=[x for x in original if x["region_id"] in REGIONS]
 if len(selected)!=240:raise ValueError("EAS5 original 240 distinct positional genes changed")
 pp={};n_sources={};sh={}
 for did in DATASETS:
  fp=permd/f"{did}.permuted.tsv.gz"
  if not fp.is_file():raise FileNotFoundError(fp)
  pp[did],n,valid=read_perm(fp)
  sh[str(fp)]=sha(fp)
  n_sources[did]={"study":DATASETS[did],"original_all_permutation_molecular_traits":n,
                   "total_genes_with_finite_p_beta":valid,
                   "genome_wide_FDR_method":"BH_all_permuted_source_gene_tests"}
 rows=[]
 for (rid,did,gene),r in sorted(keys.items()):
  ref=pp[did].get(gene)
  perm_present=ref is not None
  complete=r["restored_full_cis_n_variants_count"]=="True"
  if perm_present!=(r["official_permuted_n_variants"] not in ("",None)):
   raise ValueError("Source permutation gene test presence mismatch")
  if perm_present and int(r["official_permuted_n_variants"])!=ref["n_variants"]:
   raise ValueError("Published gene cis tested variant count inconsistent")
  if perm_present and not complete:
   raise ValueError("Merged source claims incomplete published cis window")
  if complete and r["permuted_official_lead_variant_in_extended_QTL"]!="True":
   raise ValueError("Source full cis lead genomic SNP missing")
  q=ref["BH_q"] if ref else None
  p=ref["p_beta"] if ref else None
  rows.append({
   "region_id":rid,"QTL_dataset":did,"gene_id":gene,
   "gene_symbol":r["gene_symbol"],
   "archived_646_gene_tissue_ABF_count":int(r["archived_BBJ_GTEx_ABF_gene_tissue_assays"]),
   "official_permuted_cis_gene_test_exists":perm_present,
   "full_cis_count_and_official_lead_genomic_SNP_recovered":complete,
   "official_n_cis_tested_genetic_variants":ref["n_variants"] if ref else "",
   "official_p_beta":p if p is not None else "",
   "BH_qvalue_across_ALL_STUDY_tested_genes":q if q is not None else "",
   "BH05_ebi_reprocessed_cis_eGene":q is not None and q<=.05,
   "EBI_p_beta_not_negative_for_untested_genes":True,
   "source_covariates_and_RNA_test_population_genetic_ancestry_independently_verified":False,
   "AIS_case_control_s_original_cohort_phenotype_table_independently_verified":False,
   "cohort_matched_molecular_QTL_LD_or_shared_causal_SNP_proven":False,
   "science_gate":"EBI_REPROCESSED_MOLECULAR_QTL_EGENE_ONLY_NOT_IS_CAUSAL_GENE"
  })
 gene_rows=[]
 for g in selected:
  my=[x for x in rows if x["gene_id"]==g["gene_id_stable"]]
  if len(my)!=3:raise ValueError("One of three QTL sources missing per gene")
  valid=[x for x in my if x["official_permuted_cis_gene_test_exists"]]
  significant=[x for x in my if x["BH05_ebi_reprocessed_cis_eGene"]]
  new=int(g["archived_646_coloc_assays"])==0
  gene_rows.append({
   "region_id":g["region_id"],"gene_id":g["gene_id_stable"],
   "gene_symbol":g["canonical_gene_symbol"],
   "number_archived_BBJ_GTEx_ABF_assays":int(g["archived_646_coloc_assays"]),
   "previously_ABF_untested":new,
   "number_QTL_studies_with_official_fullcis_test":len(valid),
   "number_QTL_studies_with_genomewide_BH05_eQTL":len(significant),
   "lowest_BH_q_across_three_QTL_studies":min((float(x["BH_qvalue_across_ALL_STUDY_tested_genes"]) for x in valid if x["BH_qvalue_across_ALL_STUDY_tested_genes"]!=""),default=""),
   "dataset_eGene_FDR05":";".join(x["QTL_dataset"] for x in significant),
   "QTL_source_status":(
       "EBI_REPROCESSED_EGENE_BH05_EXPLORATORY_GWAS_COLOC_PENDING" if significant else
       "EBI_CIS_TESTED_NO_STUDYWIDE_BH05_SIGNAL" if valid else
       "NO_EBI_PERMUTED_GENE_TEST_NOT_NEGATIVE"),
   "shared_GWAS_QTL_causal_variant":"NOT_ESTABLISHED"
  })
 if len(gene_rows)!=240:raise ValueError("EAS5 candidate count changed")
 manifest={
  "status":"EBI_PERMUTED_STUDYWIDE_BENJAMINI_HOCHBERG_QTL_ONLY",
  "original_2225_gene_ids_preserved":True,
  "studied_region_genes":len(gene_rows),
  "gene_tissue_pairs":len(rows),
  "study_source_permutation_counts":n_sources,
  "qtl_gene_tissue_official_permuted_tests":sum(x["official_permuted_cis_gene_test_exists"] for x in rows),
  "qtl_gene_tissue_BH_q_le_0p05":sum(x["BH05_ebi_reprocessed_cis_eGene"] for x in rows),
  "unique_genes_BH_q_le_0p05":sum(x["number_QTL_studies_with_genomewide_BH05_eQTL"]>0 for x in gene_rows),
  "unique_new_unarchived_ABF_genes_BH_q_le_0p05":sum(x["previously_ABF_untested"] and x["number_QTL_studies_with_genomewide_BH05_eQTL"]>0 for x in gene_rows),
  "these_are_not_AIS_coloc_genes":True,
  "all_QTL_genome_wide_cis_eGenes_not_EAS_matched_genotype_LD":True,
  "BH_applied_to_ALL_source_study_permuted_genes_not_240_IS_subset":True,
  "SHA256_permuted_inputs":sh
 }
 return rows,gene_rows,manifest
def main():
 ap=argparse.ArgumentParser()
 for key in ("merged","original-gene-readiness","permuted-dir","out-dir"):
  ap.add_argument("--"+key,type=Path,required=True)
 a=ap.parse_args()
 for f in (a.merged,a.original_gene_readiness):
  if not f.is_file():raise FileNotFoundError(f)
 if a.out_dir.exists() and any(a.out_dir.iterdir()):raise FileExistsError("Refuse output overwrite")
 rows,genes,manifest=integrate(rd(a.merged),rd(a.original_gene_readiness),a.permuted_dir)
 a.out_dir.mkdir(parents=True,exist_ok=True)
 write(a.out_dir/"IS_EAS5_720_OFFICIAL_STUDYWIDE_BH_EGENE_SOURCE_TESTS.tsv",rows)
 write(a.out_dir/"IS_EAS5_240_POSITIONAL_GENE_PUBLIC_CIS_EGENE_EVIDENCE.tsv",genes)
 (a.out_dir/"IS_EAS5_240_EGENES_EBI_STUDYWIDE_BH_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")
 print("IS_EAS5_PUBLIC_EGENE_BH_STUDYWIDE_PASS",json.dumps({
   k:manifest[k] for k in ("qtl_gene_tissue_official_permuted_tests",
   "qtl_gene_tissue_BH_q_le_0p05","unique_genes_BH_q_le_0p05",
   "unique_new_unarchived_ABF_genes_BH_q_le_0p05")}),flush=True)
if __name__=="__main__":main()
