#!/usr/bin/env python3
"""Summarize freeze-safe COL4A1/COL4A2 paired ABF posterior sensitivity.

Runs entirely on locally generated derivatives; cannot promote biological
causality. Separate molecular cis-eQTL FDR, prior sensitivity, SNP coverage,
and unverified ancestry-matched LD. Existing 646 ABF and 2,225 genes untouched.
"""
import argparse,csv,hashlib,json,math,statistics
from collections import defaultdict
from pathlib import Path
REQUIRED={"COL4A1":"ENSG00000187498","COL4A2":"ENSG00000134871"}
DATASETS={"QTD000131":"GTEx Artery Aorta","QTD000171":"GTEx Brain Cortex",
          "QTD000434":"ROSMAP DLPFC"}
ORIGINAL_VARIANTS={
("COL4A1","QTD000131"):9025,("COL4A2","QTD000131"):9026,
("COL4A1","QTD000171"):8218,("COL4A2","QTD000171"):8218,
("COL4A1","QTD000434"):8763,("COL4A2","QTD000434"):8764
}
def rd(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for part in iter(lambda:f.read(1024*1024),b""):h.update(part)
 return h.hexdigest()
def wr(p,rows):
 with p.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
  w.writeheader();w.writerows(rows)
def pearson(a,b):
 if len(a)!=len(b) or len(a)<2:raise ValueError("QTL per SNP correlation needs paired variants")
 return statistics.correlation(a,b)
def summarize(sens,inp,permuted,legacy):
 if len(sens)!=270 or len(inp)!=26648:raise ValueError("Source assay or SNP data volume changed")
 old=defaultdict(dict)
 for r in inp:
  did=r["QTL_dataset"];g=r["gene_symbol"];vid=r["GRCh37_fasta_REF_ALT_variant"]
  if did not in DATASETS or g not in REQUIRED:raise ValueError("Unexpected tissue/gene")
  if vid in old[(did,g)]:raise ValueError("Duplicate locus-specific gene SNP")
  old[(did,g)][vid]=r
 primary={}
 for r in sens:
  k=(r["QTL_dataset"],r["gene_symbol"],r["SNP_selection"],
     round(float(r["QTL_N_multiplier"]),2),float(r["p12"]))
  if k in primary:raise ValueError("Duplicated p12 grid sample")
  primary[k]=r
 if len(primary)!=270:raise ValueError("Incomplete 270 posterior conditions")
 for r in sens:
  posterior=[float(r[f"PP_H{i}"]) for i in range(5)]
  if not all(math.isfinite(z) and 0<=z<=1 for z in posterior) or abs(sum(posterior)-1)>1e-5:
   raise ValueError("ABF posterior mass not normalized")
  if r["co_localized_gene_claim"]!="NOT_ESTABLISHED":
   raise ValueError("Posterior cannot promote causal gene")
  if r["EAS_BBJ_independent_replication_attested"].upper()!="FALSE" or r["QTL_EAS_ancestry_matched_LD_attested"].upper()!="FALSE":
   raise ValueError("New independent verification silently added")
 per={(r["region_id"],r["gene_id"],r["QTL_dataset"]):r for r in permuted}
 if len(per)!=720:raise ValueError("Original 720 permutation universe drift")
 legacy_col4=[r for r in legacy if r["gene_base"] in REQUIRED.values()]
 if len(legacy_col4)!=160:raise ValueError("Old 16 tissue x 2 gene x 5 priors expectation changed")
 rows=[];prior_rows=[];af_rows=[]
 for did,study in DATASETS.items():
  lhs=old[(did,"COL4A1")];rhs=old[(did,"COL4A2")]
  shared=set(lhs)&set(rhs)
  if shared!=set(lhs) or shared!=set(rhs):raise ValueError("Paired inputs for two genes not identical")
  v=sorted(shared)
  z1=[float(lhs[t]["QTL_beta_alt"])/float(lhs[t]["QTL_se"]) for t in v]
  z2=[float(rhs[t]["QTL_beta_alt"])/float(rhs[t]["QTL_se"]) for t in v]
  corr=pearson(z1,z2)
  for g in REQUIRED:
   sub=old[(did,g)]
   base=primary[(did,g,"ALL_NONPAL_PAIRED",1.0,1e-5)]
   gid=REQUIRED[g]
   sources=per[("IS_XDATA_G0023",gid,did)]
   study_q=sources["BH_qvalue_across_ALL_STUDY_tested_genes"]
   qsig=sources["BH05_ebi_reprocessed_cis_eGene"]=="True"
   pperm=sources["official_p_beta"]
   if int(base["nsnps"])!=len(v) or int(sources["official_n_cis_tested_genetic_variants"])!=ORIGINAL_VARIANTS[(g,did)]:
    raise ValueError("Frozen original QTL source and posterior SNP counts do not agree")
   if g=="COL4A2" and did=="QTD000434" and not qsig:
    raise ValueError("Known full-study BH COL4A2 finding lost")
   if g=="COL4A1" and did=="QTD000434" and qsig:
    raise ValueError("COL4A1 not passing FDR improperly promoted")
   row={
    "QTL_dataset":did,"study":study,"gene_symbol":g,"gene_id_stable":gid,
    "official_cis_tested_variant_count":ORIGINAL_VARIANTS[(g,did)],
    "matched_nonpalindromic_genomic_SNPs":len(v),
    "fraction_of_official_cis_variants_in_primary_ABF":len(v)/ORIGINAL_VARIANTS[(g,did)],
    "paired_QTL_genetic_Z_correlation_with_neighbor":corr,
    "EBI_studywide_BH_eGene":qsig,
    "EBI_official_permutation_p_beta":pperm,
    "EBI_studywide_BH_q":study_q,
    "primary_coloc_p1":1e-4,"primary_coloc_p2":1e-4,"primary_coloc_p12":1e-5,
    "primary_coloc_H3":base["PP_H3"],
    "primary_coloc_H4":base["PP_H4"],
    "primary_coloc_H4_over_H3":base["H4_over_H3"],
    "main_conditional_H4_SNP":base["top_H4_conditional_SNP"],
    "molecular_QTL_source_N_from_AN_over2":base["QTL_N_AN_over2_reference"],
    "primary_QTL_expression_sdY_estimated_from_varbeta_MAF_N":base["estimated_QTL_sdY"],
    "BBJ_GIGA_independent_case_cohorts":False,
    "QTL_GWAS_ancestry_matched_LD_or_multi_signal_SuSiE":False,
    "gene_causal_status":"NOT_ESTABLISHED_EXPLORATORY_SINGLE_CAUSAL_ABF_ONLY"
   }
   rows.append(row)
   for p12 in (1e-6,3e-6,1e-5,3e-5,1e-4):
    q=primary[(did,g,"ALL_NONPAL_PAIRED",1.0,p12)]
    prior_rows.append({"QTL_dataset":did,"gene_symbol":g,"p12":p12,
      "PP_H3":q["PP_H3"],"PP_H4":q["PP_H4"],
      "H4_over_H3":q["H4_over_H3"],
      "top_conditional_H4_SNP":q["top_H4_conditional_SNP"],
      "science_status":"PRIOR_SENSITIVE_SINGLE_SIGNAL_NO_CAUSAL_ASSIGNMENT"})
   for band in ("ALL_NONPAL_PAIRED","AF_DIFF_MAX_0P2_BOTH_GENES","AF_DIFF_MAX_0P1_BOTH_GENES"):
    for multiplier in (0.8,1.0,1.2):
     q=primary[(did,g,band,multiplier,1e-5)]
     af_rows.append({"QTL_dataset":did,"gene_symbol":g,
       "QTL_AIS_AF_selection":band,"QTL_N_multiplier":multiplier,
       "SNPs":q["nsnps"],"PP_H3":q["PP_H3"],"PP_H4":q["PP_H4"],
       "status":"AF_FILTER_SENSITIVITY_AND_QTL_N_PERTURBATION_ONLY"})
 if len(rows)!=6 or len(prior_rows)!=30 or len(af_rows)!=54:raise ValueError("Summary grid count error")
 coloc_by_genes={r["gene_symbol"]:r for r in rows if r["QTL_dataset"]=="QTD000434"}
 bothp=coloc_by_genes
 if_fail=False
 if (float(bothp["COL4A2"]["primary_coloc_H4"])<.9 or
   float(bothp["COL4A1"]["primary_coloc_H4"])<.9):
  raise ValueError("Published baseline posterior changed")
 old_p12={}
 for g,gid in REQUIRED.items():
  z=[x for x in legacy_col4 if x["gene_base"]==gid and abs(float(x["conditional_p12"])-1e-5)<1e-15]
  if len(z)!=16:raise ValueError("Old GTEx BBJ 16 tissue ABF grid inconsistent")
  old_p12[g]=max(float(x["PP_H4"]) for x in z)
 manifest={
  "status":"COL4A1_COL4A2_PAIRED_ABF_PRIOR_SENSITIVE_EXPLORATORY_NOT_NEW_STROKE_CAUSAL",
  "source_270_ABF_run_count":len(sens),
  "three_tissue_matched_6_assays":len(rows),
  "new_shared_causal_genes":0,
  "QTL_SNP_COVERAGE_INCOMPLETE_FOR_COL4A2_ROSMAP":len(old[("QTD000434","COL4A2")])<8764,
  "COL4A1_and_COL4A2_share_same_ROSMAP_top_conditioned_H4_SNP":bool(
   coloc_by_genes["COL4A1"]["main_conditional_H4_SNP"]==coloc_by_genes["COL4A2"]["main_conditional_H4_SNP"]),
  "COL4A1_ROSMAP_fdr_eGene":False,"COL4A2_ROSMAP_fdr_eGene":True,
  "COL4A1_ROSMAP_ppH4_p12_1e5":float(coloc_by_genes["COL4A1"]["primary_coloc_H4"]),
  "COL4A2_ROSMAP_ppH4_p12_1e5":float(coloc_by_genes["COL4A2"]["primary_coloc_H4"]),
  "old_646_BBJ_GTEx_16_tissue_max_H4_p12_1e5_NOT_COMPARABLE_DIRECTLY_WITH_ROSMAP":old_p12,
  "high_H4_in_two_genes_same_GWAS_SNP_set_is_not_locus_gene_resolution":True,
  "confirmed_sample_ancestry_matched_QTL_LD":False,
  "confirmed_multi_signal_coloc_susie":False,
  "study_N_derived_from_allele_number_AN2_not_original_donor_metadata":True,
  "GIGASTROKE_EAS_AIS_case_fraction_source_secondary_only":True,
  "meta_conclusion":"EXPLORATORY_SHARED_REGIONAL_REGULATORY_LOCUS_NOT_UNIQUE_GENE"
 }
 return rows,prior_rows,af_rows,manifest
def main():
 p=argparse.ArgumentParser()
 for k in ("abf","paired-input","source-fdr","legacy","out-dir"):
  p.add_argument("--"+k,type=Path,required=True)
 a=p.parse_args()
 for fp in (a.abf,a.paired_input,a.source_fdr,a.legacy):
  if not fp.is_file():raise FileNotFoundError(fp)
 if a.out_dir.exists() and any(a.out_dir.iterdir()):raise FileExistsError("Refuse overwrite")
 rs,prior,af,m=summarize(rd(a.abf),rd(a.paired_input),rd(a.source_fdr),rd(a.legacy))
 a.out_dir.mkdir(parents=True,exist_ok=True)
 wr(a.out_dir/"IS_COL4A1_COL4A2_6_TISSUE_PAIRED_ABF_EXPLORATORY_CONCLUSIONS.tsv",rs)
 wr(a.out_dir/"IS_COL4A1_COL4A2_30_PRIOR_H3_H4_SENSITIVITY.tsv",prior)
 wr(a.out_dir/"IS_COL4A1_COL4A2_54_AF_N_ABF_SENSITIVITY.tsv",af)
 m["input_SHA256"]={str(fp):sha(fp) for fp in (a.abf,a.paired_input,a.source_fdr,a.legacy)}
 (a.out_dir/"IS_COL4A1_COL4A2_PAIRED_ABF_SCIENCE_GATE_MANIFEST.json").write_text(json.dumps(m,indent=2)+"\n")
 print("IS_COL4A1_COL4A2_PAIRED_ABF_270_SCENARIOS_SCIENCE_GATE_PASS",json.dumps({
  k:m[k] for k in ("COL4A1_ROSMAP_ppH4_p12_1e5","COL4A2_ROSMAP_ppH4_p12_1e5",
  "COL4A1_and_COL4A2_share_same_ROSMAP_top_conditioned_H4_SNP","new_shared_causal_genes")}),flush=True)
if __name__=="__main__":main()
