#!/usr/bin/env python3
"""Descriptive allele-exact BBJ IS vs GIGASTROKE EAS AIS sign agreement.
Correlated SNPs are NOT independent replication observations.
"""
import csv,json,argparse,math
from pathlib import Path
GENES=("FGF5","C4orf22","CALHM2","NEURL")
def run(inp,out):
 with inp.open() as f:rows=list(csv.DictReader(f,delimiter="\t"))
 if len(rows)!=20689:raise ValueError("Frozen 4-gene source drift")
 summaries=[]
 for gene in GENES:
  subset=[r for r in rows if r["symbol"]==gene and r["BBJ_exact_Fasta37"]=="True" and r["GIGA_EAS_AIS_exact_Fasta37"]=="True" and r["palindromic_review"]=="False"]
  if len(subset)<2500:raise ValueError("Unexpected GWAS overlap")
  ok=0;q=0;bbj=0;giga=0;both=0;trip=0
  for r in subset:
   b=float(r["BBJ_beta_ALT"]);g=float(r["GIGA_EAS_AIS_beta_ALT"])
   if not math.isfinite(b) or not math.isfinite(g):raise ValueError("Nonfinite beta")
   ok+=(b*g)>0
   qq=float(r["QTD000161_nominal_p"])<1e-4
   bg=float(r["BBJ_p"])<5e-8;gg=float(r["GIGA_EAS_AIS_p"])<5e-8
   q+=qq;bbj+=bg;giga+=gg;both+=bg and gg;trip+=qq and bg and gg
  summaries.append(dict(gene=gene,locus=subset[0]["locus"],
    nonpal_exact_common_snp_gene_records=len(subset),
    same_BBJ_GIGA_ALT_beta_direction=ok,descriptive_direction_fraction=ok/len(subset),
    QTD_nominal_p_lt_1e4=q,BBJ_p_lt_5e8=bbj,
    GIGA_EAS_AIS_p_lt_5e8=giga,both_stroke_GWAS_p_lt_5e8=both,
    QTD_BBJ_GIGA_triple_marginal_p_threshold=trip,
    SNPs_not_independent_due_to_LD=True,
    cohorts_not_attested_independent=True,
    no_meta_analysis_or_heterogeneity_p_performed=True,
    no_causal_coloc_claim=True))
 out.mkdir(parents=True,exist_ok=True)
 with (out/"IS_PRIORITY4_BBJ_GIGA_ALLELE_EXACT_DIRECTION.tsv").open("w") as f:
  w=csv.DictWriter(f,fieldnames=list(summaries[0]),delimiter="\t");w.writeheader();w.writerows(summaries)
 (out/"IS_PRIORITY4_BBJ_GIGA_DIRECTION_MANIFEST.json").write_text(json.dumps({
  "audit":"DESCRIPTIVE_BBJ_VS_GIGASTROKE_EAS_AIS_EFFECT_DIRECTION",
  "source_variant_gene_rows":len(rows),"summary":summaries,
  "two_genomic_loci_not_four_independent_loci":True,
  "no_independent_replication_established":True,
  "original_2225_genes_untouched":True,
  "original_646_coloc_and_new_120_ABF_untouched":True},indent=2)+"\n")
 return summaries
if __name__=="__main__":
 p=argparse.ArgumentParser()
 p.add_argument("--input",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
 a=p.parse_args()
 for r in run(a.input,a.out):print(r)
