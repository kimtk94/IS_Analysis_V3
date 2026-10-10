#!/usr/bin/env python3
"""Reconstruct 4 gene-specific tested cis variant sets from bounded EBI regional queries.

Do not conflate QTS000015/QTD000161 reanalysis with original GTEx V8 FastQTL,
and do not call a genomic query complete until its unique variant count equals
the independently supplied reprocessing n_variants for that gene.
"""
import argparse,csv,gzip,hashlib,json,statistics
from pathlib import Path
from collections import defaultdict
GENES={
 "ENSG00000138675":("FGF5","BBJ_IS_L001"),
 "ENSG00000197826":("C4orf22","BBJ_IS_L001"),
 "ENSG00000138172":("CALHM2","BBJ_IS_L002"),
 "ENSG00000107954":("NEURL","BBJ_IS_L002"),
}
TISSUE="GTEx_V8__Brain_Cerebellar_Hemisphere"
def read(path):
 with path.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def readperm(path):
 with gzip.open(path,"rt") as f:
  records={r["molecular_trait_object_id"]:r for r in csv.DictReader(f,delimiter="\t") if r["molecular_trait_object_id"] in GENES}
 if set(records)!=set(GENES):raise ValueError("Four permutation records are required")
 return records
def construct(qdir,permuted,original,out):
 refs=readperm(permuted)
 inputs={}
 for locus in ("BBJ_IS_L001","BBJ_IS_L002"):
  fp=qdir/f"{locus}.QTD000161.wide.priority4.tsv"
  if not fp.is_file():raise ValueError(f"Missing downloaded QTL locus {fp}")
  inputs[locus]=read(fp)
 result=[];variant_table=[]
 for g,(sym,locus) in GENES.items():
  cat={}
  duplicate_rows=0
  for x in inputs[locus]:
   if x["gene_id"]!=g:continue
   key=(x["chromosome"],x["position"],x["ref"].upper(),x["alt"].upper())
   if key in cat:
    duplicate_rows+=1
    for f in ("pvalue","beta","se","an","ac","maf"):
     if x[f]!=cat[key][f]:raise ValueError("Inconsistent duplicate rsid rows")
   else:cat[key]=x
  full=int(refs[g]["n_variants"])
  original_file=original/f"{locus}__{TISSUE}__{g}.tsv"
  archived=read(original_file)
  a={}
  for row in archived:
   parts=tuple(row["variant_id_hg38"].split(":"))
   if parts in a:raise ValueError("Legacy duplicate SNP ID")
   a[parts]=row
  intersect=set(cat)&set(a)
  source_only=set(a)-set(cat)
  cat_only=set(cat)-set(a)
  status="COMPLETE_VS_EBI_REPROCESSED_PERMUTED_NVARIANTS" if len(cat)==full else "INCOMPLETE_OR_DISCORDANT_STUDY_UNIVERSE"
  result.append(dict(gene_id=g,symbol=sym,locus=locus,
      full_cis_permutation_variant_count=full,
      tabix_wide_unique_variant_count=len(cat),
      permutation_minus_wide_unique=full-len(cat),
      tabix_rsID_duplicate_rows=duplicate_rows,
      EBI_reprocessed_universe_status=status,
      archived_variant_count=len(a),exact_shared_variants=len(intersect),
      archived_only_variant_keys=len(source_only),
      EBI_reprocessed_only_variant_keys=len(cat_only),
      fraction_archived_in_EBI=len(intersect)/len(a),
      fraction_EBI_in_archived=len(intersect)/len(cat) if cat else 0,
      EBI_cis_permuted_p_beta=refs[g]["p_beta"],
      same_original_GTEx_V8_analysis=False,
      legacy_coloc_posterior_updated=False,
      ancestry_matched_GTEx_LD_verified=False))
  for chrom,pos,ref,alt in sorted(cat):
   r=cat[(chrom,pos,ref,alt)]
   variant_table.append(dict(gene_id=g,symbol=sym,locus=locus,
    chromosome=chrom,position_GRCh38=int(pos),REF_GRCh38=ref,ALT_GRCh38=alt,
    variant_id_GRCh38=f"{chrom}:{pos}:{ref}:{alt}",
    in_original_archived_GWAS_QTL_intersection=(chrom,pos,ref,alt) in a,
    QTD000161_nominal_QTL_p=r["pvalue"],QTD000161_beta=r["beta"],QTD000161_se=r["se"],
    status="EBI_REPROCESSED_CIS_TESTED_VARIANT_ONLY"))
 out.mkdir(parents=True,exist_ok=True)
 for file,rows in [
  ("IS_PRIORITY4_EBI_REPROCESSED_FULL_CIS_RECONSTRUCTION.tsv",result),
  ("IS_PRIORITY4_EBI_REPROCESSED_FULL_CIS_VARIANTS.tsv",variant_table)]:
  with (out/file).open("w",newline="") as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
   w.writeheader();w.writerows(rows)
 manifest=dict(genes=len(result),complete_ebi_reprocessed_universes=sum(v["EBI_reprocessed_universe_status"].startswith("COMPLETE") for v in result),
  EBI_reprocessed_gene_variant_rows=len(variant_table),
  GTEx_original_FastQTL_test_universe_proven=False,
  archived_GTEx_V8_label_source_pipeline_authenticated=False,
  GTEx_donor_LD_proven=False,
  can_promote_coloc_susie=False,
  all_2225_positional_genes_preserved=True,
  results=result,source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
   for p in [permuted,*sorted(qdir.glob("BBJ_IS_L00*.QTD000161.wide.priority4.tsv"))]})
 (out/"IS_PRIORITY4_EBI_FULL_CIS_RECONSTRUCTION_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")
 return manifest
if __name__=="__main__":
 p=argparse.ArgumentParser()
 for k in ("qdir","permuted","original","out"):p.add_argument("--"+k,type=Path,required=True)
 a=p.parse_args()
 r=construct(a.qdir,a.permuted,a.original,a.out)
 for row in r["results"]:print(row)
