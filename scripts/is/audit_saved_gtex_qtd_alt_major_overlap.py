#!/usr/bin/env python3
"""Exact local-only raw eQTL Catalogue GTEx_v8 QTD vs archived 646-coloc comparison.

This is a presence/provenance audit of four already downloaded GRCh38 QTD
regional files (artery aorta, coronary, tibial, brain cortex) near chr12
ALDH2 region against recovered BBJ_IS_L003 gene-tissue coloc inputs.
It DOES NOT assume that missing from coloc indicates an erroneous filter:
GWAS coverage, trait/gene inclusion, locus window or rsID QC can explain it.
Never download or modify original sources.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json
from collections import Counter
from pathlib import Path

SOURCES={
 "Artery_Aorta":"QTD000131_artery_aorta_GRCh38_chr12_111650000_112370000.tsv.gz",
 "Artery_Coronary":"QTD000136_artery_coronary_GRCh38_chr12_111650000_112370000.tsv.gz",
 "Artery_Tibial":"QTD000141_artery_tibial_GRCh38_chr12_111650000_112370000.tsv.gz",
 "Brain_Cortex":"QTD000171_brain_cortex_GRCh38_chr12_111650000_112370000.tsv.gz",
}
def hfile(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for block in iter(lambda:f.read(4*1024*1024),b""):h.update(block)
 return h.hexdigest()

def raw_qtd_iter(path):
 # Source is a tabix region extraction using 19 columns from eQTL Catalogue
 # nominal QTD files. No header exists in the region-extracted sample.
 with gzip.open(path,"rt",newline="") as f:
  for row in csv.reader(f,delimiter="\t"):
   if len(row)!=19:raise ValueError("Unexpected 19-column QTD source shape")
   if row[1]!="12":raise ValueError("Wrong QTD locus chromosome")
   ac=int(float(row[12]));an=int(float(row[13]))
   if not (an>0 and 0<=ac<=an):raise ValueError("AC/AN outside expected range")
   maf=float(row[7])
   if abs(maf-min(ac/an,1-ac/an))>1e-4:raise ValueError("Original QTD MAF/AC/AN incompatible")
   variant=row[5]
   if not variant.startswith("chr12_"):raise ValueError("QTD variant is not chr12")
   is_snv=(len(row[3])==1 and len(row[4])==1 and row[3] in "ACGT" and row[4] in "ACGT")
   yield {
    "is_snv":is_snv,
    "gene_id":row[0].split(".")[0],
    "raw_variant":variant,
    "site_GRCh38":row[1]+":"+row[2]+":"+row[3]+":"+row[4],
    "ALT_major_GTEx_ac_over_an":ac/an>0.5,
    "GTEx_ALT_ac_an":ac/an,
    "original_QTD_MAF":maf,
    "QTD_beta":row[9],
    "QTD_se":row[10],
    "QTD_p":row[8],
   }

def audit(root,cache):
 observed={}
 all_locus_sites=set()
 cache_paths=list(cache.glob("BBJ_IS_L003__GTEx_V8__*__*.tsv"))
 if not cache_paths:raise ValueError("L003 recovered coloc files missing")
 for f in cache_paths:
  if not f.name.startswith("BBJ_IS_L003__GTEx_V8__") or not f.name.endswith(".tsv"):
   raise ValueError("Bad expected cached file syntax")
  with f.open(newline="") as handle:
   for r in csv.DictReader(handle,delimiter="\t"):
    all_locus_sites.add(r["match_key"])
 for tissue,fn in SOURCES.items():
  pattern="BBJ_IS_L003__GTEx_V8__"+tissue+"__*.tsv"
  paths=list(cache.glob(pattern))
  if not paths:raise ValueError("Missing L003 tissue "+tissue)
  genes=set()
  exact=set()
  for p in paths:
   with p.open(newline="") as f:
    for r in csv.DictReader(f,delimiter="\t"):
     genes.add(r["gene_base"])
     exact.add((r["gene_base"],r["variant"]))
  # Records are duplicated across gene tests. Keep both row count and variant identities.
  counts=Counter()
  unique_maj=set();major_rows=[]
  n_match_exact=0
  for row in raw_qtd_iter(root/fn):
   counts["raw_rows"]+=1
   if not row["is_snv"]:
    counts["original_QTD_non_SNV_rows_out_of_scope"]+=1
    continue
   key=(row["gene_id"],row["raw_variant"])
   if row["ALT_major_GTEx_ac_over_an"]:
    counts["raw_ALT_major_rows"]+=1
   if row["gene_id"] not in genes:
    counts["raw_gene_not_in_archived_L003_tissue"]+=1
    continue
   counts["raw_rows_eligible_gene"]+=1
   in_cache=key in exact
   in_locus=row["site_GRCh38"] in all_locus_sites
   counts["exact_raw_gene_variant_in_cached_tissue"]+=in_cache
   counts["raw_gene_variant_not_in_cached_tissue"]+=not in_cache
   if row["ALT_major_GTEx_ac_over_an"]:
    counts["ALT_major_rows_with_cached_gene"]+=1
    if in_cache:counts["ALT_major_rows_exactly_in_cached_tissue"]+=1
    if in_locus:counts["ALT_major_rows_site_present_in_any_L003_tissue"]+=1
    unique_maj.add(row["site_GRCh38"])
    major_rows.append({
      "source_file":fn,"tissue":tissue,"gene_id":row["gene_id"],
      "variant":row["raw_variant"],"site_GRCh38":row["site_GRCh38"],
      "QTD_ALT_ac_an":row["GTEx_ALT_ac_an"],
      "QTD_minor_maf":row["original_QTD_MAF"],
      "same_gene_variant_in_archived_L003_tissue":in_cache,
      "same_site_in_any_archived_L003_tissue":in_locus,
      "reason_not_in_cache":"GWAS_COVERAGE_OR_ANALYSIS_SELECTION_UNRESOLVED" if not in_cache else "PRESENT",
    })
  if counts["raw_rows"]<1000:raise ValueError("Original saved QTD source unexpectedly small")
  observed[tissue]={"file":fn,"archived_files":len(paths),
      "unique_cached_gene_variants":len(exact),
      **dict(counts),"unique_GTEx_ALT_major_sites_in_archived_gene_panel":len(unique_maj),
      "major_details":major_rows}
 return observed

def main():
 p=argparse.ArgumentParser()
 p.add_argument("--original-qtd-dir",required=True,type=Path)
 p.add_argument("--recovered-coloc-dir",required=True,type=Path)
 p.add_argument("--out-dir",required=True,type=Path)
 a=p.parse_args()
 if any(not (a.original_qtd_dir/fn).is_file() for fn in SOURCES.values()):
  raise FileNotFoundError("One or more existing original QTD source files are absent")
 report=audit(a.original_qtd_dir,a.recovered_coloc_dir)
 if a.out_dir.exists() and any(a.out_dir.iterdir()):
  raise FileExistsError("No overwrite of prior audit")
 a.out_dir.mkdir(parents=True,exist_ok=True)
 summary=[];details=[]
 for tissue,d in report.items():
  details+=d.pop("major_details")
  summary.append({"tissue":tissue,**d})
 with (a.out_dir/"IS_GTEX_QTD_L003_TISSUE_ORIGINAL_VS_COLOC.tsv").open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(summary[0]),delimiter="\t")
  w.writeheader();w.writerows(summary)
 with (a.out_dir/"IS_GTEX_QTD_ALT_MAJOR_ROW_GENE_MATCH.tsv").open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(details[0]),delimiter="\t")
  w.writeheader();w.writerows(details)
 manifest={
  "status":"ORIGINAL_SAVED_QTD_VERSUS_ARCHIVED_COLOC_REGION_SELECTION_AUDIT",
  "original_QTD_tissue_files":len(summary),
  "QTD_rows":sum(x["raw_rows"] for x in summary),
  "original_ALT_major_rows":sum(x["raw_ALT_major_rows"] for x in summary),
  "ALT_major_rows_in_archived_genes":sum(x["ALT_major_rows_with_cached_gene"] for x in summary),
  "ALT_major_gene_variant_rows_found_in_coloc":sum(x["ALT_major_rows_exactly_in_cached_tissue"] for x in summary),
  "ALT_major_rows_whose_site_was_in_any_coloc_L003_tissue":sum(x["ALT_major_rows_site_present_in_any_L003_tissue"] for x in summary),
  "original_QTD_SHA256":{fn:hfile(a.original_qtd_dir/fn) for fn in SOURCES.values()},
  "effect_allele_orientation":"OFFICIAL_EQTL_CATALOGUE_ALT_EFFECT_CONVENTION; RAW_GENOTYPE_NOT_RECHECKED",
  "cause_of_discrepant_variant_coverage":"UNKNOWN: GWAS availability, cis window, original processing, or QTL variant filters",
  "no_data_rediscovered_or_downloaded":True,
  "automatic_allele_flips_or_variant_exclusion":False,
  "causal_gene_promotion":False,
 }
 (a.out_dir/"IS_GTEX_QTD_SAVED_SOURCE_COVERAGE_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")
 print("IS_GTEX_QTD_SOURCE_AUDIT_COMPLETE",json.dumps(manifest),flush=True)
if __name__=="__main__":main()
