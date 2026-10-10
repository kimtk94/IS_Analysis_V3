#!/usr/bin/env python3
"""Consolidate 15 PUBLIC NOMINAL QTL region extracts over five IS_XDATA groups.

Five pilot regions: two technical/suggestive and three T1 outcome-associated;
gene-level inventory covers all original 2225 positional gene IDs. No allele
harmonization with GWAS is inferred; no coloc or disease-expression performed.
"""
import argparse,csv,gzip,hashlib,json,math
from collections import defaultdict,Counter
from pathlib import Path
EXPECT={"IS_XDATA_G0025":(92,61),"IS_XDATA_G0018":(53,29),
        "IS_XDATA_G0015":(71,46),"IS_XDATA_G0023":(14,9),"IS_XDATA_G0007":(10,4)}
DATASETS={"QTD000131","QTD000171","QTD000434"}
def read(path):
 with path.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def sha(path):
 h=hashlib.sha256()
 with path.open("rb") as f:
  for b in iter(lambda:f.read(4194304),b""):h.update(b)
 return h.hexdigest()
def write(path,rows):
 with path.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
  w.writeheader();w.writerows(rows)
def build(sources,ql,genes,priority,raw_dirs):
 if len(sources)!=15 or len(genes)!=2225 or len({x["gene_id_stable"] for x in genes})!=2225 or len(priority)!=30:
  raise ValueError("Expected 15 source slices and preserved 2225/30 tables")
 if len(ql)!=sum(v[0]*3 for v in EXPECT.values()):
  raise ValueError("Expected all 720 gene/dataset source coverage rows")
 group=defaultdict(list);qgenes=defaultdict(dict)
 for s in sources:
  rid,did=s["region_id"],s["dataset_id"]
  if rid not in EXPECT or did not in DATASETS or s["status"]!="PUBLIC_NOMINAL_QTL_REGION_DOWNLOADED":
   raise ValueError("Nonwhitelisted or incomplete QTL source")
  if int(s["regional_qtl_rows"])<=0:raise ValueError("Empty source rows")
  group[rid].append(s)
 for rid in EXPECT:
  if {s["dataset_id"] for s in group[rid]}!=DATASETS:
   raise ValueError("Expected three public nominal datasets per region")
 source_checks={}
 for s in sources:
  region,dataset=s["region_id"],s["dataset_id"]
  c,interval=s["queried_region_GRCh38"].split(":")
  start,end=interval.split("-")
  root=raw_dirs["T1"] if region in ("IS_XDATA_G0015","IS_XDATA_G0023","IS_XDATA_G0007") else raw_dirs["T4"]
  path=root/f"{region}__{dataset}__GRCh38_{c}_{start}_{end}.tsv.gz"
  if not path.is_file() or sha(path)!=s["raw_sha256"]:
   raise ValueError("Absent or changed source raw QTL file "+str(path))
  source_checks[str(path)]=s["raw_sha256"]
 gene_table={r["gene_id_stable"]:r for r in genes}
 region_label={r["region_id"]:r for r in priority}
 if len(region_label)!=30:raise ValueError("Duplicated priority region")
 for row in ql:
  rid,did,gene=row["region_id"],row["dataset_id"],row["gene_id"]
  if rid not in EXPECT or did not in DATASETS or gene not in gene_table or rid not in gene_table[gene]["region_ids"].split(";"):
   raise ValueError("QTL target gene is not original positional region gene")
  k=(rid,gene)
  if did in qgenes[k]:raise ValueError("Repeated QTL tissue for gene")
  qgenes[k][did]=row
 region_rows=[]
 for region,(expected_genes,expected_with_source) in EXPECT.items():
  members=[x for x in genes if region in x["region_ids"].split(";")]
  if len(members)!=expected_genes:raise ValueError("Gene positional scope drift: "+region)
  found=[]
  newly=[]
  for row in members:
   d=qgenes[(region,row["gene_id_stable"])]
   if set(d)!=DATASETS:raise ValueError("Missing gene-level raw QTL source row")
   if any(int(x["n_raw_QTL_records"])>0 for x in d.values()):
    found.append(row)
    if int(row["new_DIRECT_prior_grid_assays"])==0:newly.append(row)
  if len(found)!=expected_with_source:
   raise ValueError("Unexpected union nominal source gene coverage: "+region)
  region_rows.append({
   "region_id":region,
   "priority_tier":region_label[region]["primary_GWAS_EAS_BBJ_AIS_tier"],
   "n_original_positional_genes":len(members),
   "n_genes_with_nominal_QTL_source_records":len(found),
   "n_with_new_nominal_QTL_but_no_archived_coloc":len(newly),
   "n_with_no_rows_in_three_source_slices_not_negative":len(members)-len(found),
   "n_public_QTL_datasets":len(DATASETS),
   "n_source_gene_variant_rows":sum(int(x["target_gene_qtl_rows"]) for x in group[region]),
   "n_full_source_gene_variant_rows":sum(int(x["regional_qtl_rows"]) for x in group[region]),
   "science_gate":"NOMINAL_RAW_QTL_ONLY_NEEDS_GWAS_VARIANT_HARMONIZATION"
  })
 if sum(x["n_genes_with_nominal_QTL_source_records"] for x in region_rows)!=149 or sum(x["n_with_new_nominal_QTL_but_no_archived_coloc"] for x in region_rows)!=124:
  raise ValueError("Expected exactly 149 sourced genes, 124 without archived coloc")
 gene_rows=[]
 for r in genes:
  eas=[x for x in r["region_ids"].split(";") if x.startswith("IS_XDATA_")]
  if len(eas)>1:raise ValueError("Unexpected multiple EAS positional regions")
  rid=eas[0] if eas else ""
  study=qgenes.get((rid,r["gene_id_stable"]),{})
  vals=[x for x in study.values() if int(x["n_raw_QTL_records"])>0]
  mins=[]
  for x in vals:
   if x["min_QTL_p"]!="":
    num=float(x["min_QTL_p"])
    if not 0<=num<=1 or not math.isfinite(num):raise ValueError("Bad source QTL P")
    mins.append(num)
  status=("NOT_IN_EAS_30_GROUPS" if not eas else
    "PUBLIC_NOMINAL_QTL_REGIONAL_PILOT_NOT_EXECUTED" if rid not in EXPECT else
    "RAW_NOMINAL_QTL_RECORDS_AVAILABLE_NO_COLOC" if vals else
    "NO_NOMINAL_ROWS_IN_TARGETED_SLICES_NOT_NEGATIVE")
  gene_rows.append({
    "gene_id_stable":r["gene_id_stable"],
    "canonical_gene_symbol":r["canonical_GENCODE_v19_gene_symbol"],
    "region_id":rid,
    "EAS_source_GWAS_priority_tier":region_label[rid]["primary_GWAS_EAS_BBJ_AIS_tier"] if eas else "",
    "archived_646_coloc_assays":int(r["new_DIRECT_prior_grid_assays"]),
    "n_new_QTL_tissues_with_any_nominal_rows":len(vals),
    "source_nominal_QTL_ids":";".join(sorted(x["dataset_id"] for x in vals)),
    "n_new_nominal_QTL_gene_variant_rows":sum(int(x["n_raw_QTL_records"]) for x in vals),
    "min_source_nominal_QTL_p":min(mins) if mins else "",
    "new_nominal_QTL_gene_source_status":status,
    "GWAS_and_QTL_alleles_coharmonized":False,
    "matched_QTL_LD_verified":False,
    "coloc_recomputed":False,
    "ischemic_stroke_causal_gene":"NOT_ESTABLISHED"
  })
 manifest={
  "status":"FIVE_PUBLIC_NOMINAL_QTL_REGIONS_15_SLICES_2225_GENE_READINESS",
  "five_regions":sorted(EXPECT),
  "n_region_QTL_slices":len(sources),
  "n_public_QTL_tissue_sources":len(DATASETS),
  "n_region_QTL_raw_gene_variant_rows":sum(int(x["regional_qtl_rows"]) for x in sources),
  "n_original_positional_genes":len(gene_rows),
  "n_new_nominal_QTL_union_gene_ids":149,
  "n_nominal_QTL_union_gene_ids_without_historic_ABF":124,
  "n_T1_union_genes_with_nominal_QTL":sum(x["n_genes_with_nominal_QTL_source_records"] for x in region_rows if x["priority_tier"].startswith("T1")),
  "n_T1_nominal_QTL_genes_without_historic_ABF":sum(x["n_with_new_nominal_QTL_but_no_archived_coloc"] for x in region_rows if x["priority_tier"].startswith("T1")),
  "n_untested_archived_coloc_EAS_genes_still_826":826,
  "new_QTL_sources_not_genetic_colocalization":True,
  "molecular_QTL_population_ancestry_unverified":True,
  "raw_file_sha256":source_checks,
 }
 return gene_rows,region_rows,manifest
def main():
 ap=argparse.ArgumentParser()
 for k in ("technical-sources","technical-genes","T1-sources","T1-genes",
           "gene-universe","region-priority","technical-raw-dir","T1-raw-dir","out-dir"):
  ap.add_argument("--"+k,required=True,type=Path)
 a=ap.parse_args()
 paths=[a.technical_sources,a.technical_genes,a.T1_sources,a.T1_genes,a.gene_universe,a.region_priority]
 for path in paths:
  if not path.is_file():raise FileNotFoundError(path)
 if a.out_dir.exists() and any(a.out_dir.iterdir()):raise FileExistsError("No output overwrite")
 src=read(a.technical_sources)+read(a.T1_sources)
 source_genes=read(a.technical_genes)+read(a.T1_genes)
 gene_rows,region_rows,manifest=build(src,source_genes,read(a.gene_universe),
   read(a.region_priority),{"T1":a.T1_raw_dir,"T4":a.technical_raw_dir})
 a.out_dir.mkdir(parents=True,exist_ok=True)
 write(a.out_dir/"IS_15_PUBLIC_QTL_REGION_SOURCE_COVERAGE.tsv",src)
 write(a.out_dir/"IS_15_PUBLIC_QTL_GENE_SOURCE_COVERAGE.tsv",source_genes)
 write(a.out_dir/"IS_FIVE_EAS_REGIONS_NOMINAL_GENE_UNION.tsv",region_rows)
 write(a.out_dir/"IS_2225_GENES_PUBLIC_NOMINAL_QTL_READINESS.tsv",gene_rows)
 manifest["input_manifest_SHA256"]={str(x):sha(x) for x in paths}
 (a.out_dir/"IS_PUBLIC_NOMINAL_QTL_EXPANSION_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")
 print("IS_PUBLIC_NOMINAL_QTL_15_SLICES_124_NEW_SOURCE_GENE_IDS_PASS",json.dumps({
  k:manifest[k] for k in ("n_region_QTL_slices","n_region_QTL_raw_gene_variant_rows",
  "n_original_positional_genes","n_new_nominal_QTL_union_gene_ids",
  "n_nominal_QTL_union_gene_ids_without_historic_ABF",
  "n_T1_nominal_QTL_genes_without_historic_ABF")}),flush=True)
if __name__=="__main__":main()
