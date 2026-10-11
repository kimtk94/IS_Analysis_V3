#!/usr/bin/env python3
"""Audit complete cis-tested universe vs regional eQTL Catalogue slices.

The official *.permuted.tsv.gz n_variants is the count tested for a
molecular trait's COMPLETE cis window in its study. A preexisting
GRCh38 genomic region slice is not automatically that cis universe.
Read-only, no causal/coloc claims; provenance SHA256.
"""
import argparse,csv,gzip,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path

QTD={"QTD000131":"GTEx_Artery_Aorta","QTD000171":"GTEx_Brain_Cortex",
     "QTD000434":"ROSMAP_DLPFC"}
REGIONS=("IS_XDATA_G0015","IS_XDATA_G0023","IS_XDATA_G0007",
         "IS_XDATA_G0025","IS_XDATA_G0018")
def rd(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def write(p,rs):
 if not rs:raise ValueError("Cannot save empty audited data")
 with p.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter="\t")
  w.writeheader();w.writerows(rs)
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for b in iter(lambda:f.read(4194304),b""):h.update(b)
 return h.hexdigest()
def parse_permuted(path,target_genes):
 if not path.is_file():raise FileNotFoundError(path)
 obj={};all_records=0
 with gzip.open(path,"rt") as f:
  reader=csv.DictReader(f,delimiter="\t")
  if not {"molecular_trait_object_id","n_traits","n_variants","p_beta","p_perm"}.issubset(reader.fieldnames):
   raise ValueError("Official permutations required schema missing")
  for r in reader:
   all_records+=1
   gene=r["molecular_trait_object_id"].split(".")[0]
   if gene not in target_genes:continue
   n=int(r["n_variants"])
   if n<0:raise ValueError("Negative n_variants")
   p=float(r["p_beta"]) if r["p_beta"]!="NA" else None
   if p is not None and not (math.isfinite(p) and 0<=p<=1):
    raise ValueError("Invalid permutation p_beta")
   if gene in obj:raise ValueError("Ambiguous duplicate molecular trait object ID")
   obj[gene]=dict(n_variants=n,p_beta=p,n_traits=int(r["n_traits"]))
 return obj,all_records
def qtl_path(source,t1,tech):
 rid,did=source["region_id"],source["dataset_id"]
 chr_,positions=source["queried_region_GRCh38"].split(":")
 start,end=positions.split("-")
 root=t1 if rid in REGIONS[:3] else tech
 return root/f"{rid}__{did}__GRCh38_{chr_}_{start}_{end}.tsv.gz"
def region_variants(path,target,chrom):
 unique=defaultdict(set);nrows=Counter();duplicate_rows=Counter()
 invalid=Counter();all_n=0
 with gzip.open(path,"rt") as f:
  for ln in f:
   if not ln.strip():continue
   all_n+=1
   c=ln.rstrip("\r\n").split("\t")
   if len(c)!=19:raise ValueError("Unexpected 19-column nominal QTL schema")
   gene=c[16].split(".")[0]
   if gene not in target:continue
   nrows[gene]+=1
   if c[1]!=chrom or c[5]!=f"chr{chrom}_{c[2]}_{c[3]}_{c[4]}":
    raise ValueError("Nominal QTL SNP ID/coordinate/alleles mismatch")
   ident=(c[1],int(c[2]),c[3].upper(),c[4].upper())
   if ident in unique[gene]:duplicate_rows[gene]+=1
   else:unique[gene].add(ident)
 return unique,nrows,duplicate_rows,all_n,invalid
def infer_regional_status(full,unique):
 if full is None:return "NO_EBI_PERMUTED_GENE_RECORD_CANNOT_ASSESS"
 n=full["n_variants"]
 if unique==0:return "PERMUTED_GENE_TESTED_BUT_NO_NOMINAL_ROWS_IN_REGIONAL_SLICE"
 if n==0:return "NVARIANTS_ZERO_BUT_REGIONAL_ROWS_PRESENT_INCONSISTENT"
 if unique==n:return "FULL_CIS_EQUIVALENCE_VS_EBI_NVARIANTS_COUNT_ONLY"
 if unique<n:return "REGIONAL_SLICE_INCOMPLETE_CIS_NOT_SAFE_FOR_FULL_LOCUS_COLOC"
 return "REGIONAL_ROWS_EXCEED_PERMUTED_NVARIANTS_INCONSISTENT"
def audit(source,genecov,readiness,permd,t1,tech):
 if len(source)!=15 or len(genecov)!=720 or len(readiness)!=2225:
  raise ValueError("Frozen 15 slices, 720 gene-dataset pairs, 2225 gene IDs required")
 bysource={(x["region_id"],x["dataset_id"]):x for x in source}
 by_gene={(x["region_id"],x["dataset_id"],x["gene_id"]):x for x in genecov}
 original_genes={r["gene_id_stable"]:r for r in readiness}
 if len(bysource)!=15 or len(by_gene)!=720 or len(original_genes)!=2225:raise ValueError("Duplicate source or gene")
 by_target=defaultdict(set)
 for r in genecov:
  if r["region_id"] not in REGIONS or r["dataset_id"] not in QTD:raise ValueError("Unknown scope")
  if r["gene_id"] not in original_genes:raise ValueError("Original candidate ID missing")
  by_target[r["dataset_id"]].add(r["gene_id"])
 perm={};perm_counts={};source_sha={}
 for did,genes in by_target.items():
  p=permd/f"{did}.permuted.tsv.gz"
  perm[did],perm_counts[did]=parse_permuted(p,genes)
  source_sha[str(p)]=sha(p)
  print("PERMUTED_SOURCE",did,"all_records",perm_counts[did],
        "target_genes_with_perm",len(perm[did]),flush=True)
 rows=[]
 per_source=[]
 for (rid,did),src in sorted(bysource.items()):
  if src["status"]!="PUBLIC_NOMINAL_QTL_REGION_DOWNLOADED":
   raise ValueError("Incomplete source QTL input")
  targets={g for (r,t,g) in by_gene if (r,t)==(rid,did)}
  raw=qtl_path(src,t1,tech)
  if not raw.is_file() or sha(raw)!=src["raw_sha256"]:
   raise ValueError("QTL nominal source SHA256 mismatch")
  source_sha[str(raw)]=src["raw_sha256"]
  chr38=src["queried_region_GRCh38"].split(":")[0]
  unique,nrows,duplicates,total,invalid=region_variants(raw,targets,chr38)
  if total!=int(src["regional_qtl_rows"]):
   raise ValueError("Raw regional gene-SNP row count changed")
  st=Counter();tested_unique=0;truncated_unique=0
  for gene in sorted(targets):
   ref=by_gene[(rid,did,gene)]
   if nrows[gene]!=int(ref["n_raw_QTL_records"]):
    raise ValueError("Original nominal gene record count changed")
   observed=len(unique[gene])
   pm=perm[did].get(gene)
   label=infer_regional_status(pm,observed)
   st[label]+=1
   if pm is not None:tested_unique+=1
   if pm is not None and observed<pm["n_variants"]:truncated_unique+=1
   rows.append({
    "region_id":rid,"QTL_dataset":did,"QTL_study":QTD[did],
    "gene_id":gene,
    "gene_symbol":original_genes[gene]["canonical_gene_symbol"],
    "archived_ABF_assays":int(original_genes[gene]["archived_646_coloc_assays"]),
    "nominal_regional_raw_rows":nrows[gene],
    "nominal_regional_distinct_SNP_INDEL_variants":observed,
    "nominal_regional_rsID_duplicate_rows":duplicates[gene],
    "ebi_permuted_trait_observed":pm is not None,
    "ebi_permuted_full_cis_n_variants":pm["n_variants"] if pm is not None else "",
    "ebi_permuted_n_traits":pm["n_traits"] if pm is not None else "",
    "ebi_permuted_gene_p_beta":pm["p_beta"] if pm is not None else "",
    "fraction_ebi_cis_variant_count_in_region":observed/pm["n_variants"] if pm is not None and pm["n_variants"]>0 else "",
    "regional_nominal_source_is_full_cis":label=="FULL_CIS_EQUIVALENCE_VS_EBI_NVARIANTS_COUNT_ONLY",
    "ebi_source_test_status":label,
    "matched_QTL_donor_ancestry_LD_verified":False,
    "SNPwise_GWAS_overlap_from_prior_audit_not_coloc":True,
    "colocalization_posteriors_recomputed":False,
    "scientific_status":"COMPLETE_CIS_GATE_SEPARATE_FROM_CAUSAL_GENE_INFERENCE"
   })
  per_source.append({
   "region_id":rid,"QTL_dataset":did,"n_positional_genes":len(targets),
   "ebi_permuted_genes_with_test_records":tested_unique,
   "genes_with_no_permuted_records":st["NO_EBI_PERMUTED_GENE_RECORD_CANNOT_ASSESS"],
   "genes_with_fully_observed_cis_by_n_variants_COUNT_ONLY":st["FULL_CIS_EQUIVALENCE_VS_EBI_NVARIANTS_COUNT_ONLY"],
   "genes_with_regional_partial_cis":st["REGIONAL_SLICE_INCOMPLETE_CIS_NOT_SAFE_FOR_FULL_LOCUS_COLOC"],
   "genes_with_zero_nominal_rows_but_PERMUTED_TESTED":st["PERMUTED_GENE_TESTED_BUT_NO_NOMINAL_ROWS_IN_REGIONAL_SLICE"],
   "genes_with_count_inconsistency":st["REGIONAL_ROWS_EXCEED_PERMUTED_NVARIANTS_INCONSISTENT"]+st["NVARIANTS_ZERO_BUT_REGIONAL_ROWS_PRESENT_INCONSISTENT"],
   "regional_rows":total,
   "target_gene_nominal_rows":sum(nrows.values()),
   "source_verification":"SHA256_MATCH_PERMUTED_INDEPENDENT_QTL_CIS_TESTED_COUNT"
  })
  print("REGION_QTL_CIS_AUDIT",rid,did,"tested_perm",tested_unique,
        "complete_count",st["FULL_CIS_EQUIVALENCE_VS_EBI_NVARIANTS_COUNT_ONLY"],
        "partial",st["REGIONAL_SLICE_INCOMPLETE_CIS_NOT_SAFE_FOR_FULL_LOCUS_COLOC"],
        "count_inconsistent",st["REGIONAL_ROWS_EXCEED_PERMUTED_NVARIANTS_INCONSISTENT"],
        flush=True)
 if len(rows)!=720:raise ValueError("Candidate tests changed")
 return rows,per_source,perm_counts,source_sha
def main():
 p=argparse.ArgumentParser()
 for k in ("source-manifest","source-gene-cov","gene-readiness","permuted-dir","t1-raw","technical-raw","out-dir"):
  p.add_argument("--"+k,required=True,type=Path)
 a=p.parse_args()
 files=(a.source_manifest,a.source_gene_cov,a.gene_readiness)
 for z in files:
  if not z.is_file():raise FileNotFoundError(z)
 if a.out_dir.exists() and any(a.out_dir.iterdir()):raise FileExistsError("Refuse output overwrite")
 rows,src,n,sh=audit(rd(a.source_manifest),rd(a.source_gene_cov),
  rd(a.gene_readiness),a.permuted_dir,a.t1_raw,a.technical_raw)
 a.out_dir.mkdir(parents=True,exist_ok=True)
 write(a.out_dir/"IS_EAS5_720_GENE_QTL_FULL_CIS_UNIVERSE_AUDIT.tsv",rows)
 write(a.out_dir/"IS_EAS5_15_DATASETS_FULL_CIS_SOURCE_AUDIT.tsv",src)
 status=Counter(x["ebi_source_test_status"] for x in rows)
 summary={
  "status":"VERIFIED_EBI_PERMUTED_NVARIANTS_VS_GRCH38_REGIONAL_SOURCE_UNIVERSE",
  "n_original_gene_ids":2225,"n_pilot_regional_gene_ids":240,
  "gene_qtl_dataset_pairs":len(rows),"QTL_regional_source_slices":len(src),
  "n_permutation_records_scanned":n,
  "status_counts":dict(status),
  "full_cis_count_equivalence_does_not_guarantee_identical_variant_set":True,
  "sample_ancestry_molecular_LD_not_validated":True,
  "new_coloc_posteriors_not_calculated":True,
  "source_SHA256":{**sh,**{str(z):sha(z) for z in files}}
 }
 (a.out_dir/"IS_EAS5_PERMUTED_VS_REGIONAL_CIS_MANIFEST.json").write_text(json.dumps(summary,indent=2)+"\n")
 print("IS_EAS5_720_QTL_GENE_CIS_UNIVERSE_AUDIT_COMPLETE",json.dumps({
  "total_tests":len(rows),"status_counts":dict(status)}),flush=True)
if __name__=="__main__":main()
