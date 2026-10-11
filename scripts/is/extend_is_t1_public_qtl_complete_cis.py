#!/usr/bin/env python3
"""Read-only G0015 public eQTL extended-cis recovery against official n_variants.

For each of three EBI datasets:
* prove new GRCh38 tabix extract is a MULTISET SUPERSET of historical interval
  (including duplicated rsID rows and original beta/SE/p).
* de-duplicate by genomic allele-set per gene and compare to official
  *.permuted.tsv.gz n_variants for the same gene/tissue.
* do not infer SNP identity from count alone or run new coloc.
"""
import argparse,csv,gzip,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
EXPECTED_WINDOWS={
 "IS_XDATA_G0015":("10:101403919-105460593",71),
 "IS_XDATA_G0023":("13:108801501-111807358",14),
 "IS_XDATA_G0007":("4:78736549-81763187",10),
}
REGION="IS_XDATA_G0015"
DATASETS={"QTD000131":"QTS000015","QTD000171":"QTS000015","QTD000434":"QTS000025"}
COORDS="10:101403919-105460593"
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for b in iter(lambda:f.read(4194304),b""):h.update(b)
 return h.hexdigest()
def rd(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def save(p,rows):
 with p.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
def old_lines(path):
 cnt=Counter()
 with gzip.open(path,"rt") as f:
  for line in f:
   if line.strip():cnt[line.rstrip("\r\n")]+=1
 return cnt
def official_gene_variants(path,selected):
 genes={}
 with gzip.open(path,"rt") as f:
  reader=csv.DictReader(f,delimiter="\t")
  for row in reader:
   gene=row["molecular_trait_object_id"].split(".")[0]
   if gene not in selected:continue
   if gene in genes:raise ValueError("Duplicate official molecular trait")
   genes[gene]=int(row["n_variants"])
 return genes
def collect(oldfile,widefile,selected):
 chrom,range_s=COORDS.split(":")
 start,end=map(int,range_s.split("-"))
 pending=old_lines(oldfile)
 orig_total=sum(pending.values())
 variants=defaultdict(set);record_counts=Counter();duplicate_rows=Counter()
 new_total=0
 with gzip.open(widefile,"rt") as f:
  for ln in f:
   if not ln.strip():continue
   new_total+=1
   line=ln.rstrip("\r\n")
   if pending[line]>0:pending[line]-=1
   v=line.split("\t")
   if len(v)!=19:raise ValueError("Nominal QTL source schema drift")
   if v[1]!=chrom or not start<=int(v[2])<=end:
    raise ValueError("Outside extended source region")
   if v[5]!=f"chr{v[1]}_{v[2]}_{v[3]}_{v[4]}":
    raise ValueError("Noncanonical variant id in published QTL")
   gene=v[16].split(".")[0]
   if gene not in selected:continue
   record_counts[gene]+=1
   snp=(v[1],int(v[2]),v[3],v[4])
   if snp in variants[gene]:duplicate_rows[gene]+=1
   else:variants[gene].add(snp)
 if sum(pending.values()):raise ValueError("Extended QTL region missing "+str(sum(pending.values()))+" ORIGINAL ROWS")
 if new_total<orig_total:raise ValueError("Extended file smaller than historic region")
 return variants,record_counts,duplicate_rows,orig_total,new_total
def audit(source,genes,permd,oldroot,newroot):
 if len(source)!=15 or len(genes)!=2225:raise ValueError("Input positional catalogue drift")
 sources={s["dataset_id"]:s for s in source if s["region_id"]==REGION}
 if set(sources)!=set(DATASETS):raise ValueError("Expected three T1 qtl sources")
 targets={g["gene_id_stable"]:g for g in genes if g["region_id"]==REGION}
 if len(targets)!=EXPECTED_WINDOWS[REGION][1]:raise ValueError("Original selected T1 candidate count changed")
 rows=[];audits=[];sourcehash={}
 for did,study in DATASETS.items():
  src=sources[did]
  chrom,c=src["queried_region_GRCh38"].split(":");start,end=c.split("-")
  old=oldroot/f"{REGION}__{did}__GRCh38_{chrom}_{start}_{end}.tsv.gz"
  w_chr,w_range=COORDS.split(":")
  w_start,w_end=w_range.split("-")
  wide=newroot/f"{REGION}__{did}__GRCh38_{w_chr}_{w_start}_{w_end}.tsv.gz"
  perm=permd/f"{did}.permuted.tsv.gz"
  for p in (old,wide,perm):
   if not p.is_file():raise FileNotFoundError(p)
  if sha(old)!=src["raw_sha256"]:raise ValueError("Historical raw source checksum mismatch")
  for p in (old,wide,perm):sourcehash[str(p)]=sha(p)
  n=official_gene_variants(perm,set(targets))
  variants,records,dups,oldtotal,widetotal=collect(old,wide,set(targets))
  if oldtotal!=int(src["regional_qtl_rows"]):
   raise ValueError("Historical all gene source row count changed")
  status=Counter()
  for gene,g in sorted(targets.items()):
   m=n.get(gene)
   count=len(variants.get(gene,()))
   category=("NO_PERMUTED_GENE_TEST_ON_OFFICIAL_SOURCE" if m is None else
      "PERMUTED_GENE_PRESENT_BUT_NO_VARIANTS_IN_EXTENDED_WINDOW" if count==0 else
      "MATCHES_PERMUTED_FULL_CIS_VARIANT_COUNT" if count==m else
      "EXTENDED_CIS_STILL_TRUNCATED" if count<m else
      "EXTENDED_SOURCE_EXCEEDS_OFFICIAL_CIS_VARIANT_COUNT")
   status[category]+=1
   rows.append({
    "region_id":REGION,"QTL_dataset":did,"gene_id":gene,
    "gene_symbol":g["canonical_gene_symbol"],
    "historical_ABF_tested":int(g["archived_646_coloc_assays"])>0,
    "extended_region_GRCh38":COORDS,
    "source_nominal_QTL_records":records[gene],
    "distinct_extended_cis_variant_ids":count,
    "extended_QTL_duplicate_rsid_rows":dups[gene],
    "full_cis_official_ebi_permuted_n_variants":m if m is not None else "",
    "fraction_full_cis_count_covered":count/m if m is not None and m>0 else "",
    "cis_completeness_status":category,
    "original_public_region_is_subset_multiset_verified":True,
    "full_cis_variant_ID_identity_proven_from_count_alone":False,
    "GWAS_in_original_2M_window_only":True,
    "can_claim_coloc_result":False
   })
  audits.append({
   "region_id":REGION,"QTL_dataset":did,"original_region_GRCh38":src["queried_region_GRCh38"],
   "extended_region_GRCh38":COORDS,"source_old_all_gene_rows":oldtotal,
   "source_extended_all_gene_rows":widetotal,
   "n_71_positional_genes":len(targets),
   "n_official_permuted_gene_tests":len(n),
   "n_extended_cis_count_equal_official":status["MATCHES_PERMUTED_FULL_CIS_VARIANT_COUNT"],
   "n_still_regional_partial":status["EXTENDED_CIS_STILL_TRUNCATED"],
   "n_untested_on_official_permutation":status["NO_PERMUTED_GENE_TEST_ON_OFFICIAL_SOURCE"],
   "n_source_exceeds_permutation":status["EXTENDED_SOURCE_EXCEEDS_OFFICIAL_CIS_VARIANT_COUNT"],
   "source_validation":"FULL_MULTIDATASET_ORIGINAL_ROW_MULTISET_CONTAINMENT_AND_OFFICIAL_CIS_N_VARIANTS"
  })
  if status["EXTENDED_SOURCE_EXCEEDS_OFFICIAL_CIS_VARIANT_COUNT"]>0:
   raise ValueError("Source-wide cis count exceeds authoritative permuted n_variants")
  print("T1_EXTENDED_CIS_SOURCE_PASS",REGION,did,"wide",widetotal,
    "perm_genes",len(n),"complete_count",status["MATCHES_PERMUTED_FULL_CIS_VARIANT_COUNT"],
    "still_partial",status["EXTENDED_CIS_STILL_TRUNCATED"],flush=True)
 return rows,audits,sourcehash
def main():
 p=argparse.ArgumentParser()
 for key in ("source-manifest","gene-readiness","permuted-dir","old-root","wide-root","out-dir"):
  p.add_argument("--"+key,type=Path,required=True)
 p.add_argument("--region",required=True,choices=tuple(EXPECTED_WINDOWS))
 a=p.parse_args()
 global REGION,COORDS
 REGION=a.region
 COORDS=EXPECTED_WINDOWS[REGION][0]
 if a.out_dir.exists() and any(a.out_dir.iterdir()):raise FileExistsError("No output overwrite")
 result,summary,sh=audit(rd(a.source_manifest),rd(a.gene_readiness),
                        a.permuted_dir,a.old_root,a.wide_root)
 a.out_dir.mkdir(parents=True,exist_ok=True)
 save(a.out_dir/f"{REGION}_EXTENDED_CIS_GENE_QTL_UNIVERSE.tsv",result)
 save(a.out_dir/f"{REGION}_THREE_PUBLIC_QTL_CIS_RECOVERY_SUMMARY.tsv",summary)
 status=Counter(x["cis_completeness_status"] for x in result)
 manifest={
  "status":"G0015_CIS_SOURCE_COUNT_COMPLETENESS_VS_EBI_PERMUTED_NOT_COLOC",
  "n_target_genes":EXPECTED_WINDOWS[REGION][1],"n_gene_dataset_tests":len(result),
  "n_datasets":len(summary),"count_statuses":dict(status),
  "full_cis_count_equal_does_not_prove_variant_universe_identity":True,
  "original_region_is_exact_multiset_subset_of_extended":True,
  "joint_GWAS_QTL_SNP_harmonization_for_wide_region_not_done":True,
  "QTL_study_genotype_LD_and_case_control_ratio_not_verified":True,
  "no_new_causal_gene_promotion":True,
  "SHA256_inputs":sh
 }
 (a.out_dir/f"{REGION}_THREE_PUBLIC_QTL_EXTENDED_CIS_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")
 print("T1_CIS_SOURCE_RECOVERY_FINISHED",REGION,json.dumps({
  "n":len(result),"status":dict(status)}),flush=True)
if __name__=="__main__":main()
