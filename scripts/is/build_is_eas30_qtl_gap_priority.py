#!/usr/bin/env python3
"""Freeze 2225 positional genes and map the true East-Asian QTL test deficit.

EAS positional genes are exactly the 30 IS_XDATA loci (869 genes). The 43
archived GTEx coloc genes are only four-locus experiments; 826 other EAS
positional genes have not been screened. This is NOT evidence of no eQTL.
1356 genes outside EAS regions need distinct population/outcome pipelines.
"""
import argparse,csv,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
def read(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def write(p,rows):
 with p.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
  w.writeheader();w.writerows(rows)
def checksum(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for b in iter(lambda:f.read(4194304),b""):h.update(b)
 return h.hexdigest()
def audit(genes):
 if len(genes)!=2225 or len({r["gene_id_stable"] for r in genes})!=2225:
  raise ValueError("2225 unique positional stable IDs required")
 byregion=defaultdict(list);bygroup=Counter();out=[]
 for g in genes:
  n_eas=int(g["n_EAS_Japanese_regions"]);n_eur=int(g["n_EUR_regions"])
  n_qtl=int(g["new_DIRECT_prior_grid_assays"])
  regions=g["region_ids"].split(";")
  xdata=[x for x in regions if x.startswith("IS_XDATA_G")]
  if (n_eas==1)!=(len(xdata)==1) or n_eas not in (0,1):
   raise ValueError("EAS positional region links inconsistent")
  if n_qtl>0 and n_eas!=1:raise ValueError("Unlinked molecular screen")
  group=("EAS_TESTED_GTEX_ABF" if n_qtl>0 else
         "EAS_REGION_QTL_NOT_TESTED" if n_eas else
         "NON_EAS_REGION_OTHER_PIPELINE")
  bygroup[group]+=1
  if xdata:byregion[xdata[0]].append(g)
  out.append({
   "gene_id_stable":g["gene_id_stable"],
   "canonical_gene_symbol":g["canonical_GENCODE_v19_gene_symbol"],
   "EAS_region_id":xdata[0] if xdata else "",
   "number_EAS_positional_regions":n_eas,
   "number_EUR_positional_regions":n_eur,
   "archived_BBJ_GTEx_ABF_assays":n_qtl,
   "EAS_30_region_QTL_screen_status":group,
   "EAS_missing_QTL_is_negative_evidence":False,
   "scientific_causal_status":"NOT_ESTABLISHED"
  })
 if len(byregion)!=30 or sum(len(v) for v in byregion.values())!=869:
  raise ValueError("Expected 30 East Asian locus regions and 869 genes")
 if bygroup!={
    "EAS_TESTED_GTEX_ABF":43,
    "EAS_REGION_QTL_NOT_TESTED":826,
    "NON_EAS_REGION_OTHER_PIPELINE":1356
  }:raise ValueError("EAS and non-EAS candidate coverage count drift")
 region_rows=[]
 for name,items in byregion.items():
  tested=[x for x in items if int(x["new_DIRECT_prior_grid_assays"])>0]
  n_all=len(items);miss=n_all-len(tested)
  region_rows.append({
   "EAS_region_id":name,
   "n_positional_genes":n_all,
   "n_genes_archived_GTEx_ABF_tested":len(tested),
   "n_genes_not_tested":miss,
   "n_archived_gene_tissue_coloc_assays":sum(int(x["new_DIRECT_prior_grid_assays"]) for x in items),
   "fraction_genes_tested":len(tested)/n_all,
   "scope_gap_class":"ZERO_MOLECULAR_SCREEN_IN_ARCHIVE" if not tested else
                      "PARTIAL_MOLECULAR_SCREEN_IN_ARCHIVE",
   "proposed_science_step":"LOCUS_WIDE_EAS_OUTCOME_QTL_HARMONIZATION_AND_SINGLE_OR_MULTISIGNAL_RECHECK",
   "disposition":"EXPLORATORY_NOT_CAUSAL"
  })
 region_rows.sort(key=lambda x:(-x["n_genes_not_tested"],x["EAS_region_id"]))
 if sum(x["n_genes_archived_GTEx_ABF_tested"] for x in region_rows)!=43:
  raise ValueError("43 tested genes do not map uniquely")
 if sum(x["n_archived_gene_tissue_coloc_assays"] for x in region_rows)!=646:
  raise ValueError("646 original gene-tissue assays do not map uniquely")
 return out,region_rows,dict(bygroup)
def main():
 p=argparse.ArgumentParser()
 p.add_argument("--gene-universe",type=Path,required=True)
 p.add_argument("--out-dir",type=Path,required=True)
 a=p.parse_args()
 if not a.gene_universe.is_file():raise FileNotFoundError(a.gene_universe)
 if a.out_dir.exists() and any(a.out_dir.iterdir()):raise FileExistsError("No overwrite")
 allgenes,region,groups=audit(read(a.gene_universe))
 a.out_dir.mkdir(parents=True,exist_ok=True)
 write(a.out_dir/"IS_2225_EAS_QTL_EVIDENCE_COVERAGE.tsv",allgenes)
 write(a.out_dir/"IS_30_EAS_REGIONS_QTL_GAP_PRIORITY.tsv",region)
 s={
  "status":"EXPLORATORY_30_EAS_REGION_2225_GENE_COVERAGE_MAPPING",
  "all_positional_genes_preserved":len(allgenes),
  "n_EAS_positional_regions":len(region),
  "n_EAS_positional_genes":869,
  "n_EAS_archived_QTL_tested_genes":43,
  "n_EAS_QTL_untested_genes":826,
  "n_nonEAS_region_genes":1356,
  "n_EAS_regions_with_archived_QTL_gene_tests":sum(x["n_genes_archived_GTEx_ABF_tested"]>0 for x in region),
  "n_EAS_regions_with_no_archived_QTL_gene_tests":sum(x["n_genes_archived_GTEx_ABF_tested"]==0 for x in region),
  "n_EUR_positional_genes_overlapping_or_independent":sum(int(x["number_EUR_positional_regions"])>0 for x in allgenes),
  "archive_coloc_assays":646,
  "candidate_gene_count_is_not_number_of_causal_loci":True,
  "no_missing_eQTL_test_interpreted_as_negative":True,
  "reference_source_SHA256":checksum(a.gene_universe)
 }
 (a.out_dir/"IS_EAS30_SCOPE_COVERAGE_MANIFEST.json").write_text(json.dumps(s,indent=2)+"\n")
 print("IS_EAS30_UNTESTED_826_GENE_SCOPE_AUDIT_PASS",json.dumps(s),flush=True)
if __name__=="__main__":main()
