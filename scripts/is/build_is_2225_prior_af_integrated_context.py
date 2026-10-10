#!/usr/bin/env python3
"""Overlay 646 direct SNP-coloc prior/AF QC on all 2,225 positional gene IDs.

The exploratory ranking universe is NOT filtered. No causal assignment,
matched QTL LD claim, or double-count of the same locus/tissue/SNP is made.
Only exact (BBJ locus, stable Ensembl gene ID) keys connect QTL to GWAS.
"""
from __future__ import annotations
import argparse,csv,json,math,statistics,hashlib
from pathlib import Path
from collections import Counter,defaultdict

PRIORS=(1e-6,3e-6,1e-5,3e-5,1e-4)
KEYS=lambda r:(r["locus"],r["dataset_key"],r["gene_base"])
def read(p):
 with p.open(newline="",encoding="utf-8") as f:
  return list(csv.DictReader(f,delimiter="\t"))
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for block in iter(lambda:f.read(4*1024*1024),b""):h.update(block)
 return h.hexdigest()
def writer(path,rows):
 with path.open("w",newline="",encoding="utf-8") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
  w.writeheader();w.writerows(rows)

def integrate(genes,links,grid,af,nsens,replay):
 ids={g["gene_id_stable"] for g in genes}
 if len(genes)!=2225 or len(ids)!=2225:
  raise ValueError("Expected exactly 2225 unique stable Ensembl gene IDs")
 if len(links)!=2425 or sum(int(g["n_positional_regions"]) for g in genes)!=2425:
  raise ValueError("Expanded gene-region source changed")
 if len(grid)!=646*5 or len(af)!=646 or len(nsens)!=646 or len(replay)!=646:
  raise ValueError("Expected 646 assay × 5-prior grid and complete QC")
 gkeys={KEYS(r) for r in grid}
 if len(gkeys)!=646 or gkeys!={KEYS(r) for r in af} or gkeys!={KEYS(r) for r in nsens} or gkeys!={KEYS(r) for r in replay}:
  raise ValueError("QTL assay keys disagree across data sources")
 grouped=defaultdict(dict)
 for r in grid:
  key=KEYS(r);p=float(r["p12"])
  if p not in PRIORS or p in grouped[key]:
   raise ValueError("Duplicate/unrecognized p12 priors")
  grouped[key][p]=r
 if any(set(x)!=set(PRIORS) for x in grouped.values()):
  raise ValueError("Missing prespecified prior settings")
 af_by={KEYS(r):r for r in af}
 ns_by={KEYS(r):r for r in nsens}
 baseline={KEYS(r):r for r in replay}
 for key,r in baseline.items():
  if r["status"]!="PASS" or float(r["delta_max"])>1e-8:
   raise ValueError("New archived ABF rerun not validated for "+str(key))
  if abs(float(r["new_PP_H4"])-float(grouped[key][1e-5]["PP_H4"]))>1e-8:
   raise ValueError("p12 baseline and independently rerun H4 disagree")
 for key,x in ns_by.items():
  if abs(float(x["archived_H4"])-float(grouped[key][1e-5]["PP_H4"]))>1e-8:
   raise ValueError("QTL sample-N historic baseline differs from direct prior grid")
 assay_out=[]
 qtl_by_gene=defaultdict(list)
 for key in sorted(gkeys):
  locus,ds,gid=key
  r=grouped[key]
  ar=af_by[key]; nr=ns_by[key]
  h4=[float(r[p]["PP_H4"]) for p in PRIORS]
  h3=[float(r[p]["PP_H3"]) for p in PRIORS]
  h4_gt=[h4[i]>h3[i] for i in range(5)]
  stat={
   "locus":locus,"dataset_key":ds,"gene_base":gid,
   "PP_H4_p12_1e6":h4[0],
   "PP_H4_p12_3e6":h4[1],
   "PP_H4_p12_1e5_original":h4[2],
   "PP_H3_p12_1e5_original":h3[2],
   "PP_H4_p12_3e5":h4[3],
   "PP_H4_p12_1e4":h4[4],
   "H4_gt_H3_at_baseline":h4_gt[2],
   "H4_gt_H3_at_all_5_priors":all(h4_gt),
   "H4_ge_0p8_at_baseline":h4[2]>=0.8,
   "H4_ge_0p8_at_highest_p12":h4[4]>=0.8,
   "MAF_diff_gt_0p1_fraction":float(ar["fraction_MAF_difference_GT_0p1"]),
   "n_palindromic_recorded":int(ar["n_palindromic"]),
   "max_H4_shift_QTL_scalar_N":float(nr["max_abs_H4_shift_vs_first"]),
   "scientific_status":"PRIOR_SENSITIVE_ARCHIVED_GTEX_ABF_NOT_VALIDATED_CAUSALITY",
  }
  assay_out.append(stat)
  qtl_by_gene[gid].append(stat)
 # Every assay maps to one EAS/Japanese positional gene by stable ID + BBJ locus.
 linkset={(r["gene_id_stable"],r["legacy_original_bbj_locus"]) for r in links
          if r["ancestry"]=="EAS_AND_JAPANESE" and int(r["archived_ABF_assays"])>0}
 if any((r["gene_base"],r["locus"]) not in linkset for r in assay_out):
  raise ValueError("QTL source would be incorrectly attached to an unproven positional/EUR locus")
 merged=[]
 for g in genes:
  gid=g["gene_id_stable"];items=qtl_by_gene[gid]
  qn=int(g["archived_EAS_GTEx_ABF_assays"])
  if len(items)!=qn:
   raise ValueError("Archived GTEx assay count inconsistent for gene "+gid)
  if qn:
   p12base=max(x["PP_H4_p12_1e5_original"] for x in items)
   if abs(p12base-float(g["archived_EAS_GTEx_max_PP_H4"]))>1e-8:
    raise ValueError("One of 2225 baseline H4 values differs from original")
  else:p12base=""
  merged.append({
   **g,
   "new_DIRECT_prior_grid_assays":qn,
   "new_direct_p12_1e6_max_H4":max((x["PP_H4_p12_1e6"] for x in items),default=""),
   "new_direct_p12_1e5_baseline_max_H4":p12base,
   "new_direct_p12_1e4_max_H4":max((x["PP_H4_p12_1e4"] for x in items),default=""),
   "new_n_assays_H4_gt_H3_baseline":sum(x["H4_gt_H3_at_baseline"] for x in items) if qn else "",
   "new_n_assays_H4_gt_H3_all5":sum(x["H4_gt_H3_at_all_5_priors"] for x in items) if qn else "",
   "new_n_assays_H4_ge_0p8_highest_prior":sum(x["H4_ge_0p8_at_highest_p12"] for x in items) if qn else "",
   "new_median_fraction_snp_MAF_diff_gt_0p1":statistics.median(x["MAF_diff_gt_0p1_fraction"] for x in items) if qn else "",
   "new_total_palindromic_records_across_assays":sum(x["n_palindromic_recorded"] for x in items) if qn else "",
   "new_max_H4_shift_scalar_QTL_N":max((x["max_H4_shift_QTL_scalar_N"] for x in items),default=""),
   "science_gate":"EXPLORATORY_POSTERIOR_PRIOR_AF_QC_ONLY" if qn else "MOLECULAR_QTL_NOT_TESTED",
   "causal_status":"NOT_ESTABLISHED",
  })
 if len(merged)!=2225 or sum(x["new_DIRECT_prior_grid_assays"] for x in merged)!=646:
  raise ValueError("Full gene candidate universe or archived QTL coverage changed")
 return merged,assay_out

def main():
 p=argparse.ArgumentParser()
 for name in ("genes","gene-regions","prior-grid","af-per-assay","qtl-n","verified-replay","out-dir"):
  p.add_argument("--"+name,required=True,type=Path)
 a=p.parse_args()
 src={
  "genes":a.genes,"links":a.gene_regions,"grid":a.prior_grid,
  "af":a.af_per_assay,"nsens":a.qtl_n,"replay":a.verified_replay,
 }
 data={k:read(v) for k,v in src.items()}
 genes,assays=integrate(data["genes"],data["links"],data["grid"],data["af"],data["nsens"],data["replay"])
 if a.out_dir.exists() and any(a.out_dir.iterdir()):
  raise FileExistsError("No overwrite of derived research outputs")
 a.out_dir.mkdir(parents=True,exist_ok=True)
 writer(a.out_dir/"IS_2225_GENE_ABF_PRIOR_AF_INTEGRATED.tsv",genes)
 writer(a.out_dir/"IS_646_GENE_TISSUE_PRIOR_AF_CONTEXT.tsv",assays)
 summary={
   "status":"FULL_CANDIDATE_UNIVERSE_RETAINED_EXPLORATORY_PRIOR_AF_CONTEXT",
   "gene_ids":len(genes),"gene_region_links":2425,
   "qtl_genes_with_direct_p12_grid":sum(int(g["new_DIRECT_prior_grid_assays"])>0 for g in genes),
   "qtl_gene_tissue_assays":len(assays),
   "gene_with_high_prior_H4_0p8":sum(int(g["new_n_assays_H4_ge_0p8_highest_prior"] or 0)>0 for g in genes),
   "no_candidate_exclusions":True,
   "no_causal_genes_established":True,
   "sensitivity_p12_prior_grid":[1e-6,3e-6,1e-5,3e-5,1e-4],
   "priority_is_exploratory_not_causal":True,
   "matched_QTL_LD_not_verified":True,
   "sources_sha256":{k:sha(v) for k,v in src.items()},
 }
 (a.out_dir/"IS_2225_PRIOR_AF_INTEGRATION_MANIFEST.json").write_text(json.dumps(summary,indent=2)+"\n")
 print("IS_2225_PRIOR_AF_INTEGRATION_PASS",json.dumps({k:v for k,v in summary.items() if k in ("gene_ids","qtl_genes_with_direct_p12_grid","qtl_gene_tissue_assays","gene_with_high_prior_H4_0p8")}))
if __name__=="__main__":
 main()
