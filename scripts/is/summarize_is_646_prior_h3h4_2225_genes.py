#!/usr/bin/env python3
"""Summarize original 646 coloc.abf assays over five p12 priors, no recompute.

This is a sensitivity diagnostic, NOT a validation of causal variants,
expression mediation, matched-LD SuSiE or stroke disease-cell specificity.
Do not change the original candidate universe of 2,225 positional genes.
"""
import argparse,csv,hashlib,json,math
from collections import defaultdict
from pathlib import Path

P12=(1e-6,3e-6,1e-5,3e-5,1e-4)
LOCI=("BBJ_IS_L001","BBJ_IS_L002","BBJ_IS_L003","BBJ_IS_L004")
def read(path):
 with path.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def save(path,rows):
 with path.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
  w.writeheader();w.writerows(rows)
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for block in iter(lambda:f.read(4194304),b""):h.update(block)
 return h.hexdigest()
def audit(rows,genes):
 if len(rows)!=3230 or len(genes)!=2225 or len({r["gene_id_stable"] for r in genes})!=2225:
  raise ValueError("Requires full 646 x five grid and unchanged 2225 genes")
 groups=defaultdict(dict)
 for r in rows:
  key=(r["locus"],r["dataset_key"],r["gene_base"])
  p=float(r["p12"])
  if p not in P12 or r["locus"] not in LOCI or p in groups[key]:
   raise ValueError("Duplicated/invalid prior or locus")
  if abs(float(r["p1"])-.0001)>1e-13 or abs(float(r["p2"])-.0001)>1e-13:
   raise ValueError("Unexpected p1 p2")
  if any(not math.isfinite(float(r[f])) or not 0<=float(r[f])<=1 for f in ("PP_H0","PP_H1","PP_H2","PP_H3","PP_H4")):
   raise ValueError("Posterior numeric QC fail")
  groups[key][p]=r
 if len(groups)!=646 or any(set(r)!=set(P12) for r in groups.values()):
  raise ValueError("Incomplete 646 assay five-prior grid")
 ass=[]
 for k,prior in groups.items():
  s={p:prior[p] for p in P12}
  if len({int(x["n_snps"]) for x in s.values()})!=1:
   raise ValueError("Input SNP count changed across prior values")
  at=lambda p,f:float(s[p][f])
  ass.append({
   "locus":k[0],"dataset_key":k[1],"gene_id":k[2],
   "n_original_SNPs":int(s[1e-5]["n_snps"]),
   "original_p12_1e5_H3":at(1e-5,"PP_H3"),
   "original_p12_1e5_H4":at(1e-5,"PP_H4"),
   "original_H4_gt_H3":at(1e-5,"PP_H4")>at(1e-5,"PP_H3"),
   "high_p12_1e4_H3":at(1e-4,"PP_H3"),
   "high_p12_1e4_H4":at(1e-4,"PP_H4"),
   "high_H4_gt_H3":at(1e-4,"PP_H4")>at(1e-4,"PP_H3"),
   "high_H4_ge_0p8":at(1e-4,"PP_H4")>=.8,
   "high_only_H4_ge_0p8":at(1e-4,"PP_H4")>=.8 and at(1e-5,"PP_H4")<.8,
   "prior_sensitivity_H4_delta":at(1e-4,"PP_H4")-at(1e-5,"PP_H4"),
   "prior_sensitivity_gate":"EXPLORATORY_PRIOR_DEPENDENT",
   "disease_gtex_causal_gene_status":"NOT_ESTABLISHED"
  })
 ass.sort(key=lambda r:(-r["high_p12_1e4_H4"],r["locus"],r["gene_id"],r["dataset_key"]))
 allp=[];locp=[]
 for p in P12:
  m=[x[p] for x in groups.values()]
  h4gt=sum(float(r["PP_H4"])>float(r["PP_H3"]) for r in m)
  above=[r for r in m if float(r["PP_H4"])>=.8]
  allp.append({
   "p12":p,"n_original_gene_tissue_assays":len(m),
   "n_H4_gt_H3":h4gt,
   "n_H4_ge_0p5":sum(float(r["PP_H4"])>=.5 for r in m),
   "n_H4_ge_0p8":len(above),
   "distinct_genes_H4_ge_0p8":len({r["gene_base"] for r in above}),
   "distinct_loci_H4_ge_0p8":len({r["locus"] for r in above}),
   "n_H1_dominant_0p5":sum(float(r["PP_H1"])>.5 for r in m),
   "n_H3_dominant_0p5":sum(float(r["PP_H3"])>.5 for r in m),
   "max_H4":max(float(r["PP_H4"]) for r in m),
   "interpretation":"ABF_HYPOTHESIS_PRIOR_SENSITIVITY_ONLY"
  })
  for locus in LOCI:
   sub=[x for x in m if x["locus"]==locus]
   locp.append({
    "p12":p,"locus":locus,"n_assays":len(sub),
    "n_H4_gt_H3":sum(float(r["PP_H4"])>float(r["PP_H3"]) for r in sub),
    "n_H4_ge_0p8":sum(float(r["PP_H4"])>=.8 for r in sub),
    "n_H1_dominant":sum(float(r["PP_H1"])>.5 for r in sub),
    "n_H3_dominant":sum(float(r["PP_H3"])>.5 for r in sub)
   })
 if sum(r["original_H4_gt_H3"] for r in ass)!=40:
  raise ValueError("Original baseline H4>H3 count drifted")
 if sum(r["high_H4_gt_H3"] for r in ass)!=565:
  raise ValueError("Prior sensitivity high-H4>H3 count drifted")
 if sum(r["high_H4_ge_0p8"] for r in ass)!=12:
  raise ValueError("High-prior 12 H4>=.8 assay gate drifted")
 genes_by_id={x["gene_id_stable"]:x for x in genes}
 gene_assays=defaultdict(list)
 for a in ass:
  if a["gene_id"] not in genes_by_id:raise ValueError("Molecular trait outside 2225 positional universe")
  gene_assays[a["gene_id"]].append(a)
 allgenes=[]
 for g in genes:
  items=gene_assays[g["gene_id_stable"]]
  allgenes.append({
   "gene_id_stable":g["gene_id_stable"],
   "canonical_gene_symbol":g["canonical_GENCODE_v19_gene_symbol"],
   "n_positional_regions":g["n_positional_regions"],
   "n_archived_coloc_gene_tissue_assays":len(items),
   "baseline_max_H4":max((r["original_p12_1e5_H4"] for r in items),default=""),
   "high_prior_max_H4":max((r["high_p12_1e4_H4"] for r in items),default=""),
   "n_H4_ge_0p8_baseline":sum(r["original_p12_1e5_H4"]>=.8 for r in items),
   "n_H4_ge_0p8_highprior":sum(r["high_H4_ge_0p8"] for r in items),
   "human_reference_feature_status":g["human_reference_feature_status"],
   "causal_status":"NOT_ESTABLISHED",
   "interpretation":("NO_ARCHIVED_QTL_TEST_NOT_NEGATIVE" if not items else "PRIOR_SENSITIVITY_ONLY")
  })
 return allp,locp,ass,allgenes
def main():
 parser=argparse.ArgumentParser()
 for key in ("prior-grid","gene-universe","out-dir"):
  parser.add_argument("--"+key,type=Path,required=True)
 a=parser.parse_args()
 if a.out_dir.exists() and any(a.out_dir.iterdir()):raise FileExistsError("Refuse output overwrite")
 allp,locp,ass,genes=audit(read(a.prior_grid),read(a.gene_universe))
 a.out_dir.mkdir(parents=True,exist_ok=True)
 for filename,rows in [
  ("IS_646_PRIOR_SENSITIVITY_COUNTS.tsv",allp),
  ("IS_646_PRIOR_SENSITIVITY_FOUR_LOCI.tsv",locp),
  ("IS_646_PRIOR_ASSAY_RECLASSIFICATION.tsv",ass),
  ("IS_2225_GENE_WIDE_H3_H4_EVIDENCE_GATES.tsv",genes)
 ]:save(a.out_dir/filename,rows)
 summary={
  "status":"REPRODUCED_ARCHIVED_P12_GRID_DESCRIPTIVE_NO_POSTHOC_PROMOTION",
  "archived_assays":646,"prior_values":P12,"all_positional_genes_retained":len(genes),
  "genes_without_archived_QTL_tests":sum(r["n_archived_coloc_gene_tissue_assays"]==0 for r in genes),
  "baseline_H4_gt_H3":40,"highprior_H4_gt_H3":565,
  "baseline_H4_ge_0p8":0,"highprior_H4_ge_0p8":12,
  "highprior_distinct_genes_H4_ge_0p8":next(r for r in allp if r["p12"]==1e-4)["distinct_genes_H4_ge_0p8"],
  "no_shared_causal_variant_claims":True,
  "all_prior_analyses_are_same_646_assays_not_independent_replication":True,
  "original_prior_grid_sha256":sha(a.prior_grid),
  "original_2225_gene_universe_sha256":sha(a.gene_universe)
 }
 (a.out_dir/"IS_ABF_PRIOR_GATE_MANIFEST.json").write_text(json.dumps(summary,indent=2)+"\n")
 print("IS_PRIOR_646_AND_2225_GATES_PASS",json.dumps(summary),flush=True)
if __name__=="__main__":main()
