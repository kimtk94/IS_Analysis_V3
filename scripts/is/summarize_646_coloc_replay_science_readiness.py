#!/usr/bin/env python3
"""Descriptive QC of full historical 646 per-SNP coloc.abf reproduction.

Posterior agreement is technical reproducibility only. Never conflate with
study-matched QTL LD, robust colocalization, or a causal-gene finding.
"""
from __future__ import annotations
import argparse,csv,json,hashlib,math
from collections import Counter,defaultdict
from pathlib import Path

def read(path):
 with path.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))

def sha(path):
 h=hashlib.sha256()
 with path.open("rb") as f:
  for piece in iter(lambda:f.read(4*1024*1024),b""):h.update(piece)
 return h.hexdigest()

def numeric(x):
 try:
  n=float(x)
  return n if math.isfinite(n) else None
 except (TypeError,ValueError):return None

def assess(rows, master):
 key=lambda r:(r["locus"],r["dataset_key"],r["gene_base"])
 if len(rows)!=len(master) or len({key(r) for r in rows})!=len(rows):
  raise ValueError("Duplicate or missing assay in recomputed batch")
 raw={key(r):r for r in master}
 if set(raw)!={key(r) for r in rows}:raise ValueError("Archived and recomputed keys differ")
 decorated=[]
 for r in rows:
  m=raw[key(r)]
  mafd=numeric(r.get("maf_diff_fraction"))
  maf_n=numeric(r.get("maf_diff_gt_0p1"))
  n_var=numeric(m["nsnps"])
  n_qtl=numeric(r.get("n_qtl_sample_N_distinct"))
  delta=numeric(r.get("delta_max"))
  orig=numeric(m["PP.H4"])
  replot=numeric(r.get("new_PP_H4"))
  if r["status"]=="PASS":
   if None in (mafd,maf_n,n_var,n_qtl,delta,orig,replot):
    raise ValueError("Missing QC summary from PASS replay: "+str(key(r)))
   if not 0<=mafd<=1 or not 0<=maf_n<=n_var or n_qtl<1:
    raise ValueError("Invalid frequency/sample QC")
   if delta>1e-8 or abs(orig-replot)>1e-8:
    raise ValueError("PASS posterior differs from archived result")
  maf_gate=("NOT_RECOMPUTED" if mafd is None else
            "MAF_DIFF_GT_0P1_IN_OVER_HALF_SNPS" if mafd>0.5 else
            "MAF_DIFF_GT_0P1_IN_OVER_QUARTER_SNPS" if mafd>0.25 else
            "MAF_DIFF_GT_0P1_IN_SOME_SNPS" if mafd>0 else
            "NO_MAF_DIFF_GT_0P1")
  decorated.append({
    "locus":r["locus"],"dataset_key":r["dataset_key"],
    "gene_base":r["gene_base"],
    "new_replay_status":r["status"],
    "original_H4":m["PP.H4"],
    "recomputed_H4":r.get("new_PP_H4",""),
    "posterior_max_abs_delta":r.get("delta_max",""),
    "n_expected_variants":m["nsnps"],
    "n_maf_delta_gt_0p1":r.get("maf_diff_gt_0p1",""),
    "maf_discordance_fraction":r.get("maf_diff_fraction",""),
    "maf_QC_review":maf_gate,
    "n_distinct_QTL_sample_N":r.get("n_qtl_sample_N_distinct",""),
    "QTL_sample_N_variable_across_SNPs":n_qtl is not None and n_qtl>1,
    "source_file":r["source_file"],
    "scientific_status":"NO_MATCHED_QTL_LD_NOT_CAUSAL",
  })
 return decorated

def main():
 p=argparse.ArgumentParser()
 p.add_argument("--replay",type=Path,required=True)
 p.add_argument("--master",type=Path,required=True)
 p.add_argument("--out-dir",type=Path,required=True)
 a=p.parse_args()
 rows=read(a.replay);master=read(a.master)
 if len(rows)!=646 or len(master)!=646:
  raise ValueError("Expected complete 646-assay source and derivative")
 records=assess(rows,master)
 if a.out_dir.exists() and any(a.out_dir.iterdir()):
  raise FileExistsError("Do not overwrite previous scientific QC")
 a.out_dir.mkdir(parents=True,exist_ok=True)
 with (a.out_dir/"IS_646_COLOC_SNP_REPLAY_SCIENCE_QC.tsv").open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter="\t")
  w.writeheader();w.writerows(records)
 passed=[r for r in records if r["new_replay_status"]=="PASS"]
 delta=[float(r["posterior_max_abs_delta"]) for r in passed]
 most=[r for r in passed if float(r["original_H4"])>=0.8]
 summary={
  "status":"PER_SNP_ABF_REPLAY_TECHNICAL_REPRODUCIBILITY_ONLY",
  "original_assays":646,
  "statuses":dict(Counter(r["new_replay_status"] for r in records)),
  "total_unique_legacy_gene_ids":len({r["gene_base"] for r in records}),
  "H4_at_least_0_8_among_replayed":len(most),
  "max_abs_posterior_discrepancy":max(delta) if delta else None,
  "maf_QC_classes":dict(Counter(r["maf_QC_review"] for r in records)),
  "variable_per_SNP_QTL_sample_N_assays":sum(r["QTL_sample_N_variable_across_SNPs"] for r in records),
  "replayed_total":len(passed),
  "source_SHA256":{"replay":sha(a.replay),"archived_master":sha(a.master)},
  "N_from_original_first_SNP_QTL":"YES_AS_IN_HISTORIC_SOURCE",
  "GTEx_molecular_QTL_cohort_matched_LD_verified":False,
  "BBJ_vs_GIGASTROKE_independent_replication_claim":False,
  "gene_causality_established":False,
  "multiple_assay_correction_performed":False,
  "interpretation":"Exact H0-H4 reproduction is necessary but not sufficient for robust biological coloc evidence.",
 }
 (a.out_dir/"IS_646_COLOC_REPLAY_SCIENCE_READINESS.json").write_text(json.dumps(summary,indent=2)+"\n")
 print("IS_646_COLOC_REPLAY_SCIENCE_QC_COMPLETE")
 print("REPLAY",json.dumps(summary["statuses"]))
 print("MAF",json.dumps(summary["maf_QC_classes"]))
 print("MAX_DELTA",summary["max_abs_posterior_discrepancy"])
 print("H4_GE_0P8",summary["H4_at_least_0_8_among_replayed"])
if __name__=="__main__":main()
