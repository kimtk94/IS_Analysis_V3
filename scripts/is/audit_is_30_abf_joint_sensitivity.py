#!/usr/bin/env python3
"""Join source-grounded 30 ABF replays to previously computed N/prior sensitivities."""
import argparse,csv,json,hashlib
from collections import defaultdict,Counter
from pathlib import Path

def read(path):
 with path.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def key(r):return (r["locus"],r["dataset_key"],r.get("gene_base",r.get("gene_id")))
def run(risk,nscalar,prior,out):
 source=read(risk);ns=read(nscalar);ps=read(prior)
 if len(source)!=30 or len(ns)!=646 or len(ps)!=3230:raise ValueError("Unexpected source table cardinality")
 ni={key(x):x for x in ns}
 pp=defaultdict(dict)
 for x in ps:pp[key(x)][float(x["conditional_p12"])]=x
 grids={1e-6,3e-6,1e-5,3e-5,1e-4}
 outrows=[]
 for x in source:
  k=(x["locus"],x["dataset_key"],x["gene_id"])
  if k not in ni or set(pp[k])!=grids:raise ValueError("Missing original ABF prior/N sensitivity")
  n=ni[k];q=pp[k];base=float(x["prior_H4"])
  if abs(base-float(n["archived_H4"]))>1e-8 or abs(base-float(q[1e-5]["PP_H4"]))>1e-8:
   raise ValueError("Prior or scalar baseline mismatch")
  low=min(float(v["PP_H4"]) for v in q.values())
  high=max(float(v["PP_H4"]) for v in q.values())
  nlow=min(float(n[z]) for z in ("H4_first","H4_min","H4_median","H4_max"))
  nhigh=max(float(n[z]) for z in ("H4_first","H4_min","H4_median","H4_max"))
  outrows.append(dict(
   gene_id=x["gene_id"],locus=x["locus"],dataset_key=x["dataset_key"],source_nsnps=x["source_nsnps"],
   archived_H4=base,prior_p12_1e6_H4=float(q[1e-6]["PP_H4"]),
   prior_p12_1e5_H4=float(q[1e-5]["PP_H4"]),prior_p12_1e4_H4=float(q[1e-4]["PP_H4"]),
   prior_H4_min=low,prior_H4_max=high,
   scalar_N_H4_min=nlow,scalar_N_H4_max=nhigh,
   scalar_N_max_abs_H4_shift=float(n["max_abs_H4_shift_vs_first"]),
   source_AF_diff_gt_0p1_fraction=float(x["AF_abs_gt_0_1_fraction"]),
   source_QTL_N_unique=x["QTL_N_unique"],
   source_QTL_N_min=x["QTL_N_min"],source_QTL_N_max=x["QTL_N_max"],
   p12_variants_cross_H4_0p5=low<.5<=high,
   p12_variants_cross_H4_0p8=low<.8<=high,
   scalar_N_variants_cross_H4_0p5=nlow<.5<=nhigh,
   prior_H4_dynamic_range=high-low,
   independent_evidence_gate="SENSITIVITY_ONLY_NOT_NEW_CAUSAL_INFERENCE"))
 outrows.sort(key=lambda r:-r["archived_H4"])
 out.mkdir(parents=True,exist_ok=True)
 with (out/"IS_30_ABF_JOINT_SOURCE_N_AND_PRIOR_SENSITIVITY.tsv").open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(outrows[0]),delimiter="\t");w.writeheader();w.writerows(outrows)
 j=dict(assays=30,unique_gene_ids=len({r["gene_id"] for r in outrows}),
   p12_grid=sorted(grids),p12_cross_H4_0p5=sum(x["p12_variants_cross_H4_0p5"] for x in outrows),
   p12_cross_H4_0p8=sum(x["p12_variants_cross_H4_0p8"] for x in outrows),
   scalar_N_cross_H4_0p5=sum(x["scalar_N_variants_cross_H4_0p5"] for x in outrows),
   max_prior_H4_dynamic_range=max(x["prior_H4_dynamic_range"] for x in outrows),
   max_scalar_N_abs_H4_shift=max(x["scalar_N_max_abs_H4_shift"] for x in outrows),
   top5=[{k:r[k] for k in ("gene_id","dataset_key","archived_H4","prior_H4_min","prior_H4_max","scalar_N_H4_min","scalar_N_H4_max","source_AF_diff_gt_0p1_fraction")} for r in outrows[:5]],
   missing_inputs_unchanged=616,
   caveat="Prior sweep reweights archived summary under assumed original p12; scalar-N models do not model variant-specific N. Not independent biological replication or coloc validation.",
   source_sha256={v:hashlib.sha256(p.read_bytes()).hexdigest() for v,p in (("risk",risk),("nscalar",nscalar),("prior",prior))})
 (out/"IS_30_ABF_JOINT_SENSITIVITY_SUMMARY.json").write_text(json.dumps(j,indent=2)+"\n")
 return j
if __name__=="__main__":
 p=argparse.ArgumentParser()
 for name in ("risk","nscalar","prior","out"):p.add_argument("--"+name,type=Path,required=True)
 a=p.parse_args();print(json.dumps(run(a.risk,a.nscalar,a.prior,a.out),indent=2))
