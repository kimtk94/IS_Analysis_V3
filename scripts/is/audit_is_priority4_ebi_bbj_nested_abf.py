#!/usr/bin/env python3
"""Paired exact-source SNP-universe control, no new causal claim."""
import argparse,csv,json,hashlib,math
from pathlib import Path
from collections import defaultdict
GENES=("FGF5","C4orf22","CALHM2","NEURL")
def read(path):
 with path.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def write(path,data):
 with path.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter="\t");w.writeheader();w.writerows(data)
def run(full,nested,out):
 a=read(full);b=read(nested)
 if len(a)!=60 or len(b)!=60:raise ValueError("both complete 60-scenario grids required")
 def index(rows):
  ret={}
  for r in rows:
   key=(r["gene_symbol"],int(r["QTL_N_assumption"]),float(r["assumed_p12"]))
   if key in ret or r["status"]!="EXPLORATORY_ABF_SUCCESS":
    raise ValueError("duplicate/failed scenario")
   if not all(0<=float(r["PP_H"+str(i)])<=1 for i in range(5)):
    raise ValueError("posterior out of bounds")
   ret[key]=r
  return ret
 aa=index(a);bb=index(b)
 if aa.keys()!=bb.keys():raise ValueError("paired sensitivity cases do not align")
 result=[]
 for key in sorted(aa):
  gene,n,prior=key;x=aa[key];y=bb[key]
  if gene not in GENES:raise ValueError("unexpected gene")
  if int(x["snp_intersection"])<int(y["snp_intersection"]):raise ValueError("nested > full")
  if not math.isclose(float(x["assumed_p1"]),float(y["assumed_p1"])) or not math.isclose(float(x["assumed_p2"]),float(y["assumed_p2"])):
   raise ValueError("prior discrepancy")
  fx=float(x["PP_H4"]);ny=float(y["PP_H4"])
  result.append(dict(gene_symbol=gene,QTL_N_assumption=n,assumed_p12=prior,
   n_full=int(x["snp_intersection"]),n_nested=int(y["snp_intersection"]),
   full_ebi_bbj_PP_H4=fx,nested_ebi_bbj_PP_H4=ny,delta_full_minus_nested_PP_H4=fx-ny,
   full_ebi_bbj_PP_H3=float(x["PP_H3"]),nested_ebi_bbj_PP_H3=float(y["PP_H3"]),
   source_stats_matched_within_variant=True,
   abf_method_and_priors_fixed=True,
   causal_gene_validated=False))
 out.mkdir(parents=True,exist_ok=True)
 write(out/"IS_PRIORITY4_SAME_SOURCE_FULL_VS_NESTED_60_SCENARIOS.tsv",result)
 summary=[]
 for g in GENES:
  rows=[r for r in result if r["gene_symbol"]==g]
  if len(rows)!=15:raise ValueError("missing paired gene scenarios")
  base=next(r for r in rows if r["QTL_N_assumption"]==175 and r["assumed_p12"]==1e-5)
  d=[r["delta_full_minus_nested_PP_H4"] for r in rows]
  summary.append(dict(gene_symbol=g,n_nested=base["n_nested"],n_full=base["n_full"],
    baseline_nested_PP_H4=base["nested_ebi_bbj_PP_H4"],
    baseline_full_PP_H4=base["full_ebi_bbj_PP_H4"],
    baseline_PP_H4_full_minus_nested=base["delta_full_minus_nested_PP_H4"],
    delta_PP_H4_across_15_min=min(d),delta_PP_H4_across_15_max=max(d),
    scenarios_full_H4_lower_than_nested=sum(v<0 for v in d),
    scenario_count=15,allele_matched_same_QTL_BBJ_beta_se_sources=True,
    original_v8_data_used_in_either_model=False))
 write(out/"IS_PRIORITY4_SAME_SOURCE_FULL_VS_NESTED_SUMMARY.tsv",summary)
 manifest=dict(status="SAME_SOURCE_NESTED_SNP_UNIVERSE_CONTROL_ABF_ONLY",
    paired_scenarios=60,genes=4,study_sources_fixed_between_paired_runs=True,
    identical_ancestry_QTL_effective_N_assumptions=True,
    p12_p1_p2_priors_identical_paired=True,
    full_minus_nested_H4_all_negative=all(r["delta_full_minus_nested_PP_H4"]<0 for r in result),
    cannot_attribute_historical_vs_new_H4_to_variant_coverage_alone=True,
    original_historical_GTEx_v8_model_not_authenticated=True,
    source_GTEx_QTL_donor_LD_not_verified=True,
    no_SuSiE_finemapping_or_causal_gene_promotion=True,
    summary=summary,
    source_SHA256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (full,nested)})
 (out/"IS_PRIORITY4_NESTED_SNP_UNIVERSE_CONTROL_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")
 return manifest
if __name__=="__main__":
 p=argparse.ArgumentParser()
 for x in ("full","nested","out"):p.add_argument("--"+x,type=Path,required=True)
 a=p.parse_args();r=run(a.full,a.nested,a.out)
 for v in r["summary"]:print("GENE",v)
