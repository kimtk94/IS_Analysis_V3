#!/usr/bin/env python3
"""Exhaustively QC 4x3x5 *exploratory* coloc.abf sensitivity and create digest.
No confirmation of original GTEx_v8 model, independence, effect-size N or donor LD.
"""
import argparse,csv,json,math,hashlib
from collections import defaultdict
from pathlib import Path
GENES=("FGF5","C4orf22","CALHM2","NEURL")
COUNTS={"FGF5":2773,"C4orf22":2772,"CALHM2":2766,"NEURL":2739}
PRIORS=(1e-6,3e-6,1e-5,3e-5,1e-4)
NS=(125,150,175)
def read(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def analyze(grid,qc,ready,out):
 rows=read(grid);q=read(qc);inp=read(ready)
 if len(rows)!=60 or len(q)!=4 or len(inp)!=11050:raise ValueError("4x3x5 scenario/full source counts drift")
 cells={}
 by_gene=defaultdict(list)
 for r in rows:
  gene=r["gene_symbol"]
  if gene not in GENES or r["status"]!="EXPLORATORY_ABF_SUCCESS" or r["original_GTEx_v8_model_reproduced"]!="FALSE":
   raise ValueError("Posterior status or provenance contract mismatch")
  n=int(r["QTL_N_assumption"]);prior=float(r["assumed_p12"])
  if n not in NS or prior not in PRIORS:raise ValueError("Unknown sensitivity scenario")
  k=(gene,n,prior)
  if k in cells:raise ValueError("Duplicated scenario")
  cells[k]=r
  post=[float(r["PP_H"+str(i)]) for i in range(5)]
  if not all(math.isfinite(x) and 0<=x<=1 for x in post) or abs(sum(post)-1)>1e-6:
   raise ValueError("Posterior not a probability distribution")
  if int(r["snp_intersection"])!=COUNTS[gene]:raise ValueError("SNP denominator mismatch")
  by_gene[gene].append((prior,n,post[4],post[3]))
 if len(cells)!=60 or set(by_gene)!=set(GENES):raise ValueError("Missing scenario")
 result=[]
 for gene in GENES:
  rr=by_gene[gene]
  base=cells[gene,175,1e-5]
  prior175=[(p,float(cells[gene,175,p]["PP_H4"])) for p in PRIORS]
  nbase=[(n,float(cells[gene,n,1e-5]["PP_H4"])) for n in NS]
  # Prior H4 increasing monotonically with more cross-trait prior for fixed other priors.
  for i in range(1,len(prior175)):
   if not prior175[i][1]>prior175[i-1][1]:raise ValueError("Unexpected p12 ABF monotonicity")
  result.append(dict(gene_symbol=gene,matched_nonpal_SNPs=COUNTS[gene],
    H4_reference_N175_p12_1e5=float(base["PP_H4"]),
    H3_reference_N175_p12_1e5=float(base["PP_H3"]),
    all_scenarios_H4_min=min(t[2] for t in rr),
    all_scenarios_H4_max=max(t[2] for t in rr),
    H4_at_N125_p12_1e5=nbase[0][1],
    H4_at_N150_p12_1e5=nbase[1][1],
    H4_at_N175_p12_1e5=nbase[2][1],
    H4_N175_p12_1e6=prior175[0][1],
    H4_N175_p12_1e4=prior175[-1][1],
    H4_crosses_0p8_among_scenarios=min(t[2] for t in rr)<.8<max(t[2] for t in rr),
    original_v8_pph4_retained_unchanged=True,
    fullcis_gtex_canonical_source_authenticated=False,
    QTL_expression_effective_N_attested=False,
    source_dataset="EBI_EQTL_CATALOGUE_QTS000015_QTD000161_REPROCESSED"))
 out.mkdir(parents=True,exist_ok=True)
 tab=out/"IS_PRIORITY4_EBI_BBJ_ABF_60_SCENARIO_QC_SUMMARY.tsv"
 with tab.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(result[0]),delimiter="\t");w.writeheader();w.writerows(result)
 j=dict(audited_grid=60,genes=4,reference_cis_is_qtl_source="EBI_QTD000161_REPROCESSING",
    QTL_N_scenarios=list(NS),assumed_p12_scenarios=list(PRIORS),
    expected_BBJ_N_assumption=174686,
    original_historical_coloc_H4_modified=False,
    qtl_effective_N_is_unverified=True,
    QTL_donor_signed_LD_unavailable=True,
    coloc_susie_performed=False,
    source_status="EXPLORATORY_SINGLE_SIGNAL_ABF_PRIOR_AND_N_SENSITIVITY_ONLY",
    summary=result,
    SHA256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (grid,qc,ready)})
 (out/"IS_PRIORITY4_EBI_BBJ_ABF_60_SCENARIO_AUDIT.json").write_text(json.dumps(j,indent=2)+"\n")
 return j
if __name__=="__main__":
 p=argparse.ArgumentParser()
 for a in ("grid","qc","ready","out"):p.add_argument("--"+a,type=Path,required=True)
 a=p.parse_args()
 x=analyze(a.grid,a.qc,a.ready,a.out)
 for v in x["summary"]:print("GENE",v)
