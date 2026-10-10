#!/usr/bin/env python3
"""Review historical ABF replay gene claims, separating numerical PASS from biological validity."""
import csv,json,argparse,hashlib
from pathlib import Path
from collections import Counter,defaultdict

def records(path):
 with path.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def run(queue,replay,out):
 candidates=[x for x in records(queue) if x["workflow_priority_tier"]=="A_EXISTING_SNP_REPLAY_READY"]
 rows=records(replay)
 if len(candidates)!=14 or len(rows)!=646:raise ValueError("Input counts changed")
 passed=[x for x in rows if x["status"]=="PASS"]
 if len(passed)!=30:raise ValueError("PASS source drift")
 indexed=defaultdict(list)
 for x in passed:indexed[x["gene_base"]].append(x)
 outrows=[]
 for g in candidates:
  xs=indexed.get(g["gene_id_stable"],[])
  if len(xs)!=int(g["archived_replay_PASS"]):raise ValueError("Gene-linked replay count disagree: "+g["symbol"])
  errors=[float(x["max_abs_hypothesis_delta"]) for x in xs]
  if any(not (0<=z<=1e-8) for z in errors):raise ValueError("ABF numerical mismatch")
  outrows.append(dict(gene_id_stable=g["gene_id_stable"],symbol=g["symbol"],region_ids=g["region_ids"],
   numerical_ABF_replay_PASS=len(xs),source_ABF_assays=int(g["archived_ABF_assays"]),
   missing_original_SNP_inputs=int(g["archived_replay_MISSING_INPUT"]),
   max_replay_nsnps=max([int(x["actual_nsnps"]) for x in xs],default=0),
   max_replayed_PP_H4=max([float(x["PP_H4"]) for x in xs],default=None),
   min_replayed_PP_H4=min([float(x["PP_H4"]) for x in xs],default=None),
   max_abs_ABF_delta=max(errors,default=None),
   total_assay_maf_difference_gt_0_1=sum(int(x["maf_diff_gt_0_1"]) for x in xs),
   any_assay_variable_qtl_n=any(int(x["observed_qtl_n_range"])>0 for x in xs),
   replayed_tissues=";".join(sorted({x["dataset_key"] for x in xs})),
   claim_status="NUMERICAL_REPLAY_ONLY_NOT_COLOC_VALIDATED",
   next_step="VERIFY_QTL_EFFECT_AF_N_BY_SNP_AND_HYPOTHESIS_PRIORS"))
 extra=[x for x in passed if x["gene_base"] not in {g["gene_id_stable"] for g in candidates}]
 out.mkdir(parents=True,exist_ok=True)
 with (out/"IS_14_REPLAY_READY_GENE_EVIDENCE_QC.tsv").open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(outrows[0]),delimiter="\t");w.writeheader();w.writerows(outrows)
 summary=dict(gene_universe_reviewed=14,archived_assays=646,assays_numerically_replayed=30,
               linked_pass_in_14=sum(x["numerical_ABF_replay_PASS"] for x in outrows),
               unlinked_pass_assays=len(extra),
               missing_assays=616,source_replay_pass_details_not_causal=True,
               high_H4_is_not_replicated_causality=True,
               genes_with_replay_maf_discordance=sum(x["total_assay_maf_difference_gt_0_1"]>0 for x in outrows),
               genes_with_variable_QTL_sample_N=sum(x["any_assay_variable_qtl_n"] for x in outrows),
               original_queue_sha256=hashlib.sha256(queue.read_bytes()).hexdigest(),
               original_replay_sha256=hashlib.sha256(replay.read_bytes()).hexdigest(),
               no_causal_gene_confirmed=True)
 (out/"IS_14_REPLAY_EVIDENCE_QC_SUMMARY.json").write_text(json.dumps(summary,indent=2)+"\n")
 return summary
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--queue",type=Path,required=True);p.add_argument("--replay",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
 print(json.dumps(run(a.queue,a.replay,a.out),indent=2))
