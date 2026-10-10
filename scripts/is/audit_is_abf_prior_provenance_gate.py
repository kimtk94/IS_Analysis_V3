#!/usr/bin/env python3
"""Fail-closed provenance gate for legacy IS coloc p12 interpretation."""
import argparse,csv,json,hashlib,re
from pathlib import Path
def read(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def run(master,joint,replay_script,out):
 original=read(master);scenarios=read(joint)
 if len(original)!=646 or len(scenarios)!=30:raise ValueError("unexpected source count")
 src=replay_script.read_text()
 if not re.search(r"p1\s*=\s*1e-4\s*,\s*p2\s*=\s*1e-4\s*,\s*p12\s*=\s*1e-5",src):
  raise ValueError("replay parameter signature changed")
 statuses={x["status"] for x in original}
 if statuses!={"PASS"}:raise ValueError("archived model has unexpected status")
 flagged=[]
 for s in scenarios:
  a=float(s["archived_H4"]);lo=float(s["prior_H4_min"]);hi=float(s["prior_H4_max"])
  if not lo<=a<=hi:raise ValueError("prior grid excludes archive")
  flagged.append({"gene_id":s["gene_id"],"locus":s["locus"],"dataset_key":s["dataset_key"],
   "archived_pp_h4":a,"p12_grid_min_h4":lo,"p12_grid_max_h4":hi,
   "p12_grid_crosses_h4_0_8":lo<.8<=hi,
   "prior_origin":"REPLAY_SCRIPT_CONFIRMED; HISTORICAL_ORIGINAL_NOT_INDEPENDENTLY_DOCUMENTED",
   "report_status":"CONDITIONAL_ONLY_NO_VALIDATED_CAUSAL_GENE"})
 result={"assays":len(flagged),"archived_assays":len(original),
  "replay_script_p1":1e-4,"replay_script_p2":1e-4,"replay_script_p12":1e-5,
  "historical_original_p12_verified_from_original_execution_provenance":False,
  "numerical_baseline_replay_confirms_algorithm_and_assumED_parameter_set_only":True,
  "p12_grid_cross_h4_0_8":sum(x["p12_grid_crosses_h4_0_8"] for x in flagged),
  "claim_gate":"CONDITIONAL_PRIOR_ASSUMPTION_NOT_INDEPENDENT_SOURCE_VERIFIED",
  "original_master_sha256":hashlib.sha256(master.read_bytes()).hexdigest(),
  "joint_sensitivity_sha256":hashlib.sha256(joint.read_bytes()).hexdigest(),
  "replay_script_sha256":hashlib.sha256(replay_script.read_bytes()).hexdigest()}
 out.mkdir(parents=True,exist_ok=True)
 with (out/"IS_30_P12_PRIOR_PROVENANCE_GATE.tsv").open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(flagged[0]),delimiter="\t");w.writeheader();w.writerows(flagged)
 (out/"IS_30_P12_PRIOR_PROVENANCE_GATE.json").write_text(json.dumps(result,indent=2)+"\n")
 return result
if __name__=="__main__":
 p=argparse.ArgumentParser()
 for x in ("master","joint","replay-script","out"):p.add_argument("--"+x,type=Path,required=True)
 a=p.parse_args();print(json.dumps(run(a.master,a.joint,a.replay_script,a.out),indent=2))
