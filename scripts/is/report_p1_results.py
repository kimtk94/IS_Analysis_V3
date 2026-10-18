#!/usr/bin/env python3
"""Read-only P1 full universe audit summary (no causal claims)."""
import csv,json,hashlib,sys
from collections import Counter
from pathlib import Path
ROOT=Path("/srv/is-analysis/results/is/audits/is_p1_full_universe_gh_39124908298_attempt_1")
FILES=("IS_2225_NONCAUSAL_WORKFLOW_PRIORITY.tsv","IS_TOP40_NON_G0022_EXISTING_DATA_QUEUE.tsv","IS_2225_WORKFLOW_PRIORITY_MANIFEST.json","P1_GITHUB_RUN_PROVENANCE.json")
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
 if not all((ROOT/n).is_file() for n in FILES):raise RuntimeError("REQUIRED_P1_FILES_MISSING")
 manifest=json.loads((ROOT/FILES[2]).read_text())
 provenance=json.loads((ROOT/FILES[3]).read_text())
 if manifest["total_source_unique_genes"]!=2225 or manifest["total_source_gene_region_pairs"]!=2425:raise RuntimeError("MANIFEST_UNIVERSE_DRIFT")
 if provenance.get("status")!="P1_COMPLETE_NONCAUSAL":raise RuntimeError("RUN_PROVENANCE_NOT_COMPLETE")
 for filename,digest in provenance["output_hashes"].items():
  if sha(ROOT/filename)!=digest:raise RuntimeError("OUTPUT_HASH_DRIFT: "+filename)
 def rows(name):
  with (ROOT/name).open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
 allrows=rows(FILES[0]);top=rows(FILES[1])
 if len(allrows)!=2225 or len({r["gene_id_stable"] for r in allrows})!=2225:raise RuntimeError("GENE_UNIVERSE_DRIFT")
 if len(top)!=40 or len({r["gene_id_stable"] for r in top})!=40:raise RuntimeError("TOP40_DRIFT")
 if any("G0022" in r["region_ids"] for r in top):raise RuntimeError("G0022_IN_TOP40")
 counts=dict(Counter(r["workflow_priority_tier"] for r in allrows))
 if counts!=manifest["tier_counts"]:raise RuntimeError("TIER_COUNT_DRIFT")
 summary={"status":"P1.1_AUDITED_NONCAUSAL","source_run":39124908298,"source_path":str(ROOT),
 "genes":len(allrows),"gene_region_pairs":2425,"top40_count":len(top),
 "tier_counts":counts,"top40":[{"gene":r["symbol"],"id":r["gene_id_stable"],
 "tier":r["workflow_priority_tier"],"readiness_score":int(r["workflow_priority_score_NONCAUSAL"]),
 "next_action":r["recommended_next_action"]} for r in top],
 "source_hashes":provenance["inputs"],"causal_gene_claims":"NOT_ESTABLISHED"}
 print("P1_SUMMARY_JSON="+json.dumps(summary,ensure_ascii=False,sort_keys=True),flush=True)
 print("P1_SUMMARY_AUDIT=PASS",flush=True)
if __name__=="__main__":
 try:main()
 except Exception as e:
  print("P1_SUMMARY_ERROR="+repr(e),file=sys.stderr,flush=True);sys.exit(1)
