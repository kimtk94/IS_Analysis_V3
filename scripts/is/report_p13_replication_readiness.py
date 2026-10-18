#!/usr/bin/env python3
"""Read-only P1.3 replication readiness, not evidence of independent replication."""
import csv,json,sys,hashlib
from pathlib import Path
P1=Path("/srv/is-analysis/results/is/audits/is_p1_full_universe_gh_39124908298_attempt_1")
P0=Path("/srv/is-analysis/results/is/audits/is_priority4_evidence_scorecard_gh_39116626473_attempt_1")
def rows(path):
 with path.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def main():
 manifest=json.loads((P1/"P1_GITHUB_RUN_PROVENANCE.json").read_text())
 full=P1/"IS_2225_NONCAUSAL_WORKFLOW_PRIORITY.tsv"
 top=P1/"IS_TOP40_NON_G0022_EXISTING_DATA_QUEUE.tsv"
 for f,h in manifest["output_hashes"].items():
  p=P1/f
  if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=h:raise RuntimeError("P1_HASH_DRIFT "+f)
 universe=rows(full);top40=rows(top)
 if len(universe)!=2225 or len({x["gene_id_stable"] for x in universe})!=2225:raise RuntimeError("UNIVERSE_INVALID")
 if len(top40)!=40:raise RuntimeError("TOP40_INVALID")
 score=rows(P0/"IS_PRIORITY4_INTEGRATED_EVIDENCE_SCORECARD.tsv")
 if len(score)!=4:raise RuntimeError("P0_SCORECARD_INVALID")
 score_by={x["gene"]:x for x in score}
 detailed=[]
 for r in top40:
  gene=r["symbol"];p0=score_by.get(gene)
  detailed.append(dict(gene=gene,gene_id_stable=r["gene_id_stable"],tier=r["workflow_priority_tier"],
   bbj_giga_same_SNP_GWS=int(p0["same_SNP_BBJ_and_GIGA_GWS"]) if p0 else None,
   qtl_nominal_p_lt_1e4=int(p0["QTL_nominal_p_lt_1e4"]) if p0 else None,
   allele_harmonization_status="P0_DESCRIPTIVE_ONLY" if p0 else "NOT_CHECKED_FOR_P1",
   independent_cohort_status="PENDING_COHORT_OVERLAP_ATTESTATION",
   ancestry_matched_effect_status="PENDING_STUDY_LEVEL_METADATA",
   sample_N_status="PENDING_STUDY_LEVEL_METADATA",
   qtl_ld_and_sample_N_status="PENDING_QTL_STUDY_LD_AND_N",
   formal_replication="NOT_ESTABLISHED",
   causal_gene_status="NOT_ESTABLISHED"))
 summary=dict(status="P1_3_REPLICATION_READINESS_ONLY",genes_total=2225,top40=len(top40),
  top40_with_p0_descriptive=len([x for x in detailed if x["allele_harmonization_status"]=="P0_DESCRIPTIVE_ONLY"]),
  top40_independent_replication_established=0,
  cohort_overlap_attested=0,
  source_p1_run=39124908298,source_p0_run=39116626473,
  checks="P1_SHA256_AND_P0_CARDINALITY_PASS",
  highest_readiness=detailed[:15])
 print("P13_SUMMARY_JSON="+json.dumps(summary,sort_keys=True),flush=True)
 print("P13_REPLICATION_READINESS=PASS_NONCAUSAL",flush=True)
if __name__=="__main__":
 try:main()
 except Exception as exc:
  print("P13_AUDIT_ERROR="+repr(exc),file=sys.stderr,flush=True);sys.exit(1)
