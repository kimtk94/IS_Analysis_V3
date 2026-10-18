#!/usr/bin/env python3
"""P1.4 source-level metadata sufficiency audit; no inferred sample overlap."""
import csv,json,hashlib,sys
from pathlib import Path
P1=Path("/srv/is-analysis/results/is/audits/is_p1_full_universe_gh_39124908298_attempt_1")
P0=Path("/srv/is-analysis/results/is/audits/is_priority4_evidence_scorecard_gh_39116626473_attempt_1")
def read(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def main():
 provenance=json.loads((P1/"P1_GITHUB_RUN_PROVENANCE.json").read_text())
 for name,expected in provenance["output_hashes"].items():
  f=P1/name
  if not f.is_file() or hashlib.sha256(f.read_bytes()).hexdigest()!=expected:raise ValueError("P1_SOURCE_HASH_DRIFT")
 universe=read(P1/"IS_2225_NONCAUSAL_WORKFLOW_PRIORITY.tsv")
 if len(universe)!=2225 or len({x["gene_id_stable"] for x in universe})!=2225:raise ValueError("P1_UNIVERSE_DRIFT")
 score=read(P0/"IS_PRIORITY4_INTEGRATED_EVIDENCE_SCORECARD.tsv")
 if len(score)!=4 or {x["gene"] for x in score}!={"CALHM2","FGF5","NEURL","C4orf22"}:raise ValueError("P0_SCOPE_DRIFT")
 rows=[]
 for x in score:
  rows.append({"gene":x["gene"],"locus":x["locus"],
    "allele_exact_descriptive_snp_count":int(x["allele_exact_GWAS_SNPs"]),
    "BBJ_GIGASTROKE_shared_gws_snp_count":int(x["same_SNP_BBJ_and_GIGA_GWS"]),
    "independent_cohort_overlap":"UNVERIFIED_NOT_IN_AUDITED_SCORECARD",
    "BBJ_study_case_control_N":"UNVERIFIED_NOT_IN_AUDITED_SCORECARD",
    "GIGASTROKE_EAS_study_N":"UNVERIFIED_NOT_IN_AUDITED_SCORECARD",
    "QTL_study_N":"UNVERIFIED_NOT_IN_AUDITED_SCORECARD",
    "QTL_study_matched_LD":"NOT_AVAILABLE_IN_AUDITED_SCORECARD",
    "formal_multisignal_coloc":x["formal_multisignal_coloc_status"],
    "independent_replication":"NOT_ESTABLISHED",
    "causal_gene":"NOT_ESTABLISHED"})
 result={"status":"P1_4_AUDITED_METADATA_GAPS_NOT_REPLICATION",
   "p1_gene_count":len(universe),"p0_gene_count":len(score),
   "source_p1_run":39124908298,"source_p0_run":39116626473,
   "unverified_cohort_overlap":len(rows),"qtl_N_not_verified":len(rows),
   "matched_qtl_LD_not_verified":len(rows),
   "required_followup":"RETRIEVE_OFFICIAL_BBJ_GIGASTROKE_QTL_STUDY_DESIGN_N_AND_OVERLAP_DOCUMENTS",
   "interpretation":"Absence from audited scorecard is NOT proof that information does not exist elsewhere.",
   "genes":rows}
 print("P14_METADATA_JSON="+json.dumps(result,sort_keys=True),flush=True)
 print("P14_METADATA_GATE=PASS_WITH_UNRESOLVED_EVIDENCE",flush=True)
if __name__=="__main__":
 try:main()
 except Exception as exc:
  print("P14_AUDIT_ERROR="+repr(exc),file=sys.stderr,flush=True);sys.exit(1)
