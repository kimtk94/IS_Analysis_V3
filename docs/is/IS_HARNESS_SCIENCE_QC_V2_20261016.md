# IS Harness science QC v2 — 2026-10-16

## Updated execution contract

The server checks its input-file metadata at each five-minute polling bucket. A delayed timer can still detect the most recent five-minute bucket; it no longer has to run during its first minute.

Five-minute checks do not automatically create full audit Issues. When source file size/mtime and audit schema match the latest validated manifest, the last full audit is younger than one hour, and there is no pending job, the worker skips redundant scientific scanning and Issue creation. Changed inputs or an audit older than one hour trigger full SHA-256 hashing and evidence generation. File size and mtime are only an optimization, NOT proof of identical content.

GitHub Actions scheduled execution is retained as a best-effort hourly failover; it checks every five minutes but does not enqueue when a recently successful full audit exists. Neither timer guarantees exact wall-clock execution. Original scientific input files and analysis worktree are not edited.

## Generated science QC artifacts

| Server output (under results/is/github_actions_runs/evidence/) | Purpose | Limitation |
| --- | --- | --- |
| p1_ld_gate.json | Source cohort ancestry and matched molecular-QTL LD readiness by region | Does not compute LD |
| p2_region_qc.json | Existing EAS5 region-gene-QTL dataset groups, p-value missingness and numeric validity | Does not recompute FASTA/allele harmonization |
| p4_nominal_candidate_triage.json | Individual nominal-QTL available genes and their evidence from the frozen 2,225-gene table | Not an MR/coloc causal gene ranking |
| research_manifest.json | Compact provenance and P1/P2/P4 result counts sent to web through a sanitized GitHub Issue record | Run PASS is not biological validation |

## First real server acceptance

GitHub Issue #99 (https://github.com/kimtk94/IS_Analysis_V3/issues/99) completed with PASS after 21.228 seconds.

- P1: 393 EAS5 region-gene-QTL dataset readiness rows; source ancestry verified 0 and matched molecular-QTL LD verified 0.
- P2: 393 region-gene-QTL dataset groups; invalid nonmissing QTL p 0, BBJ p 0, EAS p 0.
- P4: 149 nominal-QTL-available genes in a separate individual-level triage output.

The 149 candidates have not been established as ischemic stroke causal genes. The next substantive scientific gate is sourcing and validating matched cohort ancestry and QTL LD provenance, then performing qualified allele/QTL checks and controlled coloc re-analysis.

Audit metadata and displayed conclusions must be kept separate: 5-minute input checks, full audit completion, new data publication, and validated scientific discoveries are different events.
