# IS Harness Autonomous Research Stages — v3 (2026-10-16)

## Objective and evidence rules

The five-minute server watcher checks source file metadata and stores a local heartbeat. When audited inputs have changed or the prior full evidence pass is at least 60 minutes old, `research-refresh` performs detailed P1/P2/P3/P4 checks and regenerates the validated web summary. `PASS` means that a reproducible *procedure* completed, **not** that a new causal gene was established.

## Source and output locations

Research sources are READ-ONLY within `/srv/is-analysis/results/is/audits`. Existing source files are not edited. Science reports are written atomically under `/srv/is-analysis/results/is/github_actions_runs/evidence/science/`:

| File | Output | Completion claim |
| --- | --- | --- |
| `p1_source_catalog.json` | QTL source identifiers, region and gene-assay counts, missing provenance | No ancestry assignment is inferred from a dataset name |
| `p1_cohort_ld_review.json` | Per gene/dataset ancestry and LD source-cohort gate | Missing verified LD is **BLOCKED**, not negative biology |
| `p2_independent_variant_qc.json` | Row-level aggregate SNP-key parser, source exact-match flags, QTL/EAS AF consistency and optional FAI-backed GRCh37 REF check | Source booleans are separate from independent REF/ALT validation |
| `p3_coloc_gates.json` | Candidate-level signed-LD/sample-size/allele/FASTA/coloc-input gate | No coloc until approved source files, validated matrix and tested executor are available |
| `p4_candidate_cross_evidence.json` | Nominal-QTL genes joined to BBJ/EAS source-flag evidence and independently established SNP/FASTA evidence | Not an FDR/MR causal target ranking |
| `science_pipeline_summary.json` | Sanitized stage-level progress summary and blocked reasons | No claims beyond validated input scope |
| `../harness_heartbeat.json` | Latest watcher observation, slot, and skip/enqueue state | A watcher heartbeat is not a scientific result |

## Critical correction: BBJ_exact_SNP and EAS_AIS_exact_SNP

The original EAS5 source TSV fields are **boolean exact-match flags**, not variant identifier strings. Counting nonempty strings as true counted `False` values erroneously. Issue [#111](https://github.com/kimtk94/IS_Analysis_V3/issues/111) confirms 809,245 total rows, BBJ flag true 805,037, EAS-AIS flag true 427,259, and BOTH true **423,051**. This does NOT establish 423,051 independently harmonized alleles. The earlier displayed 809,245 *exact matches* was an overcount and must not be used in manuscripts.

Independent SNP-key checks compare literal genomic identifiers, where present; a source `True` flag can never substitute for a reference FASTA or original GWAS variant/allele table. For source-flag-only data, independent literal-ID overlap remains 0. The full-GWAS identifier tables and approved GRCh37 FASTA must be obtained for independent allele validation.

## Human-reviewed source inputs (never infer, never fake)

Templates are located in `docs/is/templates/`. Copy the corresponding template under `/srv/is-analysis/results/is/github_actions_runs/evidence/science_inputs/` only when authoritative provenance is available:

- `P1_APPROVED_QTL_LD.tsv`: dataset/source cohort, source ancestry, signed LD source cohort/ancestry, actual on-server LD matrix, full SHA-256, reviewer, date and evidence URL. A cohort-matched LD file is NOT implied by an EAS/JPT population reference. File hash matching alone is not sufficient to prove valid SNP row order.
- `P2_APPROVED_GRCH37_REFERENCE.json`: approved GRCh37 reference genome FASTA and exact `.fai` index with recorded source and reviewer. The stage reads it without modifying it. Until approved, all independent reference checks are BLOCKED.
- `P3_APPROVED_COLOC_INPUTS.tsv`: reviewed, checked independent QTL/GWAS input file with variant alleles, sample sizes, assembly and exact file hash. Even if this registry exists, P3 will remain `AWAITING_VALIDATED_COLOC_EXECUTOR` until the LD-aligned SuSiE/coloc runner and independent tests are approved.

Never put access tokens, controlled individual-level genotypes, unreviewed private patient data or arbitrary shell commands in GitHub Issues/Logs.

## P1–P5 stage states

| Research stage | Autonomous safe action | Scientific gate |
| --- | --- | --- |
| **P1** | Enumerate distinct QTL datasets; check linked cohort/ancestry and actual LD SHA-256 | Correct source cohort and real matrix/variant-order proof |
| **P2** | Independently parse mapped genomic variants; distinguish booleans from SNP IDs; compute AF review flags for source dual matches; verify GRCh37 REF if approved FASTA exists | Independent allele ID, FASTA, strand, build and EAF checks |
| **P3** | Evaluate per-assay gate; do not run coloc with unmatched LD or missing sample sizes | Verified inputs, signed LD, allele order and separately tested SuSiE/coloc R executor |
| **P4** | Join 149 nominal-QTL candidate genes to source-flag and independent SNP/Fasta evidence | FDR/causality never inferred from nominal p |
| **P5** | Preserve 629 unexecuted regional pilot cohort labels; do not label these negative results | Pilot research design and QTL source inputs |

## Verifications and next science tasks

- Test and installation are managed by `.github/workflows/is-harness-research-jobs-ci.yml` and `.github/workflows/is-harness-safe-worker-update.yml`.
- Baseline acceptance [#99](https://github.com/kimtk94/IS_Analysis_V3/issues/99) established 393 assay groups, 149 nominal-QTL candidates and zero matched molecular QTL LD.
- Staged acceptance [#109](https://github.com/kimtk94/IS_Analysis_V3/issues/109) established the P1/P2/P3/P4 report paths and explicit P3 blocks.
- Source-flag correction [#111](https://github.com/kimtk94/IS_Analysis_V3/issues/111) fixed the previous exact-match overcount.
- P1 blocking task: [#100](https://github.com/kimtk94/IS_Analysis_V3/issues/100).

**Next actual advancement:** identify QTL provenance for the discovered source dataset identifiers; secure authoritative per-cohort source/LD information; independently map official GRCh37 reference alleles; then open the strictly gated P3 multi-signal colocalization path. No further candidate causality should be claimed before those gates are satisfied.
