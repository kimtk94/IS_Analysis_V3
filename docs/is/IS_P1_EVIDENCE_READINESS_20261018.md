# IS P1 evidence-readiness handoff — 2026-10-18

## Immutable baseline
- Discovery: 2,225 distinct stable ENSG gene IDs, 2,425 gene-region pairs. The four-gene P0 group is not the discovery universe.
- P0 GitHub Actions run: https://github.com/kimtk94/IS_Analysis_V3/actions/runs/39116626473 (P0_AUDIT_STATUS=COMPLETE).
- P0 research source commit: 8fb747d4bdf68b0495c2b77391c25cf89f11c9f5.
- P0 scorecard server-local path: /srv/is-analysis/results/is/audits/is_priority4_evidence_scorecard_gh_39116626473_attempt_1/IS_PRIORITY4_INTEGRATED_EVIDENCE_SCORECARD.tsv.
- P0 nominal QTL SNP counts (p<1e-4): CALHM2=3, FGF5=5, NEURL=0, C4orf22=0. None is an established causal gene.
- Formal multisignal colocalization status: BLOCKED_QTL_STUDY_LD_AND_N. Do not promote exploratory evidence to a causal claim.

## P1 sequential gates
| Gate | Output | Acceptance criteria |
| --- | --- | --- |
| P1.0 | Provenance freeze | input paths, SHA256, source commit, allele/build/ancestry metadata |
| P1.1 | 2,225-gene evidence-readiness audit | exactly 2,225 unique stable ENSG IDs and 2,425 gene-region pairs retained |
| P1.2 | Non-G0022 top-40 workflow queue | derived from full source universe, not a hard whitelist |
| P1.3 | Independent GWAS replication ledger | cohort/ancestry/effect allele/direction/N/overlap explicitly verified; UNKNOWN preserved |
| P1.4 | Molecular QTL coverage matrix | tissue/assay/N/cis-window/LD source recorded; missing LD/N blocks coloc |
| P1.5 | Figure/report publication gate | distinguish nominal QTL, reference LD, and formal coloc/MR |

## Existing implementation to reuse
- scripts/is/audit_existing_2225_evidence_priority_queue.py (options: --root DIRECTORY --out NEW_DIRECTORY).
- Expected files in --root: IS_2225_GENE_EVIDENCE_MATRIX.tsv, IS_2425_GENE_REGION_EVIDENCE.tsv, IS_GENE_EVIDENCE_MATRIX_MANIFEST.json.
- Outputs: IS_2225_NONCAUSAL_WORKFLOW_PRIORITY.tsv, IS_TOP40_NON_G0022_EXISTING_DATA_QUEUE.tsv, IS_2225_WORKFLOW_PRIORITY_MANIFEST.json.
- scripts/is/build_is_priority4_evidence_scorecard.py is descriptive and limited to four P0 genes; never use it to replace the full-universe analysis.

## Read-only server handoff (not yet executed)
1. Locate the directory containing all three source matrix files on the server; do not guess or download substitutes.
2. Verify manifest counts and source SHA256 before any research run.
3. Execute existing full-universe audit only with a fresh, uniquely named directory under /srv/is-analysis/results/is/audits/; never overwrite existing output.
4. Save start/end time, input hashes, checked-out git SHA, exact CLI and exit status in a run record.
5. Verify all 2,225 unique IDs and 2,425 region links before viewing or publishing the top-40 queue.
6. Compare the P0 four-gene scorecard with newly checked replication metadata, marking unverified replication as PENDING.
7. Do not trigger LIVE or mutation workflows, and do not infer causal genes, independent signals or replication from high reference-panel LD.

## Interpretation
- CALHM2 three nominal QTL variants form one high-LD cluster in EAS/JPT public reference genotypes. This is not a fine-mapped credible set.
- CALHM2 tags BBJ lead but poorly tags GIGASTROKE lead; this discrepancy alone is not proof of distinct causal variants.
- FGF5 shows one common BBJ/GIGASTROKE genome-wide-significant SNP in P0, but this does not prove gene mediation.
- Missing QTL coverage/N/LD, or absent nominal hits, is not evidence of a true biological negative.

## Automation acceptance
Only enable automatic P1 refresh after finding the authoritative 2,225-gene input root, freezing its checksums and observing a read-only non-overwriting P1 run with a provenance manifest. No assumption of completed P1 validation.
