# IS 2,225-gene whole-universe workflow readiness audit — 2026-10-11

The existing reference tables were reconciled without removing any gene: **2,225 distinct stable Ensembl gene IDs, 2,425 gene–region combinations**, 869 East Asian/Japanese pairs and 1,556 European pairs. No gene has causality established; all current sQTL/pQTL gene-level fields remain NOT_ASSESSED and must not be interpreted as negative biology.

The newly computed operational tiers are:
- A: 14 genes with source-linked SNP-level QTL replay context ready.
- B: 29 genes with historical ABF assays and SNP-level replay source pending.
- C: 826 East-Asian positional GWAS candidates without attached QTL replay evidence.
- D: 1,356 European-only positional GWAS candidates, which should not borrow EAS LD as matched population LD.

No causal priority claim is justified by the purely operational tier/score. Preserve original 2,225 discovery denominator rather than truncate to the 14 or historical Top 8. The historical 646 ABF assay records had 30 PASS and 616 MISSING_INPUT; 30 are assay records, NOT 30 unique genes.

We generated the complete 2,225-row TSV and **40 candidate existing-data follow-ups outside G0022**. The first queue includes SH3PXD2A, COL4A2, CALHM2, COL4A1, NEURL, FGF5, INA, C4orf22, SLK, COL17A1, PRDM8, ANTXR2 and CALHM1. These names are a workflow sort based on available earlier data, not genetic probability, novel effect replication, or validated drug targets.

Inputs, read-only:
- `/srv/is-analysis/results/is/audits/is_existing_2225_evidence_matrix_20261011_v1/IS_2225_GENE_EVIDENCE_MATRIX.tsv`
- `/srv/is-analysis/results/is/audits/is_existing_2225_evidence_matrix_20261011_v1/IS_2425_GENE_REGION_EVIDENCE.tsv`
- original manifest JSON with source SHA256 and one legacy NEURL1 stable gene alias.

Outputs:
- `/srv/is-analysis/results/is/audits/is_2225_workflow_priority_20261011_v1/IS_2225_NONCAUSAL_WORKFLOW_PRIORITY.tsv`
- `/srv/is-analysis/results/is/audits/is_2225_workflow_priority_20261011_v1/IS_TOP40_NON_G0022_EXISTING_DATA_QUEUE.tsv`
- `/srv/is-analysis/results/is/audits/is_2225_workflow_priority_20261011_v1/IS_2225_WORKFLOW_PRIORITY_MANIFEST.json`

Executable generator: `scripts/is/audit_existing_2225_evidence_priority_queue.py` and fail-closed regression tests `tests/test_is_2225_workflow_priority_queue.py`.

This uses no extra controlled database, no new GWAS sample data, no changes to existing canonical results and no hard-gene filtering. Suggested next research unit: recover SNP-source inputs for the 29 previously computed ABF genes and inspect 14 fully replay-ready genes' LD, priors and coloc reproducibility first; contrast by ancestry rather than combining EUR and EAS inference.
