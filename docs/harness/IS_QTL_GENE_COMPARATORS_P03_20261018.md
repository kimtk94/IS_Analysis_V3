# IS Harness P0.3 — gene-matched molecular evidence controls

**Prepared:** 2026-10-18 KST. **Status:** `GITHUB_CODE_READY / SERVER_NOT_EXECUTED`.

## Research aim
The P0.2 analysis targets pairwise LD among three CALHM2 QTL SNPs; this separate P0.3 study compares molecular evidence across **CALHM2, NEURL, FGF5, C4orf22** using the exact same SNP/allele/p-value filters. A zero nominal QTL count for NEURL or C4orf22 is not proof that the gene is biologically irrelevant: the assay, tissue and p-value cutoff condition this observation.

## Locked provenance and gates
- Input: `/srv/is-analysis/results/is/audits/is_priority4_fullcis_bbj_giga_snp_exact_20261011_v1/IS_PRIORITY4_FULLCIS_EBI_BBJ_GIGASTROKE_SNP_EVIDENCE.tsv` (20,689 gene×SNP records; source SHA256 recorded on execution).
- Restrict to GWAS allele-verified, nonpalindromic SNPs common to BBJ IS and GIGASTROKE EAS AIS, preserving per-gene QTL evidence.
- Prespecified QTL nominal cutoff: `p<1e-4`. Prior observed per-gene counts: CALHM2 3, NEURL 0, FGF5 5, C4orf22 0. These serve as a fail-closed consistency gate.
- Count same-variant genome-wide-significant BBJ/GIGASTROKE records and triple QTL/GWAS marginal-p overlaps, plus **descriptive** genomic ALT effect-direction agreement.
- Important: FGF5 and C4orf22 share chr4; CALHM2 and NEURL share chr10. Nearby GWAS SNPs are correlated and must not be treated as independent replication. EBI molecular N is not independently authenticated. No multicausal colocalization inferred.

## New source code and target artifacts
- `scripts/is/audit_is_priority4_qtl_gene_comparators.py`
- `tests/test_is_priority4_qtl_gene_comparators.py`
- Target new output `/srv/is-analysis/results/is/audits/is_priority4_qtl_gene_comparison_20261018_v1/` with `IS_QTL_GENE_COMPARATOR_SUMMARY.tsv`, `IS_QTL_GENE_COMPARATOR_SIGNIFICANT_VARIANTS.tsv`, `IS_QTL_GENE_COMPARATOR_MANIFEST.json`.
- Do not overwrite any existing dataset/result folder. If source count or expected QTL counts drift, stop and inspect input.

## Exact execution (when server access is restored)
```bash
set -euo pipefail
cd /srv/is-analysis/worktrees/is-broad-discovery-20261009
git status --short
git fetch origin research/is-broad-discovery-20261009
# Preserve unrelated or concurrent untracked files before integrating remote changes.
python3 -m unittest discover -s tests -p 'test_is_priority4_qtl_gene_comparators.py' -v
python3 scripts/is/audit_is_priority4_qtl_gene_comparators.py \
 --source /srv/is-analysis/results/is/audits/is_priority4_fullcis_bbj_giga_snp_exact_20261011_v1/IS_PRIORITY4_FULLCIS_EBI_BBJ_GIGASTROKE_SNP_EVIDENCE.tsv \
 --out /srv/is-analysis/results/is/audits/is_priority4_qtl_gene_comparison_20261018_v1
python3 -m unittest discover -s tests -p 'test_is_*.py' -q
```

## Completion criteria
1. Audit code runs without invalid allele orientation or row-count mismatch, output folder contains exactly 3 expected artifacts, source digest and date/version recorded.
2. Audit summaries give the expected gene-wise eQTL nominal-hit numbers (3, 0, 5, 0) and each hit's locus, chromosome/build, beta/SE/p and GWAS evidence.
3. Integrate result with the **separately executed** P0.1 exact 3-SNP ledger and P0.2 3×3 SNP-pair LD matrix. Do not infer 3 distinct causal SNPs from three significant tags.
4. Compare CALHM2/NEURL to FGF5/C4orf22 and update manuscript Figure and scorecard with evidence-level labels.
5. Execute and report the complete IS regression suite, commit verified results and check research web status if deployed.

**Current blocker:** Remote Desktop Commander monthly command quota; GitHub code creation itself is **not** a successful experimental run. Keep original 2,225 broad IS gene universe, 646 archived coloc and 120 exploratory ABF sensitivity records unchanged.
