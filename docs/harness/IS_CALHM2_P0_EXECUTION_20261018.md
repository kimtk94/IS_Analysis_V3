# IS Harness P0 — CALHM2 3-SNP QTL audit (execution-ready)

**Prepared:** 2026-10-18 KST. **Status:** `CODE_COMMITTED_SERVER_EXECUTION_PENDING`. Server-side analysis and tests were **not executed in this step** because the Remote Desktop Commander quota is exhausted. The existing verified checkpoint is commit `9c1a641` with 329 passing IS tests; later report `4f7974d` contains the original Harness plan. Do not claim newly calculated CALHM2 SNP coordinates or LD from this prepared code until real execution.

## Mission
Generate a **three-SNP variant-level ledger**, not another inferred causal statement. Extract all exact-allele non-palindromic CALHM2 EBI QTD000161 SNPs with nominal p<1e-4 on chr10, with GRCh37/38 positions, QTL beta/SE/p and BBJ/GIGASTROKE beta/p. Join each SNP to 1000G Phase 3 EAS504, JPT104, and otherEAS400 signed r/r² with **both BBJ and GIGASTROKE regional leads**. The code checks that all three variants are present and the earlier LD tagging counts (BBJ lead 3, GIGASTROKE lead 0 in EAS and JPT) are reproduced.

## Inputs already on research server
- Evidence: `/srv/is-analysis/results/is/audits/is_priority4_fullcis_bbj_giga_snp_exact_20261011_v1/IS_PRIORITY4_FULLCIS_EBI_BBJ_GIGASTROKE_SNP_EVIDENCE.tsv`
- Previous per-variant LD tags: `/srv/is-analysis/results/is/audits/is_chr10_two_lead_ld_tags_20261011_v1/IS_CHR10_EAS_JPT_TWO_LEAD_PER_SNP_LD.tsv`
- Target **new** output (non-overwrite): `/srv/is-analysis/results/is/audits/is_calhm2_three_qtl_snp_ledger_20261018_v1/`
- Script: `scripts/is/audit_is_calhm2_three_qtl_snps.py`
- Tests: `tests/test_is_calhm2_three_qtl_snps.py`

## Resume commands (run in authenticated server environment)
```bash
set -euo pipefail
cd /srv/is-analysis/worktrees/is-broad-discovery-20261009
git status --short
git fetch origin research/is-broad-discovery-20261009
# IMPORTANT: preserve unrelated untracked files; resolve remote divergence
# without force-push or overwriting concurrent work.
python3 -m unittest discover -s tests -p 'test_is_calhm2_three_qtl_snps.py' -v
python3 scripts/is/audit_is_calhm2_three_qtl_snps.py \
 --evidence /srv/is-analysis/results/is/audits/is_priority4_fullcis_bbj_giga_snp_exact_20261011_v1/IS_PRIORITY4_FULLCIS_EBI_BBJ_GIGASTROKE_SNP_EVIDENCE.tsv \
 --ld /srv/is-analysis/results/is/audits/is_chr10_two_lead_ld_tags_20261011_v1/IS_CHR10_EAS_JPT_TWO_LEAD_PER_SNP_LD.tsv \
 --out /srv/is-analysis/results/is/audits/is_calhm2_three_qtl_snp_ledger_20261018_v1
python3 -m unittest discover -s tests -p 'test_is_*.py' -q
```

## Expected outputs and gates
- `CALHM2_THREE_QTL_SNP_LEDGER.tsv`: 3 records, QTL effects and two-GWAS marginal values plus 6 ancestry×lead signed-LD/r² measures per SNP.
- `CALHM2_THREE_QTL_SNP_MANIFEST.json`: source paths and `EXPLORATORY_VARIANT_LEDGER_NOT_CREDIBLE_SET` evidence designation.
- Fail closed if fewer or more than 3 qualifying molecular SNPs, reference tagging cannot reproduce prior counts, or requested output directory exists.
- Before extending to pairwise-QTL SNP clustering, obtain actual **pairwise ALT genotype** r² from EAS/JPT rather than treating their separate correlations with a lead as the pairwise matrix.
- Add **NEURL** as negative molecular comparator and **FGF5** as positive comparison in a subsequent audited batch; no results generated for those comparators in the current implementation.
- Upon actual execution, check exact SNP coordinates/allele order and results; amend report with values and the execution commit SHA. Test/CI should pass from repository root, not from another current working directory.

## Scientific boundaries and invariants
The three SNPs are not independent causal observations and r²≥0.8 is not a credible set. 1000G EAS504/JPT104 are reference proxies, **not** BBJ/GIGASTROKE/GTEx study LD. The molecular tested sample size and multi-signal coloc remain unverified. No valid COJO, SuSiE multi-signal credible sets, or independent causal replication is generated. Preserve all original **2,225 IS candidate genes**, **646 historical coloc** records and **120 ABF sensitivity** runs. Do not touch unrelated untracked files `scripts/is/audit_is_col4a_eas_eur_ais_cross_ancestry_transfer.py` or `server/is_audit_col4a_eur_ais_v1.sh`.
