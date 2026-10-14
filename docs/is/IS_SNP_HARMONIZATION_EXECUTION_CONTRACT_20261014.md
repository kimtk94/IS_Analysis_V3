# IS SNP-level harmonization execution contract — 2026-10-14

## Scope
Run the versioned `scripts/is/audit_gwas_qtl_snp_harmonization.py` only on **one region, one GWAS phenotype, one molecular-QTL study/tissue and one ancestry at a time**. The public 149-gene readiness table is a queue, not a SNP-level summary-statistics source.

## Required standardized input
Tab-delimited GWAS and QTL tables each require:
`chrom pos ref alt effect_allele other_allele beta pvalue`.

Record, alongside the input tables:
- source file and SHA256, accession, study, phenotype/tissue, ancestry and genome build;
- genome build harmonization/liftover provenance, if performed upstream;
- QTL source ID and region ID;
- sample size, sample overlap, effect scale, and LD reference source/ancestry;
- the mapping from original source column names to normalized names.

**Do not assume** that variant positions are in GRCh37 or that `beta` is harmonized. Do not combine EUR and EAS molecular-QTL effects into one row. Do not overwrite canonical GWAS/QTL files.

## Commands (Codex terminal / research server; no set -e)
```bash
cd /srv/is-analysis/IS_Analysis_V3 || exit 1
python3 -m unittest discover -s tests -p 'test_is_*' -v
python3 scripts/is/validate_public_nominal_qtl_readiness.py \
  docs/data/is/IS_2225_GENES_PUBLIC_NOMINAL_QTL_READINESS.tsv
python3 scripts/is/build_public_nominal_qtl_region_queue.py \
  docs/data/is/IS_2225_GENES_PUBLIC_NOMINAL_QTL_READINESS.tsv \
  --out /srv/is-analysis/results/is/audits/public_nominal_qtl_region_queue_20261014
# Only after the two real, provenance-verified, normalized SNP inputs exist:
python3 scripts/is/audit_gwas_qtl_snp_harmonization.py \
  --gwas /PATH/TO/ONE_REGION_GWAS_NORMALIZED.tsv \
  --qtl /PATH/TO/ONE_REGION_QTL_NORMALIZED.tsv \
  --gwas-build GRCh37 --qtl-build GRCh37 \
  --out /srv/is-analysis/results/is/audits/one_region_one_qtl_allele_audit
```
The example input paths are placeholders, **not** claims that such files exist. If Codex runs on a different machine, SSH to the research server or use its actual local checkout path.

## Scientific gates
- `MATCHED` and `SWAPPED` denote same-reference biallelic SNPs with an aligned effect allele. They are **not** LD validation.
- `PALINDROMIC_UNRESOLVED`, `STRAND_UNRESOLVED`, `MISMATCH` must not enter coloc without additional independent resolution.
- Single-position keying is deliberately fail-closed for duplicate/multiallelic coordinates; normalize upstream with a documented variant key before broader use.
- Check per-variant GWAS/QTL allele frequencies and effect scales before coloc; frequency and SE are not yet in the minimal schema.
- For multi-signal loci, select a validated conditional/SuSiE workflow with ancestry-matched LD; do not infer causality from a single nominal-QTL p-value.
- No candidate gene is promoted to causal status from this audit.

## Status
- Source code: COMMITTED
- Source-level test suite: COMMITTED
- GitHub Actions for commit 6e78ee14: no associated runs returned at inspection
- Real server execution: NOT VERIFIED
- Real GWAS–QTL SNP matching: NOT VERIFIED
