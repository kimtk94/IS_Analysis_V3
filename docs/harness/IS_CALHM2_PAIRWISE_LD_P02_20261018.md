# IS Harness P0.2 — CALHM2 three-SNP pairwise LD

**2026-10-18 KST · Status: CODE_COMMITTED / NOT_RUN_ON_RESEARCH_SERVER**

## Purpose
Extend the earlier `scripts/is/audit_is_calhm2_three_qtl_snps.py` ledger by measuring actual SNP-to-SNP genotype correlations in 1000 Genomes Phase 3. Correlations to the same lead **do not determine pairwise correlation between the three QTL SNPs**. This distinction is essential before calling them one LD cluster.

## Inputs
- Three-SNP QTL evidence: `/srv/is-analysis/results/is/audits/is_priority4_fullcis_bbj_giga_snp_exact_20261011_v1/IS_PRIORITY4_FULLCIS_EBI_BBJ_GIGASTROKE_SNP_EVIDENCE.tsv`
- 504-person full-cis PLINK ALT-dosage export: `/srv/is-analysis/results/is/audits/is_priority4_eur503_eas504_conditional_20261011_v1/BBJ_IS_L002.EAS504.fullcis.raw`
- 1000G population panel: `/srv/is-analysis/data/is/ld_reference/1kg_eas_grch37/integrated_call_samples_v3.20130502.ALL.panel`

## New code
- `scripts/is/audit_is_calhm2_pairwise_ld.py`: validates three allele-matched CALHM2 QTL SNPs at nominal p<1e-4; checks ALT dosage; selects JPT104 inside EAS504 and disjoint otherEAS400; calculates all 10 distinct pairwise correlations among 3 QTL and 2 stroke lead SNPs for 3 groups (30 rows).
- `tests/test_is_calhm2_pairwise_ld.py`: synthetic 504-person inputs, positive genotype and provenance check, rerun non-overwrite gate, three-SNP cardinality error gate.
- Output: `/srv/is-analysis/results/is/audits/is_calhm2_pairwise_ld_20261018_v1/CALHM2_THREE_SNP_PAIRWISE_LD.tsv` and `CALHM2_PAIRWISE_LD_MANIFEST.json`.

## Server execution when authorized environment is available

```bash
set -euo pipefail
cd /srv/is-analysis/worktrees/is-broad-discovery-20261009
git status --short
git fetch origin research/is-broad-discovery-20261009
# Rebase only after protecting any concurrent untracked or modified files.
python3 -m unittest discover -s tests -p 'test_is_calhm2_pairwise_ld.py' -v
python3 scripts/is/audit_is_calhm2_pairwise_ld.py \
  --evidence /srv/is-analysis/results/is/audits/is_priority4_fullcis_bbj_giga_snp_exact_20261011_v1/IS_PRIORITY4_FULLCIS_EBI_BBJ_GIGASTROKE_SNP_EVIDENCE.tsv \
  --genotypes /srv/is-analysis/results/is/audits/is_priority4_eur503_eas504_conditional_20261011_v1/BBJ_IS_L002.EAS504.fullcis.raw \
  --panel /srv/is-analysis/data/is/ld_reference/1kg_eas_grch37/integrated_call_samples_v3.20130502.ALL.panel \
  --out /srv/is-analysis/results/is/audits/is_calhm2_pairwise_ld_20261018_v1
python3 -m unittest discover -s tests -p 'test_is_*.py' -q
```

## Interpretation requirements and success criteria
- Confirm **exactly three** QTL SNPs after genomic allele harmonization, and 10 pair comparisons × 3 population groups.
- Report the 3×3 QTL-only r² matrix separately from each QTL SNP's correlations with both stroke leads. If a target is monomorphic or missing, block rather than impute.
- A fully connected r²≥0.8 group is an operational high-LD cluster, **not** a SuSiE credible set nor three independent functional signals.
- Summarize how strong EAS/JPT SNP-to-SNP LD compares to reference-LD uncertainty; do not claim study-specific QTL or BBJ cohort LD.
- Require SNP allele/strand cross-check, tests, versioned TSV/manifest and manuscript-grade table before `P0.2=COMPLETE`.
- `NEURL` negative comparator and `FGF5` positive comparator are still follow-up work, not included in the pairwise script.
- Preserve existing 2,225 candidate genes, 646 historical coloc runs, 120 ABF sensitivity results. No GCTA, SuSiE or formal colocalization inference was executed.

**Why not executed now:** authorized Remote Desktop Commander returned a monthly-quota stop. Its instructions explicitly prohibit retry/reconnect; GitHub code commits can still be made but are not evidence that the server analysis or tests ran.
