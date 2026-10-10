# IS Phase10/11 source provenance & SNP-level replay (2026-10-10)

Status: SOURCE_CODE_PRIORS_VERIFIED; 8/646 SNP_REPLAY_PASS; CAUSAL_CLAIM_BLOCKED.

## Source
Archived Colab: https://drive.google.com/file/d/1HM2j0kS_Twk5zBF0DWZPpafy-Xk4k0iu/view

The original notebook explicitly specifies coloc.abf(p1=1e-4,p2=1e-4,p12=1e-5), binary GWAS N=174686, cases=22664, quantitative QTL N taken from the first qtl_n_scalar entry. This resolves prior *code* provenance; the original run logs/package versions remain unavailable. The input index and raw master on Google Drive have SHA256 1fe044f1ffd8e9aa6b678718b6772063bd11edd462d7bede4104843cdb0c7414 and 0160a01fc6433f9bec839113ac3a2390927fc07db9252d014b0cc4b6b11716ea respectively.

All 646 index rows are READY and match 646 original PASS rows one-to-one by locus/dataset/gene. SNP count and chosen QTL N agree for all 646. Every input has SNP-specific QTL N variation (median max-minus-min 25; IQR 17–30; maximum 81); the scalar N is a historical approximation.

## Initial three-pair SNP replay (retained audit)
Eight original per-SNP GWAS/eQTL inputs fetched from Drive and SHA256-verified on server. Original first three: FGF5 cerebellar hemisphere 1694 SNPs, PP.H4=0.7786969180; SH3PXD2A tibial artery 1583, PP.H4=0.3508155847; COL4A2 brain putamen 3015, PP.H4=0.3307779941; five additional pairs are detailed in the extended table below. All original H0–H4 vectors reproduced within 1e-8. **120 SNP-level runs completed: 8 gene–tissue pairs x 5 p12 values x 3 scalar-QTL-N assumptions.** Other **638 pairs** have NOT been independently replayed.

QTL minimum N sensitivity at p12=1e-5: FGF5 0.777030; SH3PXD2A 0.350831; COL4A2 0.313262. This is an assumption-sensitivity diagnostic, not a SNP-specific sample-size model.

## Cross-ancestry MAF diagnostic
A SNP-level sensitivity excluded variants with |GWAS MAF minus GTEx QTL MAF| >0.2 or >0.1. Results are *exploratory* and are NOT correct ancestry harmonization.

| Gene | Original SNPs | MAF difference >0.2 | Original H4 | H4 with <=0.2 | H4 with <=0.1 |
|---|---:|---:|---:|---:|---:|
| FGF5 | 1694 | 444 | 0.7787 | 0.8332 | 0.8705 |
| SH3PXD2A | 1583 | 121 | 0.3508 | 0.3535 | 0.3770 |
| COL4A2 | 3015 | 331 | 0.3308 | 0.3418 | 0.3645 |

At <=0.1, the strongest eQTL variant is **lost** for FGF5 and COL4A2, despite the leading GWAS variant remaining. Higher H4 after a post hoc SNP filter therefore cannot be claimed as stronger evidence. Large MAF discrepancies are expected for some alleles when EAS GWAS and GTEx reference ancestries differ and do not independently establish allele error.

## Reproducibility
- R source rerun: scripts/is/replay_is_phase10_abf_snp_inputs.R
- MAF diagnostic: scripts/is/audit_is_legacy_crossancestry_maf.R
- SNP-by-SNP source input audit: results/is/audits/legacy_coloc_snp_replay_20261010_v1/
- Original files unchanged; separate output subdirs replay/ and maf_qc/.
- Source original archived notebook SHA256 c8695136c62002e1f5e362a561e32f2c6e7e217bb966ce815d11992d9a908945.
- Baseline full prior sensitivity 646 x 5 = 3230 rows from Phase9C remains valid under code-confirmed original priors, but *not* a variant-level rerun of all 646 pairs.

## Remaining publication blockers
1. Recover authenticated runtime logs/version for original 646-test source.
2. Rerun the remaining 638 input assays in separate storage and include multiplicity, tissue selection, and failed-pair denominator.
3. Audit GWAS and GTEx effect-allele/build matching, source GWAS ancestry, SNP coverage, sample-specific missingness, and population-matched LD. A MAF filter alone cannot solve this.
4. Use LD-aware multi-signal coloc with cohort-matched QTL LD and control for heterogeneity/pleiotropy.
5. Obtain truly orthogonal disease-cell regulatory or functional evidence: R3_7 healthy normal brain donor CPM cannot validate a stroke-specific GWAS mechanism.
6. Keep all 80 expanded provisional genomic windows and 2225 unique positional gene IDs eligible, not narrowed solely by this historical four-locus screen.

**Interpretation:** This stage validates computational reproducibility and source prior settings, not any causal target or disease-specific cellular mechanism.

## Direct SNP-level posterior vs summary-only prior adjustment

For the initial three-pair subset, each SNP input was run at five p12 values (with original FIRST scalar QTL N) and compared against coloc 5.2.3 prior.adjust outputs for the **same gene/locus/tissue**.

- 3 genes × 5 prior settings = 15 matched scenarios and 75 H0–H4 posterior comparisons.
- At the original p12=1e-5, all H0–H4 posteriors match within 1e-8.
- At changed priors, small numerical discrepancies occur: maximum absolute difference over the 75 H0–H4 values **5.69513e-5**; largest H4-specific difference 4.68054e-5 (COL4A2 at p12=1e-4).
- The 646 × 5 prior.adjust grid is therefore *approximately* concordant with independent per-SNP coloc recomputation for the initial three replayed pairs, not numerically identical at every changed prior.
- Prefer directly recomputed SNP-level posterior when both exist; do not interpret posterior-reweighted H4 as if every historical pair had a SNP-level rerun.

Script: `scripts/is/audit_is_snp_vs_summary_prior_grid.py`.
JSON output: `results/is/audits/legacy_coloc_snp_replay_20261010_v1/IS_DIRECT_VS_SUMMARY_PRIOR_GRID_AUDIT.json`.


## Extended independent 8-pair replay — supplemental (2026-10-10)

Eight prespecified gene/tissue pairs were reproduced using original shared-SNP GWAS+eQTL beta/SE/MAF; coloc 5.2.3, original p1/p2/p12 and BBJ sample specification.

| Gene | Tissue | Shared SNPs | Original H4 = SNP replay | H4 with minimum scalar QTL N |
|---|---|---:|---:|---:|
| FGF5 | Cerebellar hemisphere | 1694 | 0.778697 | 0.777030 |
| SH3PXD2A | Artery tibial | 1583 | 0.350816 | 0.350831 |
| COL4A2 | Brain putamen | 3015 | 0.330778 | 0.313262 |
| CALHM2 | Cerebellar hemisphere | 1583 | 0.797684 | 0.766643 |
| NEURL1 | Cerebellar hemisphere | 1583 | 0.677586 | 0.665234 |
| C4orf22 | Cerebellar hemisphere | 1694 | 0.659936 | 0.655550 |
| INA | Anterior cingulate BA24 | 1581 | 0.536369 | 0.493376 |
| COL4A1 | Brain amygdala | 3016 | 0.117146 | 0.114673 |

The 8-pair grid has 120 direct SNP-level `coloc.abf` runs, all eight historical H0–H4 posterior vectors exactly reproduced within 1e-8 at original settings. `FIRST/MIN/MAX` are scenario assumptions for scalar QTL N; these are **not** correct per-SNP effective sample sizes. INA, CALHM2 and COL4A2 are notably sensitive to this assumption.

Summary-only `prior.adjust` 8-pair matched grid (40 direct SNP scenarios using historical QTL FIRST N) agrees with actual per-SNP rerun within maximum absolute H0–H4 difference of **0.000157301** across nonbaseline priors. Direct SNP calculations take precedence.

MAF-delta diagnostics (8 x 3 = 24 scenarios) condition on cross-population differences. At |EAS GWAS MAF − GTEx QTL MAF| <= 0.1, leading eQTL SNPs were lost for **FGF5, COL4A2, NEURL1 and C4orf22** despite H4 increasing. Therefore increased H4 after post hoc filtering is *not* independent evidence; removed SNP sets change the model and posterior denominator.

**Source/version boundary:** archived source code p1/p2/p12 now recovered and eight targeted source SNP-replays pass; the original execution runtime is not attested; 638 pairs remain outside SNP replay. No expanded 80-window causal claim supported.

The separate audited 8-gene ledger in `results/is/audits/is_priority8_source_verified_20261010_v2/` logs SHA256 checksums of **all eight original per-SNP source TSV files**, the posterior vector verification, 120-run scenario count, MAF filtering caveats and the exact 638-pair unverified denominator. Raw SNP inputs are never included in Git.

Extended script names:
- `scripts/is/replay_is_phase10_abf_priority8_snp_inputs.R`
- `scripts/is/audit_is_legacy_crossancestry_maf_priority8.R`
- `scripts/is/audit_is_priority8_snp_vs_summary_prior_grid.py`

Separated immutable output: `results/is/audits/legacy_coloc_snp_replay_priority8_20261010_v1/`.

## Source-frozen eight-pair scientific reporting figures and tables

The priority-eight SNP-level report was generated independently of the original results and includes per-input SHA256 provenance and the full five-H0–H4 replay check for each selected pair:

- [Figure — direct SNP-level PP.H4 under p12 sensitivity](figures/is/FIG_IS_PRIORITY8_DIRECT_SNP_P12.svg).
- Audit output: `results/is/audits/is_priority8_genetics_summary_20261010_v1/IS_PRIORITY8_SOURCE_VERIFIED.json`.
- Eight gene–tissue summary: `results/is/audits/is_priority8_genetics_summary_20261010_v1/IS_PRIORITY8_SOURCE_VERIFIED.tsv`.
- Human donor cross-evidence (9-gene status ledger): `results/is/audits/is_p0_integrated_genetics_cell_20261010_v3/`.

All eight SNP-level analyses numerically replicate the original baseline posterior. The plot shows five *direct* SNP-level H4 values per gene, rather than summary-only p12 posterior adjustment. Note these selected gene–tissue pairs came from the earlier four-locus GTEx screen and do not replace broad 80-window discovery, causal variant identification, or multiple-testing correction.

The original code settings are recovered from the archived v3 notebook; its runtime was not attested. There are **638 historical tests remaining** for independent SNP-level replay. A complete original 646-test result table is not equivalent to 646 independent source reruns.
