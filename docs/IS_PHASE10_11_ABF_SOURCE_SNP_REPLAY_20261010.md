# IS Phase10/11 source provenance & SNP-level replay (2026-10-10)

Status: SOURCE_CODE_PRIORS_VERIFIED; 3/646 SNP_REPLAY_PASS; CAUSAL_CLAIM_BLOCKED.

## Source
Archived Colab: https://drive.google.com/file/d/1HM2j0kS_Twk5zBF0DWZPpafy-Xk4k0iu/view

The original notebook explicitly specifies coloc.abf(p1=1e-4,p2=1e-4,p12=1e-5), binary GWAS N=174686, cases=22664, quantitative QTL N taken from the first qtl_n_scalar entry. This resolves prior *code* provenance; the original run logs/package versions remain unavailable. The input index and raw master on Google Drive have SHA256 1fe044f1ffd8e9aa6b678718b6772063bd11edd462d7bede4104843cdb0c7414 and 0160a01fc6433f9bec839113ac3a2390927fc07db9252d014b0cc4b6b11716ea respectively.

All 646 index rows are READY and match 646 original PASS rows one-to-one by locus/dataset/gene. SNP count and chosen QTL N agree for all 646. Every input has SNP-specific QTL N variation (median max-minus-min 25; IQR 17–30; maximum 81); the scalar N is a historical approximation.

## Direct SNP replay
Three original per-SNP GWAS/eQTL inputs fetched from Drive and SHA256-verified on server. Re-executed coloc 5.2.3 using source notebook settings: FGF5 cerebellar hemisphere 1694 SNPs, PP.H4=0.7786969180; SH3PXD2A tibial artery 1583, PP.H4=0.3508155847; COL4A2 brain putamen 3015, PP.H4=0.3307779941. All original H0–H4 vectors reproduced within 1e-8. 45 SNP-level runs completed: 3 genes x 5 p12 values x 3 scalar-QTL-N assumptions. Other 643 pairs have NOT been independently replayed.

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
2. Rerun all selected 646 input assays in separate storage and include multiplicity, tissue selection, and failed-pair denominator.
3. Audit GWAS and GTEx effect-allele/build matching, source GWAS ancestry, SNP coverage, sample-specific missingness, and population-matched LD. A MAF filter alone cannot solve this.
4. Use LD-aware multi-signal coloc with cohort-matched QTL LD and control for heterogeneity/pleiotropy.
5. Obtain truly orthogonal disease-cell regulatory or functional evidence: R3_7 healthy normal brain donor CPM cannot validate a stroke-specific GWAS mechanism.
6. Keep all 80 expanded provisional genomic windows and 2225 unique positional gene IDs eligible, not narrowed solely by this historical four-locus screen.

**Interpretation:** This stage validates computational reproducibility and source prior settings, not any causal target or disease-specific cellular mechanism.

## Direct SNP-level posterior vs summary-only prior adjustment

Each of the three independently replayed SNP inputs was run at five p12 values (with original FIRST scalar QTL N) and compared against coloc 5.2.3 prior.adjust outputs for the **same gene/locus/tissue**.

- 3 genes × 5 prior settings = 15 matched scenarios and 75 H0–H4 posterior comparisons.
- At the original p12=1e-5, all H0–H4 posteriors match within 1e-8.
- At changed priors, small numerical discrepancies occur: maximum absolute difference over the 75 H0–H4 values **5.69513e-5**; largest H4-specific difference 4.68054e-5 (COL4A2 at p12=1e-4).
- The 646 × 5 prior.adjust grid is therefore *approximately* concordant with independent per-SNP coloc recomputation for the three replayed pairs, not numerically identical at every changed prior.
- Prefer directly recomputed SNP-level posterior when both exist; do not interpret posterior-reweighted H4 as if every historical pair had a SNP-level rerun.

Script: `scripts/is/audit_is_snp_vs_summary_prior_grid.py`.
JSON output: `results/is/audits/legacy_coloc_snp_replay_20261010_v1/IS_DIRECT_VS_SUMMARY_PRIOR_GRID_AUDIT.json`.
