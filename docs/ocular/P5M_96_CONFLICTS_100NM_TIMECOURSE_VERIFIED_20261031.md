# Ocular P5m: official 100nM TM time-course and 96-sign-conflict hypothesis tiers — VERIFIED 2026-10-31

## Provenance and exact dataset
- **PASS GitHub Actions [#41160429359](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41160429359)**, [aggregate artifact #13027385114](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41160429359/artifacts/13027385114): `P5m_96_conflicts_timecourse_priorities.tsv` and `P5m_timecourse_tier_summary.json`.
- Publication: Mehrotra et al. *IOVS* 2026, DOI 10.1167/iovs.67.2.13, PMCID PMC12882097. Official unmodified Supplementary XLSX `iovs-67-2-13_s003.xlsx` via Europe PMC official `supplementaryFiles`; **TM 1h** `SuppTable64_TM1Hr_All_DEGs` (17,077 genes), **TM 6h** `SuppTable41_TM6Hrs_All_DEGs` (17,132 genes), **TM 48h** `SuppTable4_TM2Days_All_DEGs` (13,313 genes). These are 100 nM DEX published gene-level DESeq2 contrasts (not reprocessed FASTQ).
- Input original GSE313736 57 paired high-confidence BeamA donor pseudobulks at nominal 500nM 48h vs Vehicle, corrected edgeR primary [#41085982185](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41085982185); conservative 9-paired-source-library pool model [#41090827791](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41090827791); 96 genes called opposite-sign BH-FDR significant in both studies [P5l #41158660082](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41158660082).
- Code `scripts/ocular/ocular_p5m_timecourse_and_tiers.py`, orchestration `.github/workflows/ocular-p5m-timecourse-tiers.yml` (manual only after successful execution). No donor IDs, individual-level expression data, or clinical phenotypes saved.

## Exact 96-gene hypothesis tier results (exploratory, post-hoc)
| Tier (based on primary and external effect size plus pool-level significance) | Genes |
|---|---:|
| A: 9-pool model BH-FDR<0.05 with matching primary sign, primary abs(log2FC)>=1 and external abs(log2FC)>=0.5 and pool abs(log2FC)>=0.5 | **3** |
| B: as above but primary abs(log2FC)>=0.5 (and not Tier A) | **11** |
| C: smaller but stable: primary abs(log2FC)>=0.25 and external >=0.5 and pool >=0.25 (and not Tier A/B) | **13** |
| D: pool-level FDR<0.05 but weaker primary/external magnitude than C | **37** |
| E: pool-level FDR not significant or sign differs from primary | **32** |
| **Total** | **96** |

Tiers are *hypothesis-prioritization thresholds chosen after the original P5 results*, **not** new inferential significance tiers or preregistered claims. Tier A **ERFE, CHAC1, PARM1**; Tier B 11 examples PYGM, DDIT4L, ADAMTS15, LPAR1, STK17B, EGR3, DEPDC1B, ABHD3, RAB32, TFPI, SDC1.

## Strongest three: 100nM published time course vs 500nM primary
| Gene | Primary 500nM BeamA 48h log2FC | Primary FDR | Pool-unit 48h log2FC | Pool-unit FDR | External TM 100nM 1h log2FC/FDR | External 6h log2FC/FDR | External 48h log2FC/FDR |
|---|---:|---:|---:|---:|---:|---:|---:|
| **ERFE** | -1.206625 | 1.54e-14 | -1.430800 | 2.06e-5 | +0.170374 / 0.999988 | +0.863640 / **0.04026** | +0.804425 / **0.01242** |
| **CHAC1** | -1.146898 | 1.01e-12 | -1.253334 | 1.03e-4 | -0.298936 / 0.999988 | -0.305432 / 0.999628 | +0.575081 / **0.01174** |
| **PARM1** | +1.143010 | 3.59e-9 | +1.155957 | 0.00521 | -0.088726 / 0.999988 | +0.374992 / 0.999628 | -0.709311 / **0.000614** |

ERFE is significantly up at both 6h and 48h in external 100nM bulk TM but decreased in our 500nM selected BeamA. CHAC1 external 1/6h signs are negative but not statistically significant, switching to positive at 48h; PARM1 external 6h positive is also nonsignificant, 48h negative is significant. **Do not treat raw sign flips at nonsignificant early times as demonstrated within-study regulatory reversal.**

## Time-course coverage and pooled patterns of all 96 discordances
- All **96/96** measurable at 1h and 48h; **95/96** measurable at 6h.
- Relative to external 48h effect sign, **50/96** have opposite raw 1h signs, **50/95** have opposite raw 6h signs. Only **2** genes FDR-significant at 1h and **3** at 6h among available 96, so **raw sign flips are mainly statistical noise/weak signals until further testing**, not a confirmed dynamic process.
- **Critical gap:** P5b's 9 formal edgeR QL leave-one-pool-out refits were restricted to its original **top 100** FDR-ranked genes. The current 96 external opposite-direction genes overlap **0/96** with that P5b top100. It was incorrect to infer that the 96 have each passed formal 9-pool exclusion testing. P5d 9-paired-library-unit edgeR did cover 96/96 and found 64 FDR-significant, but this is a **different sensitivity model**.

## Correct interpretation
1. **GO** signed gene-level externally reproduced contrasts and hypothesis tiers with actual independent publication's 100nM gene expression trajectories.
2. **HOLD** dose-dependent reversal: primary 500nM vs external 100nM is confounded with 10x vs bulk, donor cohort and BeamA selection; within-external 1/6/48h contrasts share 100nM dose and have no significant early effect for CHAC1/PARM1.
3. **HOLD** targeted 9-pool dropout/biological pairing stability of 96 (P5b formally tested none); design a new 96-target analysis with proper experiment units. The 9-pool model is relevant but does not replace leave-one-pool-out testing.
4. **HOLD** new mechanism: CHAC1 is a glutathione/ER-stress-response gene in general biology; TM GC-induced ER stress is well established (Kasetti et al. *Scientific Reports* 2017; DOI 10.1038/s41598-017-14938-0). ERFE and PARM1 functions and TM effects need independent literature and protein validation. No causal mechanism was tested.
5. **STOP** clinical steroid-induced IOP/genetic risk prediction from these culture-only data with no matched clinical delta-IOP.

## Next strict analysis
Run a **targeted, transparent** 96-gene donor-paired pool-sign and leave-one-pool-out *descriptive* sensitivity using raw counts/normalization, and for Tier-A genes assess pool-level confidence intervals. The descriptive check **must not** be labeled a nine-fold edgeR QL refit or independent cohort replication. Fully controlled biological validation would require matched TM cultures at 0/100/500nM DEX and preplanned 1h/6h/48h sampling across independently prepared libraries.
