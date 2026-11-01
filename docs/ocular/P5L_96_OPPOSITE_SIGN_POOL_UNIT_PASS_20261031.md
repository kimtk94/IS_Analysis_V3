# Ocular P5l — external opposite-sign DEG pool-unit crosscheck (VERIFIED 2026-10-31)

## Original evidence and purpose
- [P5l GitHub Actions #41158660082](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41158660082) **PASS**, [artifact #13025744482](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41158660082/artifacts/13025744482).
- Input 1: P5j [Actions #41158343447](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41158343447), original Mehrotra et al. *IOVS* 2026 TM 100 nM DEX 48h published Table 4 and GSE313736 BeamA sc-pseudobulk 500nM 48h log2FC: **96 genes FDR<0.05 in both but opposite signed directions**.
- Input 2: P5d [Actions #41090827791](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41090827791), conservative **9 paired 10x pool-library** edgeR model, testing 21,081 genes without individual donor pseudoreplication.
- `scripts/ocular/ocular_p5l_discordant_poolunit.py` matches exact Ensembl IDs, verifies both-significance/opposite-sign condition, and writes gene-level aggregate-only results. **No donor identifiers or individual count matrices**.

## Verified empirical results
| Gate | Count |
|---|---:|
| Original P5j opposite-sign genes with FDR<0.05 in both publications | **96** |
| Present in P5d 9-pool experimental-library model | **96/96** |
| Pool model log2FC matches primary BeamA direction | **95/96** |
| Pool model BH-FDR<0.05 | **64/96** |
| Pool model BH-FDR<0.05 **and** primary direction retained | **64/96** |
| Primary BeamA abs(log2FC)>=1 | **3/96** |
| Primary BeamA abs(log2FC)<0.25 | **33/96** |
| Only pool-unit direction flip vs primary | **BMPR2** (+0.0797 primary vs -0.00248 pool model, pool FDR=0.979; *not a stable reversal*) |

Examples of directions independently confirmed by 9-pool pooled-library model within original study: **ERFE**, **PARM1**, **ADAMTS15**, **DEPDC1B**, **PYGM**, **CHAC1**, **EGR3**, **NR4A3**, **GASK1B**, **BAALC**. These remain **opposite to the published external 2026 bulk TM effect direction**; this is a reproducible *between-study discrepancy*, not a molecular novelty or clinical outcome result.

## Interpretation and constraints
- 64/96 opposite-sign genes remain statistically significant in the conservative 9-pool library-unit model, so the discordance is not trivially explained by treating cells or donors as independent 10x libraries.
- Conversely, 33/96 originally have small absolute log2FC<.25; only 3/96 have abs(log2FC)>=1. The direction of such small effects is especially sensitive to expression context, dose, bulk-vs-single-cell cell-state differences and technical normalization.
- **Neither study** separates DEX exposure from all potentially correlated treatment-specific culture/library preparation differences. The external 2026 TM study and GSE313736 are **different donor cohorts and methods**. Published aggregate concordance does not provide sample-level matched-dose, randomized-reprep validation.
- The findings can motivate **dose-/time-/cell-state response heterogeneity**, particularly ERFE/PARM1/CHAC1, but no mechanism is established by signed gene-level results alone.

## Next exact research acceptance
1. Prioritize a smaller mechanistic hypothesis set from 64 direction-discordant pool-unit significant genes **only** after existing literature / cell-type expression / dose-time checks; avoid overfitting to FDR.
2. Obtain matched biological TM cultures with 100/500nM DEX and independent vehicle controls, matched time point and preparation/batch technical replicates, for orthogonal qPCR/proteomic/phenotypic validation. A new experiment is required to demonstrate genuine dose/cell-state reversal.
3. Continue external public raw-data matching only after verifying per-run treatment/time/biological donor (ENA PRJNA1243143 record alone does not supply enough attributes).
4. **STOP** clinical steroid-induced IOP predictor claims without linked clinical ΔIOP outcomes.

**Decision:** GO actual *between-study sign-discordance* and conservative within-study pool-unit robustness for 64 genes; HOLD biological explanation/novelty and external raw-data replication; STOP clinical prediction.
