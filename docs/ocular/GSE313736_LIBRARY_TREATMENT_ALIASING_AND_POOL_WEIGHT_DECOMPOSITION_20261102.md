# GSE313736: exact DEX–10x-library aliasing and observed BeamA pool-weight decomposition

**Date:** 2026-11-02 (KST). This is an original-study design constraint and descriptive public-aggregate QC; it is **not** an experiment identifying DEX pharmacological causality, nor new clinical evidence.

## Executed evidence and reproducibility

- **Actual [GitHub Actions #41430965326](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41430965326) 2/2 SUCCESS** — independent synthetic algebra/decomposition self-test and actual source aggregate reinspection.
- [Small source-validated aggregate artifact #13133928461](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41430965326/artifacts/13133928461) contains JSON with full per-pool source counts and symmetric decomposition, hash of three source files, source model numbers, and explanatory Markdown report. No donor identities, per-donor raw UMI counts, original 10x matrices, genotype or clinical records were downloaded in this CI.
- Reproducible script: `scripts/ocular/ocular_p5_technical_confounding_identifiability_audit.py`; workflow `.github/workflows/ocular-p5-technical-confounding-identifiability.yml`.
- Input A: actual cell annotation selection summary from [P5e #41091146021](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41091146021), artifact #12994476349.
- Input B: original source pool leave-one-out and anonymous MDS summary [P5b #41085982185](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41085982185), artifact #12993824564.
- Input C: original 9-pool library-unit robust edgeR [P5d #41090827791](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41090827791), artifact #12995175346.

## 1. Exact structural identifiability gate

The original nine source donor pools A–I were separately prepared in **Vehicle** and **DEX** 10x libraries, giving 18 distinct, condition-specific original libraries. The paired-pool model with 9 pool indicator columns plus 1 DEX column has **rank 10/10**, with 18 observations and 8 residual degrees of freedom. It can fit an arm contrast **under the assumed model**.

However, suppose each original 10x library may have its own unrestricted nuisance intercept. An 18-library-indicator matrix has rank 18 across the 18 pool-arm profiles. The DEX indicator vector is **exactly the sum of the nine DEX-specific library indicator columns**. Consequently, adding the DEX column to the 18 unrestricted library columns gives **19 design columns but still rank 18** (rank increment **zero**). No comparison in this source alone can distinguish a DEX-specific drug coefficient from arbitrary condition-associated source-library nuisance intercepts.

This does **not** assert the library differences explain all observed biology. The rank result means that a DEX-specific causal effect is **not identifiable without additional design assumptions or independently balanced library preparation**. Using 57 paired donors or leaving out one of nine pools does not create cross-condition replication in the same original 10x library. Neither high donor-level MDS treatment separation nor successful edgeR FDR values change this algebraic fact.

## 2. Pooled-cell selection proportion versus pool-composition arithmetic

The source P5e metadata had **169,398 SNG cells**, 59 donors with both arms, a BeamA confidence threshold ≥0.75, and full coverage of all nine pools (validated aggregate totals). The original global DEX−Vehicle BeamA fraction difference is **+2.4064 percentage points**, whereas the paired 59-donor median difference is **−0.421 percentage points** (29 positive, 30 negative, two-sided sign-test p=1). These are **different weighted estimands** and not contradictory.

Using *the original published pool-arm SNG counts* and BeamA pool proportions (rounded to six decimals in P5e), the new audit reconstructs the global difference to within the required source-precision tolerance and applies an **exact symmetric Kitagawa arithmetic decomposition** on the published rounded rates:

\[
\sum_p w_{DEX,p}r_{DEX,p}-\sum_p w_{Veh,p}r_{Veh,p}
=\sum_p \frac{w_{DEX,p}+w_{Veh,p}}{2}(r_{DEX,p}-r_{Veh,p})
+\sum_p \frac{r_{DEX,p}+r_{Veh,p}}{2}(w_{DEX,p}-w_{Veh,p}).
\]

| Term | Measured percentage-point contribution |
|---|---:|
| Observed global DEX−Vehicle fraction | **+2.406339 pp** |
| Mean-arm-weighted **within-pool BeamA fraction** differences | **+0.868638 pp** |
| Mean-arm-proportion-weighted **pool capture/composition weights** | **+1.537701 pp** |
| Sum of components minus observed total | **~0 pp** |
| Matched donor-level median | **−0.421 pp**, 29 vs 30 donor signs |

The pool-weight term is **63.90% of the global descriptive difference** (1.537701/2.406339). **It is NOT an estimate that batch artifact causes 63.9% of transcriptional variation or of the drug response.** Cell-capture depth, cell recovery, donor composition and/or post-treatment classifier relabeling may contribute; this decomposition does not separate such mechanisms. Pool-level rates were rounded to 6 decimals in the original P5e archive. For precise per-pool SNG counts, relative pool weights, rates and contributions, use the archived JSON.

## 3. Existing actual edgeR and MDS evidence remains an association

| Original frozen statistic | Reverified value | Causal interpretation |
|---|---:|---|
| Paired donor-model independent donor pairs | **57** | Biological pairing, not 57 distinct 10x preparation pairs |
| Source pool-level paired 10x library profiles | **9 pairs, 18 libraries** | Pool sensitivity, not randomized cross-condition library replicates |
| Donor-paired and pool-unit common tested genes | **16,049** | Shared test universe |
| BH-FDR <0.05 in both on common universe | **11,079/13,716** original significant | Model-unit sensitivity only |
| Donor/pool model log2FC Pearson | **0.99624485** | High within-source numerical agreement |
| Nine leave-one-pool refit minimal Pearson | **0.99726488** | Not independent cross-library control |
| MDS axis 1 regression on DEX alone, unadjusted R² | **0.95557721** | Not `95.6% DEX causal variance explained` |

This is stronger evidence for a **specific study-design inference boundary** than merely stating that "batch may influence treatment." But it does not estimate how much of each gene's apparent DEX log2FC is a true pharmacological effect, because the decomposition of arm signal into drug and library nuisance parameters is not determined by the original experimental design.

## 4. Defensible thesis decision, with no author outreach

- **GO:** source-code rank check, observed pool capture/composition decomposition, paired-vs-pooled statistical sensitivity and explicit provenance. The previously observed candidate direction consistency (MYOC↑, TGFBI↓, MMP2↓) remains a **transcript association** and cross-study prioritization hypothesis.
- **HOLD:** pharmacological DEX causality *within GSE313736* unconfounded by source-library effects, biological selection versus post-treatment cell-label artifact, gene-level absolute effects attributed to DEX, unrelated clinical outcomes.
- **HOLD / non-blocking:** unavailable original author DEG80 list; there is no request to contact investigators.
- **STOP:** clinical steroid-induced ΔIOP risk prediction, direct ECM fibrosis/aqueous outflow phenotype or new molecular mechanism claim.
- **Next feasible task:** integrate the rank and 2.406 pp decomposition into manuscript Methods, Results and central Limitations. Compare any independent 100 nM/500 nM sources only at defensible gene-level estimands; a randomized batch-balanced wet-lab design is needed for stronger mechanistic claims.

**Study novelty ceiling:** transparent identification of technical identifiability boundaries, not claiming a corrected DEX causal effect without balanced controls.
