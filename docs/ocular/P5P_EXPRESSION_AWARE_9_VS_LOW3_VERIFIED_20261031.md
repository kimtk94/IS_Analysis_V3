# Ocular P5p — expression-aware external-discordance shortlist (2026-10-31)

## Verified evidence
[GitHub Actions #41161006241 PASS](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41161006241), [artifact #13027111338](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41161006241/artifacts/13027111338). Pure gene-level aggregate crosswalk of P5b primary 57 paired-donor edgeR [#41085982185](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41085982185) and P5m official published 100nM TM 1/6/48h [#41160429359](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41160429359). No donor-level information read or written.

## Direct issue
The original Tier-A hypothesis set in P5m (ERFE, CHAC1, PARM1, |log2FC|>=1) is **low-expression in our BeamA pseudobulk**, and thus the largest log2FC alone is a poor priority rule.

| High-effect target | Mean logCPM, primary 57-donor | Primary BeamA 500nM log2FC | External TM 100nM 48h log2FC | Pool-unit primary FDR |
|---|---:|---:|---:|---:|
| ERFE | **0.2394** | -1.2066 | +0.8044 | 2.06e-5 |
| CHAC1 | **1.3315** | -1.1469 | +0.5751 | 1.03e-4 |
| PARM1 | **1.0045** | +1.1430 | -0.7093 | 0.00521 |

All three still meet gene-count filtering and are statistically significant; **low logCPM is a QC risk flag, not a rejection of the result**, and their measured effect can be real. However, a validated treatment/protein assay must check sufficient starting expression and technical measurement reliability before mechanistic prioritization.

## Exploratory adequate-expression nine-gene set
Filters *defined for P5p exploratory screening*, not independent inferential claims: primary mean logCPM >=2, primary abs(log2FC)>=0.5, official independent 100nM TM 48h abs(log2FC)>=0.5, pooled-library edgeR abs(log2FC)>=0.5, pooled-library BH-FDR<0.05 and same primary sign; P5j source genes also had opposite external signs and were FDR<0.05 in both studies.

| Gene | Primary mean logCPM | Primary BeamA 500nM 48h log2FC | Conservative 9-library-pool log2FC | External TM 100nM 48h log2FC |
|---|---:|---:|---:|---:|
| **DDIT4L** | 5.30 | +0.774 | +0.717 | -0.572 |
| **ADAMTS15** | 2.50 | +0.771 | +0.601 | -1.068 |
| **LPAR1** | 7.22 | +0.746 | +0.670 | -0.660 |
| **STK17B** | 6.91 | +0.745 | +0.684 | -0.531 |
| **EGR3** | 4.29 | +0.740 | +0.790 | -0.825 |
| **ABHD3** | 3.52 | +0.711 | +0.624 | -0.656 |
| **RAB32** | 7.05 | -0.648 | -0.646 | +0.601 |
| **TFPI** | 6.82 | +0.615 | +0.567 | -0.723 |
| **SDC1** | 4.43 | -0.554 | -0.726 | +0.631 |

The source counts are **9 adequate-expression**, **3 high-effect low-expression** and **84 remaining direction conflicts**. All 96 are in-vitro expression hypotheses; none is proven clinically causal. These 9 can include previously known cellular steroid biology, so **they are not nine novel genes**.

## Important time-course subtleties
- P5m official external 100nM TM results are 1h, 6h and 48h. Early 1h/6h effects are often nonsignificant even if raw log2FC sign differs from 48h. Time or dose-dependent **reversal must not be claimed** without a proper interaction test in independently prepared matched cell cultures.
- Our primary GSE313736 500nM BeamA and published external IOVS 100nM bulk TM use different donors, protocols, and cell state definitions, so one cannot infer concentration effect by comparing their published logFCs.
- Prior P5b nine formal QL leave-one-pool-out fits were **restricted to original top100** and covered **0/96** conflict genes. The P5d 9-source-pool experimental-unit model does cover all 96 and found 64 significant, but does not substitute for a formal 9x leave-out QL refit of these genes.
- P5o has been initiated to compute 96-gene **descriptive** donor-paired TMM log2CPM deltas, 9 pool-average signed effects and 9 leave-pool-out signed mean effects with a cluster bootstrap. This is **not** equivalent to 9 robust edgeR fits or independent validation. Any claims need its actual completed artifact.
- Technical treatment-vs-library preparation confounding remains. Post-treatment BeamA high-confidence cell selection is another source of potential bias (see P5e).
- No clinical steroid-induced ΔIOP phenotype for these cultured cells; clinical glaucoma response prediction remains STOP.

## Next decision
If expression-aware candidates and Tier-A low-abundance candidates retain uniform 9-pool signed effects in P5o, prioritize **orthogonal qPCR/protein** across independently cultured human TM donors using balanced Vehicle/100nM/500nM at 1h/6h/48h, with randomized library preparations and predeclared assay quality. No mechanistic/clinical GO before those independent tests.
