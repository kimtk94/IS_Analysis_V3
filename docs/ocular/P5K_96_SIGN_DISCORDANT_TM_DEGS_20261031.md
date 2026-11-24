# Ocular P5k — 96 significant sign-discordant TM DEX genes (VERIFIED 2026-10-31)

## Evidence
- This is the **executed** P5j official published bulk human TM versus GSE313736 BeamA pseudobulk gene-level Ensembl-ID comparison. [CI PASS #41158343447](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41158343447), aggregate [artifact #13026128279](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41158343447/artifacts/13026128279).
- `P5j_sign_discordant_both_FDR_genes.tsv` contains exactly **96** genes with BH-FDR<0.05 in both studies **but opposite treatment log2FC signs**. This file is derived quantitative research evidence, **not** a statement that either experiment is invalid.
- Official source: Mehrotra et al., *IOVS* 2026 PMID 41636427, PMCID PMC12882097 `iovs-67-2-13_s003.xlsx` worksheet `SuppTable4_TM2Days_All_DEGs` (13,313 original gene results); human TM 100nM DEX at 48h. Primary GSE313736: BeamA-scRNA donor-paired edgeR, nominal 500nM DEX at 48h, 57 pairs.

## High-priority direction-discordant genes (largest log2FC gap)
| Gene | Primary BeamA 57 donors log2FC | External TM 2026 bulk log2FC | Primary FDR | External FDR |
|---|---:|---:|---:|---:|
| ERFE | -1.207 | +0.804 | 1.54e-14 | 0.0124 |
| PARM1 | +1.143 | -0.709 | 3.59e-9 | 6.14e-4 |
| ADAMTS15 | +0.771 | -1.068 | 8.98e-13 | 1.31e-4 |
| DEPDC1B | +0.731 | -1.101 | 5.81e-10 | 0.0374 |
| PYGM | +0.783 | -1.027 | 1.02e-11 | 5.15e-4 |
| LCNL1 | -0.329 | +1.442 | 8.19e-4 | 0.0112 |
| CHAC1 | -1.147 | +0.575 | 1.01e-12 | 0.0117 |
| EGR3 | +0.740 | -0.825 | 3.66e-14 | 0.00993 |
| NR4A3 | +0.488 | -1.034 | 1.49e-5 | 0.0298 |
| GASK1B | +0.499 | -1.020 | 3.07e-8 | 5.79e-9 |

All ten and all 96 were mechanically identified by strict gene ID harmonization, **not selected after inspecting the clinical IOP phenotype**. Other relevant contrasting examples in this output include LPAR1, GRIK2, DDIT4L, PCK2 and TFPI.

## Correct interpretation
- Both studies' FDR overlap = **755** of 12,251 common tested IDs; **659** both-significant/same-direction (87.3%), **96** both-significant/opposite-direction (12.7%). Overall logFC Pearson r=0.4613, general same-sign fraction 56.74%.
- Opposite-sign results may reflect true concentration-dependent glucocorticoid response, 100nM vs 500nM, exposure kinetics, broad bulk TM versus selected BeamA cellular substate, donor genetics, different experimental handling/10x libraries and mapping/normalization discrepancies. **Do not assume a specific explanation without a sensitivity test.**
- This opposing-effect class may eventually support an interesting **dose-/cell-state-dependence** story, but only if the direction difference survives technical confounding controls and independent dose/time confirmation. The contrary findings **do not independently validate clinical steroid-induced ΔIOP**.

## Required validation before manuscript GO
1. Compare robust donor-paired logFC and 9-pool-library-unit logFC specifically for the 96; flag any sign flips or weak magnitudes in the conservative model.
2. Confirm gene-ID and gene-symbol mapping, avoid mixed HG38 gene transcripts/ID versions and inspect adjusted p-value missingness.
3. Assess treatment-independent TM annotation cutoffs and pooled label distribution (P5e indicated large global-vs-donor compositional discrepancy).
4. Obtain an independent DEX dose/time dataset with true 100 vs 500nM matched design and independent biological replicates. The 2026 raw ENA run-level treatment map remains unavailable, although the **published 857-gene aggregate table is fully verified and compared**.
5. Do not use the 96-gene conflict set as a novel glaucoma biomarker or steroid-induced IOP prediction panel without functional and clinical validation.

## Gate
**GO** exact cross-study gene-level disagreement audit; **HOLD** causal explanation of discordance, independent raw-data reanalysis and manuscript mechanistic novelty; **STOP** clinical adverse-risk prediction from cultured cells alone.
