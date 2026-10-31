# P5q — formal leave-one-source-pool-out edgeR inference for 12 target TM genes

Date: 2026-10-31. **Status: IMPLEMENTED / NOT YET EXECUTED.** Do not label PASS until the actual workflow result, artifact, and per-gene table are inspected.

## Why P5q is necessary

P5b performed nine formal donor-paired leave-one-source-pool-out edgeR QL refits but published only the top100 and selected canonical genes; **0/96** of the P5j opposite-sign candidates appeared in that export. P5o reconstructed all 18 original 10x source matrices and found descriptive sign stability for the 96; **pool-level mean logCPM deltas and bootstrap intervals are not edgeR QL inference**.

The P5p post-hoc expression-aware list defines exactly 12 targeted hypotheses, which P5q holds fixed using P5l 96-gene Ensembl IDs rather than symbol-only joins:

- High effect, low expression: **ERFE, CHAC1, PARM1** (average logCPM <2).
- Adequate expression: **ADAMTS15, EGR3, LPAR1, ABHD3, DDIT4L, TFPI, STK17B, RAB32, SDC1** (average logCPM >=2).

The two external study conditions must remain distinct: our selected BeamA state, nominal **500 nM DEX, 48 h** versus published bulk TM, **100 nM DEX, 48 h** (and 1/6 h supplementary measures). Sign discordance is not proof of a drug-dose interaction.

## Execution implementation

Manual-only Actions workflow: `.github/workflows/ocular-p5q-12-target-formal-loo.yml`.

Runner: `scripts/ocular/ocular_p5q_12_target_formal_loo.R`.

1. Download verified aggregate P5l source [Actions #41158660082](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41158660082), with a **96 unique Ensembl ID / exact 12-target** schema guard.
2. Rebuild **57 paired donor pseudobulks** (114 columns; 36,601 integer genes; 1,927,677,307 UMIs) from all **18 original matrices**. Keep donor-level metadata/counts in transient `RUNNER_TEMP` only.
3. For each excluded source pool A–I, refit **both** the paired donor design `~ pair_id + treatment` (robust=FALSE, matching prior P5b bounded sensitivity method) and **pooled 10x-library unit design** `~ pool + treatment` (robust=TRUE, 8 source-pool pairs / 16 condition libraries). TMM normalization, treatment-group `filterByExpr`, estimated dispersion and quasi-likelihood testing are recomputed for **every exclusion**.
4. Apply BH-FDR **over all genes retained in each model**, not just the 12 post-hoc hypotheses. Report target filtering dropout explicitly, never silently impute significance.
5. Emit only 216 target/model/exclusion gene-level rows and 24 target/model summary rows, plus a non-sensitive source QC JSON. Never export sample IDs, individual count matrices or donor trajectories.

Public expected artifact name: `ocular-p5q-12-target-formal-loo` (retained 30 days).

## Predefined numerical readouts, without unearned GO labels

For each of 12 genes *in each model*, report tested count (of 9), original-contrast sign agreement count (of 9), FDR<0.05 count, **FDR<0.05 AND correct sign** count, worst FDR and smallest absolute refit log2FC. A filtering dropout is automatically not a 9/9 PASS.

- Strongest **within-study numerical robustness** criterion: 9/9 tested, 9/9 signed concordance, and 9/9 model-specific BH FDR<0.05 **in the conservative pooled-library model**. Report donor-model results separately, never as 57 independent 10x libraries.
- 8/9 significant/sign-stable: **sensitivity warning**, not a pooled 9/9 pass.
- 9/9 same sign but only 0–7/9 FDR: effect direction may be stable while inferential stability fails.
- 0–8/9 tested: **filtering-dependent / HOLD**, do not interpret absent results as nonsignificance.
- An executed numerical PASS means successful 18 model fits and transparent readouts, **not** that all 12 genes pass biological reproducibility.

## Research-level gates

| Question | Verdict now | Evidence required to change it |
|---|---|---|
| Source-count reconstruction, 96 discordances, 3+9 stratification | GO / internally verified | Existing P5l/m/o/p reports and artifacts |
| Formal 12-gene 9x+9x exclusion robustness | **PENDING** | P5q actual Actions run and per-gene artifact |
| Biological DEX dose-dependent reversal | **HOLD** | Independently prepared matched cultures, same donor and cell-state assay, Vehicle/100nM/500nM × 1h/6h/48h, randomized/balanced culture and library batches, donor-blocked dose×time interaction |
| Mechanism / protein translation | **HOLD** | Orthogonal qPCR or targeted proteomics and perturbation, cell-state markers, sufficient low-expression detection |
| Clinical steroid IOP risk prediction | **STOP** | Independent, prospective patient-level treatment and ΔIOP measures, external validation and clinical model calibration |

### Avoid interpreting post-hoc selection as independent replication

All 12 were picked because their original and external study directions opposed; leave-one-pool-out inference reuses the original dataset. Multiple comparisons and subgroup analyses are already conditioned on prior discovery. Neither stable source-pool effects nor pooled-library QL significance resolves the structural confounding: condition-specific 10x library preparation was not randomized independently of DEX. Do not claim a drug concentration switch, clinical causal mechanism, or nine novel genes.
