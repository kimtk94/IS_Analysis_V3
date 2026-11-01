# P5q and P5s — formal 12-gene nine-pool-out robustness and independent numerical QA (2026-10-31)

## Verified real GitHub Actions evidence

- **P5q** [GitHub Actions SUCCESS #41170641379](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41170641379), [aggregate-only artifact #13031830547](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41170641379/artifacts/13031830547): original GSE313736 18 10x count matrices reconstructed for 57 paired BeamA donors. Nine source-pool omissions, with **both** 57-donor paired edgeR `~pair_id+treatment` (robust=FALSE sensitivity) and conservative source pooled-library `~pool+treatment` (robust=TRUE), for **18 formal QL refits**.
- **P5s** [independent artifact reconciliation SUCCESS #41172431153](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41172431153), [artifact #13031665988](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41172431153/artifacts/13031665988): separately recomputed **216/216** target/model/omitted-pool raw rows and **24/24** target/model summary rows; cross-checked Ensembl identities, reference logFC against P5r, exact threshold baseline, BH-FDR/sign counts and worst FDR.
- P5r earlier [3/3 parallel jobs SUCCESS #41172214095](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41172214095): see `P5R_THREE_PARALLEL_AUDITS_ACTUAL_PASS_20261031.md` and its mode-specific artifacts.
- P5q public files have only aggregate gene-level tabulations, source-library QC, and no individual donor trajectories or count matrices. Runner-private counts removed in the successful workflow. The published source GSM accessions in `source_qc.json` are public dataset identifiers, not donor IDs.

## Formal 9/9 result in BOTH edgeR models

**12/12 target genes tested and significant (BH-FDR<0.05, original-treatment sign stable) for each of the nine omitted source pools, in both model families. All 216 individual gene/model/omission entries were tested and FDR significant.** Each FDR is adjusted across the full filtered gene universe *in that refit*, not 12 tests alone.

| Gene | Mean logCPM, original primary | Worst FDR across nine conservative pooled-library fits | Pool-level descriptive signs |
|---|---:|---:|---:|
| ABHD3 | 3.515 | 0.00011731 | 9/9 |
| ADAMTS15 | 2.502 | 0.00025723 | 9/9 |
| DDIT4L | 5.302 | 0.00046348 | **8/9** |
| EGR3 | 4.288 | 0.00168641 | 9/9 |
| LPAR1 | 7.223 | 0.00009101 | 9/9 |
| RAB32 | 7.047 | 0.00045521 | **8/9** |
| SDC1 | 4.433 | 0.00023447 | 9/9 |
| STK17B | 6.910 | 0.00002883 | 9/9 |
| TFPI | 6.820 | 0.00137195 | 9/9 |
| CHAC1 (low expression) | 1.331 | 0.00049348 | 9/9 |
| ERFE (low expression) | 0.239 | 0.00010613 | 9/9 |
| PARM1 (low expression) | 1.004 | **0.01449488** | 9/9 |

Why P5o 8/9 pool signed counts for DDIT4L and RAB32 are **not a contradiction** of formal P5q 9/9 results: P5o describes each individual unadjusted source-pool TMM logCPM signed average; P5q estimates a treatment coefficient in eight pooled pairs after omitting a source pool, re-normalizing and refitting the quasi-likelihood model. Different estimands and normalization procedures must remain separate.

## Stress tests and publication-level significance

P5r post-hoc abundance/effect screen retains **9** at mean logCPM >=2 and all three contrasts absolute log2FC >=0.5; **8** at logCPM>=3, **7** at >=4, **5** at >=5. Increasing the three-contrast minimum effect-size requirement to 0.75 while holding logCPM>=2 retains **zero**. Thus this nine-gene shortlist is **threshold-sensitive** even though all 12 pass formal source-pool exclusion significance.

In official external 100nM bulk TM *published* timepoint contrasts, **EGR3** had 1h log2FC -0.357/FDR .00930 and 48h -0.825/FDR .00993, while **ERFE** had 6h +0.864/FDR .04026 and 48h +0.804/FDR .01242. No formal matched time-course interactions were performed. Across 12 candidates, **8** had at least one early raw sign different from 48h, generally with early FDR not significant; this is not evidence of transcriptional reversal.

## Research classification

- **GO**: original within-study 57 paired donor count analysis, 9 paired-source-library model and the targeted formal 18 refits; all 12 exhibit strong internal numerical robustness. This does **not** strengthen the degrees of independence or the causal interpretation.
- **HOLD**: original DEX vs Vehicle treatment was coupled to distinct 10x library preparation; donor blocking or LOO cannot identify pharmacological causality against batch effects. Post-treatment high-confidence BeamA cell selection may bias contrasts. External 100nM bulk TM and 500nM BeamA data differ in donor cohort, dose, cell definition, sample processing and technologies.
- **HOLD**: novel TM molecular mechanism and protein function, particularly 3 low-expression targets. Do not advertise original 96 discordances or these 12 as novel genes without separate literature mapping and orthogonal experiments.
- **STOP**: patient-level steroid-induced delta-IOP or glaucoma adverse-event prediction from these cultured TM datasets; no clinical outcome observed.

## Next decisive experiment — not executable in GitHub Actions alone
Within the same donor-derived TM culture, independently prepare and randomize Vehicle, 100nM DEX and 500nM DEX at 1/6/48 hours; balance extraction/library batches against all conditions; preregister dose×time contrasts and a primary endpoint; test gene expression with orthogonal qPCR and proteins, particularly low-abundance target assay detection limits. Only actual independent experimental evidence can move mechanism or clinical claims out of HOLD.

## Operational CI policy
P4 barcode preflight, P4 external metadata and P5c ENA metadata workflows were repeatedly invoked on every new ocular PR commit because **PR path filters are evaluated against the base branch**, even when the immediate commit modifies unrelated ocular files. Their earlier source/metadata checks were already executed. These three workflows are now **manual dispatch only**; `ocular-pr-smoke.yml` continues to run lightweight automatic PR validation. New expensive or redundant source transfer is NOT scheduled on each commit.
