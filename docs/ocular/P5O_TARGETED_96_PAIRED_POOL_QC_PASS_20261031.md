# Ocular P5o — 96-gene paired-donor and 9-pool descriptive sensitivity VERIFIED (2026-10-31)

## Evidence and scope
- **[GitHub Actions run #41160840453](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41160840453) SUCCESS**, [artifact ID 13027396135](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41160840453/artifacts/13027396135).
- 18 published GSE313736 source 10x libraries streamed from NCBI, gene ordering and raw UMI checked; 57 high-confidence BeamA **biologically paired donors** (114 pseudobulks), 36,601 raw gene features, **1,927,677,307 raw UMIs**.
- Actual runner-private files held gene-level counts and pseudonymous donor mappings only during analysis and were deleted before artifact upload; public files `P5o_targeted_pool_sensitivity_summary.json`, `P5o_96_targeted_paired_pool_sensitivity.tsv`, `P5o_three_targets_pool_mean_directions.png`, `source_qc.json` provide only aggregate evidence.
- Code `scripts/ocular/ocular_p5o_targeted_paired_pool.R`, action `.github/workflows/ocular-p5o-targeted-pool.yml` (**manual dispatch after first run**).
- **Important: not 9 edgeR QL refits**. TMM normalization of full count matrix, per-donor delta log2CPM, nine equal-weight source-pool mean deltas, nine leave-one-pool-out pooled-mean effects and 10,000 seeded bootstrap resamplings of the nine pool means. No new FDR computed. This is an *orthogonal descriptive sensitivity* complementing a separate, already-validated 9 paired-library-unit edgeR model (P5d).

## Actual full 96 results
| Criterion | Observed |
|---|---:|
| External-publication FDR-significant but sign-opposite genes reanalyzed | **96/96** |
| Both original 500nM and external 100nM 48h published DEG | 96 by prior P5j/P5l |
| Same signed effect in **all 9 original donor pools** | **33/96** |
| Same signed equal-pool mean effect after **each of the 9 leave-one-pool-out exclusions** | **96/96** |
| Nine-pool bootstrap 95% interval excludes zero | **75/96** |
| Formal 9-fold robust edgeR QL re-fit of these 96 genes | **NOT PERFORMED** (original P5b LOO only covered top100, zero overlap) |

These results strengthen signed within-study robustness, especially for strong candidates, but the source 9 pool pairs were each processed as DEX and Vehicle in distinct 10x libraries. An across-pool shared preparation artifact could also produce the consistent signed effects. Bootstrap intervals are **conditional descriptive estimates across the nine source pool units**, not independent causal experiment confidence intervals.

## Largest-effect low-abundance discordant Tier-A candidates
| Gene | Primary donor-paired edgeR log2FC | Primary mean logCPM | 57 donor-pair delta same sign | All nine pools same sign | All nine omit-one-pool means same sign | Pool-mean delta 95% bootstrap interval |
|---|---:|---:|---:|---:|---:|---|
| **ERFE** | -1.207 | **0.239** | **50/57** | **9/9** | **9/9** | **[-1.299, -0.773]** |
| **CHAC1** | -1.147 | **1.331** | **50/57** | **9/9** | **9/9** | **[-1.474, -0.718]** |
| **PARM1** | +1.143 | **1.004** | **48/57** | **9/9** | **9/9** | **[+0.643, +1.762]** |

Within the published 100nM TM bulk time-course (IOVS 2026 supplementary workbook), ERFE has significant **positive** log2FC at 6h (+0.864, FDR .0403) and 48h (+0.804, FDR .0124), opposite our negative BeamA effect at nominal 500nM 48h. CHAC1 and PARM1 external 100nM 48h effects have the opposite sign, but their 1h/6h timepoint effects are nonsignificant. These are dose/cell-state/technical **hypotheses**, **not proven concentration-dependent reversal**. All three have relatively low primary mean expression, raising orthogonal detection concerns.

## Nine adequately expressed opposite-sign candidates
P5p pre-filter defined *post hoc*: primary logCPM>=2, primary abs(logFC)>=.5, external 48h abs(logFC)>=.5, conservative nine-library-pair edgeR BH-FDR<.05, abs(logFC)>=.5 and same as primary sign.
| Candidate | Sign-concordant among nine source pools | Nine leave-one-pool pooled means stable | Pool bootstrap interval excludes 0 |
|---|---:|---:|---|
| **DDIT4L** | 8/9 | 9/9 | YES |
| **ADAMTS15** | 9/9 | 9/9 | YES |
| **LPAR1** | 9/9 | 9/9 | YES |
| **STK17B** | 9/9 | 9/9 | YES |
| **EGR3** | 9/9 | 9/9 | YES |
| **ABHD3** | 9/9 | 9/9 | YES |
| **RAB32** | 8/9 | 9/9 | YES |
| **TFPI** | 9/9 | 9/9 | YES |
| **SDC1** | 9/9 | 9/9 | YES |

**7/9** consistent across every source pool; 2/9 across 8 pools, but all nine maintain omitted-pool mean effect sign and descriptive 95% pool bootstrap interval excluding zero. Do not call them nine novel genes or independently replicated dose-dependent phenomena. Existing IOVS gene-level result clearly implies they were already tested and FDR significant in an external human TM publication but with opposite direction.

## Research GO/HOLD/STOP
- **GO** source raw-count integrity, donor pairing, 9-pool source-libraries, 96 independent published TM gene-ID comparisons, and a new validated *descriptive within-source pool robust sign screen*.
- **HOLD** formal 9x edgeR leave-one-pool QL refit for these 96, post-treatment BeamA selection bias, treatment-vs-library preparation artifacts, exact original GSE313736 manuscript's 63 up/17 down full gene list, mechanistic novelty and independently randomized 100/500nM validation.
- **STOP** clinical steroid-induced IOP/POAG risk prediction using only cultured TM cells lacking observed clinical donor-level ΔIOP.
- Pre-registered next intervention hypothesis design should use donor-blocked *same-culture* Vehicle/100/500nM DEX at preplanned 1h/6h/48h, randomize extraction and 10x library prep independently of dose, and measure qPCR/protein + ER stress/ECM readouts if the biological mechanism is supported; no inferential claims before such a study.
