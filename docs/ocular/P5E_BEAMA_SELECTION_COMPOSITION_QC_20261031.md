# Ocular P5e: actual BeamA state-selection / compositional shift audit — 2026-10-31

## Execution provenance
- GitHub Actions [run 41091146021](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41091146021) PASS, [artifact 12994476349](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41091146021/artifacts/12994476349), executed `scripts/ocular/ocular_p5e_cellstate_bias.py`.
- Public dataset GSE313736, **169,398 SNG cells**, 59 donors with data in both DEX and no-DEX arms. Prespecified annotation `celltypist_chen2025_fine == BeamA`, score >=0.75.
- **No individual donor IDs or counts exported**: only global and 9-pool aggregate metrics, paired-donor median and exact binomial sign test.

## Observed global-versus-donor discrepancy
| Metric | Vehicle | DEX |
|---|---:|---:|
| SNG cells | **80,495** | **88,903** |
| BeamA any classification | 67,047 | 74,932 |
| BeamA score>=0.75 | **40,139** | **46,471** |
| BeamA >=0.75 fraction of all SNG | **49.8652%** | **52.2716%** |
| BeamA >=0.75 among all predicted BeamA | **59.8670%** | **62.0176%** |
| BeamA classification score median (within predicted BeamA) | **0.86667** | **0.89185** |
| scpred independently assigned TM_Fibro (any score) | **13,274** | **19,428** |
| Both BeamA (any) and scpred TM_Fibro | 10,483 | 15,473 |

Pooled fraction DEX minus Vehicle is **+2.4064 percentage points**. However, across **59 paired donors** the median per-donor BeamA score>=.75 fraction difference is **−0.421 percentage points**, with **29 positive**, **30 negative**, exact two-sided sign-test **p=1.000**. For BeamA without confidence threshold: 32 positive/27 negative, p=0.603. These contrasts indicate pooled cell-count weights, capture and/or pool composition explain why an apparent global shift may not represent the median donor.

This is **not** proof that glucocorticoid causes zero change: the sign test addresses only the direction balance; small effects, technical capture differences and gene-expression-driven label assignment remain possible. Sign-test nonsignificance does not establish equivalence.

## Pool-level heterogeneity (BeamA score>=.75 / all SNG)
| Pool | Vehicle % | DEX % |
|---|---:|---:|
| A | 60.10 | 62.04 |
| B | 61.32 | 64.88 |
| C | 47.64 | 47.98 |
| D | 57.43 | 55.94 |
| E | 49.56 | 51.80 |
| F | 54.29 | 54.35 |
| G | 50.53 | 50.08 |
| H | 29.10 | 34.35 |
| I | 26.95 | 21.11 |

Pools H and I have markedly lower BeamA high-confidence fraction than A/B. **Cell-state composition is heterogeneous between donor pools**. This can induce selection/collider bias when selecting post-treatment cells by classifier score; compare gene effects under alternative score thresholds and broader class definitions.

## Manuscript decision
- **PASS diagnostic**: class-label and score composition QC based on real published cell metadata.
- **HOLD cell-state causal claim**: the 2.41-pp global fraction difference does not appear as a shared direction across paired donors. Do not present global cell-level counts as independent biological replicates.
- **HOLD broad pathway inference**: primary 57-donor edgeR DEG signal (13,717 FDR<0.05 of 16,050) may reflect actual strong DEX response, treatment-specific library handling, and/or post-treatment BeamA classifier selection. Need adjusted stability analyses, 9-pool library-unit QL test and external replication.
- **STOP clinical steroid-induced IOP prediction** from these cultured cells: no clinical outcomes.
