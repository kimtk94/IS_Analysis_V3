# P5b pre-validation: primary gene-level result diagnostics (2026-10-31)

## Evidence
This table is a **secondary re-aggregation of real P5 primary results**, not a new edgeR fit. Primary run [41074586446](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41074586446), artifact [12987425927](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41074586446/artifacts/12987425927), file `p5_all_gene_results.tsv.gz` (16,050 genes).

## Threshold stratification
| Filter | Genes |
|---|---:|
| Total genes tested | 16,050 |
| BH FDR < 0.05 | **13,717** (85.46%) |
| BH FDR < 0.05 and abs(log2FC) >= 0.25 | **10,070** |
| BH FDR < 0.05 and abs(log2FC) >= 0.5 | 5,340 |
| BH FDR < 0.05 and abs(log2FC) >= 1.0 | **1,365** |
| BH FDR < 0.05 and abs(log2FC) >= 2.0 | 246 |
| BH FDR < 0.05 and abs(log2FC) >= 3.0 | 83 |
| BH FDR < 0.05 but abs(log2FC) < 0.25 | **3,647** |

The tested-gene median absolute log2FC is 0.3488; median signed log2FC is 0.0342. Extreme statistical significance does not imply a very large gene-level effect; sample size, technical confounding and regularity across donors can all yield small-P phenomena. FDR is not an effect-size filter.

## Discovery rate by average log-CPM (among genes surviving primary filter)
| logCPM stratum | n | Share FDR<0.05 |
|---|---:|---:|
| <0 | 948 | 71.62% |
| 0–2 | 3,454 | 78.98% |
| 2–4 | 4,219 | 84.47% |
| 4–6 | 4,689 | 90.34% |
| 6–8 | 2,187 | 91.40% |
| >=8 | 553 | 92.41% |

Even high-abundance genes show broad significance, so low-abundance stochastic noise alone cannot explain the phenomenon. **This does not rule out DEX response, general library confounding or improper normalization.**

## Unresolved before claim of manuscript-grade biology
- DEX versus vehicle were sequenced in separate 10x libraries; donor blocking protects donor baseline but does *not* separate treatment and library preparation artifacts.
- Nine formal leave-one-pool-out models are running in [Actions 41085982185](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41085982185); values should only be reported when the artifact exists and all 9 fits pass. Robust primary edgeR and non-robust sensitivity fits use different dispersion estimators; compare signed effects first, FDR second.
- Corrected anonymous MDS point plot is part of the same workflow; validate actual rendered PNG, not just an exit code.
- PRJNA1243143 external 46-run metadata cannot identify treatment/time; separate formal independent replication remains HOLD [Actions 41086796126](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41086796126).
- HIPK2, PTK2B, FKBP5, ZBTB16, ATP10A and EDIL3 already occur in preceding DEX experiments, thus their expression direction is **not a novel gene discovery**; see `P5B_LITERATURE_GENE_CROSSWALK_20261031.md`.
- No clinical steroid IOP phenotype available in this dataset, therefore no clinical adverse-effect prediction claim.

**Decision:** primary DE pipeline executed; manuscript-level robustness and novelty remain HOLD until pool sensitivity and independent evidence.
