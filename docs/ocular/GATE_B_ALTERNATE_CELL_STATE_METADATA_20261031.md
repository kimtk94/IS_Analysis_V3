# Gate B — alternate cell-state annotation donor feasibility verified (2026-10-31)

## Verified execution

- [GitHub Actions #41201790738](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41201790738): **SUCCESS**, [artifact #13043246012](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41201790738/artifacts/13043246012) downloaded, unzipped and inspected (public aggregate JSON + TSV).
- Initial run [#41201695447](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41201695447) failed early due truncated GSE313736 public metadata gzip. Reused previously successful bounded official GEO retry + full CRC code and reran. Failure history retained.
- This is **metadata feasibility only**, not a re-fit of donor-level expression and not a claim of biological independence. No donor-level barcodes, genotypes or phenotype data were committed/uploaded.
- Exact script \`scripts/ocular/ocular_gate_b_alt_annotations.py\`; reused on-demand \`.github/workflows/ocular-gate-b-eligibility.yml\`. Source official GSE313736 cell metadata 169,398 rows.

## Counts at ≥30 annotated cells in each paired condition

| Post-treatment cell-state definition | Eligible donor pairs | Overlap with primary BeamA score≥0.75 donor set | Selected metadata Vehicle/DEX cells | Number of source pools |
|---|---:|---:|---:|---:|
| CellTypist BeamA ANY score | 59 | 57 | 67,047 / 74,932 | 9 |
| **BeamA score≥0.75** primary | **57** | 57 | **40,115 / 46,454** | 9 |
| **BeamA score≥0.50** verified paired edgeR | **57** | 57 | **48,903 / 55,652** | 9 |
| **scpred_vanzyl2022 = TM_Fibro** | **57** | **55** | **13,233 / 19,259** | 9 |
| BeamA score≥0.50 ∩ scPred TM_Fibro | **44** | 44 | 6,550 / 9,380 | 9 |
| BeamA score≥0.75 ∩ scPred TM_Fibro | **37** | 37 | 4,781 / 6,851 | 9 |
| CellTypist JCT ANY | **56** | 55 | 8,186 / 8,611 | 9 |
| CellTypist JCT score≥0.75 | **8** | 8 | 481 / 564 | **6** |

Cell counts above are **metadata barcodes**, not joined raw 10x-filtered matrix cells, and are not biological sample sizes. Stricter scPred×BeamA intersections materially reduce eligible donors and cell depth. JCT score≥0.75 **n=8** is not suitably powered for the 57-donor model, do not force the same statistical method/conclusion.

## Decision on next raw-count analysis

**Prioritize scPred TM_Fibro** rather than attempting every annotation and threshold: n=57 eligible paired donors, 9 source pools, different cell-label system, and 55 donors overlap the frozen primary BeamA score≥0.75 set. This tests dependence of MYOC/TGFBI/MMP2 on the BeamA classification pipeline.

Important caveats:
- 2 scPred-eligible donors were NOT in the original BeamA57; strict within-donor annotation comparison must use the **55 overlapping donors** on both sides (or explicitly separate donor change from cell annotation). A simple scPred 57 vs BeamA 57 differential comparison is *partially confounded by donor composition*, though useful as descriptive robustness.
- All labels are inferred from cells **after** DEX; this cannot remove collider/selection bias completely. Treatment-specific 10x library preparation remains perfectly condition-aliased, so alternatives do **not** establish causality.
- scPred TM_Fibro has less per-arm selected cell depth, and barcode intersection may drop some donor pairs below 30. The metadata n=57 does not guarantee 57 surviving the actual matrix join; the existing P4 code requires at least 50 but cannot assume exact 57 without evidence.
- Before costly 18-library reruns, produce one single-pass multi-cell-state aggregation design to avoid repeated ~6GB GEO fetches and track donor sets, controls and gene-universe; only execute it with guaranteed separation of donor eligibility and a fixed common donor set.

**GO:** full primary raw UMI score 0.50 vs 0.75 sensitivity (verified [#41199532259](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41199532259)), and alternative cell-annotation feasibility (this task).

**HOLD:** actual scPred TM_Fibro gene-level edgeR and fixed 55-pair direct comparison; author's original 63-up/17-down individual gene list; causal treatment vs batch; functional protein/IOP effects.

**STOP:** treating n=8 high-confidence JCT as a powered 57-pair replication, or treating any label sensitivity as a patient steroid-response classifier.
