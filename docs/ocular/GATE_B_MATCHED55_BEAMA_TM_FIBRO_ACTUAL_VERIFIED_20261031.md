# Gate B — EXACT 55 same-donor BeamA versus scPred TM_Fibro raw pseudobulk validation (2026-10-31)

## Immutable actual CI evidence

This is **an internal cell-annotation sensitivity analysis**, not independent biological replication. Both cell annotation approaches use the **same 55 biological donor pairs**, all 18 original source 10x libraries and partially overlapping post-exposure cells. Treatment is condition-aliased with separate 10x libraries, so this analysis does not identify a causal DEX effect free from preparation/selection bias.

- **Full 18-matrix one-pass dual-label aggregation and two donor-blocked edgeR fits**: GitHub [Actions #41206289766](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41206289766). The overall workflow is **red/FAIL**, **solely because its final comparator required all 16 originally preselected genes to be testable in both cell-label subsets**. The actual raw-UMI aggregation **PASS** and **both edgeR model steps PASS**. [Full public aggregate artifact #13044762180](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41206289766/artifacts/13044762180) downloaded, decompressed and inspected. No raw donor count matrices were uploaded.
- A prior first run [#41206177550](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41206177550) failed immediately on a Python metadata counter type error, which was corrected. No invented biological results from that run.
- **Fixed lightweight aggregate-only reviewer**: [Actions #41207583296](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41207583296) **SUCCESS**, [review artifact #13045350639](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41207583296/artifacts/13045350639) actually downloaded, JSON and 16-row TSV inspected. Historical initial reviewer [#41207507021](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41207507021) failed on a code variable typo `NameError`, corrected without re-downloading any original matrix.
- Reproducible scripts:
  - `scripts/ocular/ocular_gate_b_55_dual_counts.py`: one-pass sparse coordinate decoding, accumulate both cell-state pseudobulks simultaneously, complete source gzip/NNZ QC, 18 library source hashes, strict joined 55 common donor check.
  - `scripts/ocular/ocular_gate_b_55_edger.R`: **two independent fits** to the same 55 paired donors, identical `~ pair_id + treatment` TMM and robust edgeR quasi-likelihood model; different *post-treatment cell annotation labels*, not independent donor data.
  - `scripts/ocular/ocular_gate_b_55_compare.py`: compare only shared gene ID test universe, report filtered-out genes explicitly, never synthesize untested FDR.
  - Expensive GitHub workflow `.github/workflows/ocular-gate-b-dual55.yml` and lightweight reviewer `.github/workflows/ocular-gate-b-dual55-review.yml` both manual `workflow_dispatch` after one-shot PR launches.
- **Privacy and sample unit**: donor IDs and 36,601 × 110 per-annotation raw count matrices existed only in runner temporary storage and were deleted in the `if:always()` step. Only aggregate gene statistics, anonymous MDS/MA figures, pool QC, source SHA256 and pseudonymous cell count summaries are in the artifacts.

## Actual source and modelling counts

| QC | BeamA score≥0.75, fixed 55 | scPred TM_Fibro, same 55 |
|---|---:|---:|
| Matrix-joined Vehicle cells | **38,924** | **12,876** |
| Matrix-joined DEX cells | **44,439** | **18,891** |
| Raw UMI counts Vehicle | **990,968,235** | **361,740,480** |
| Raw UMI counts DEX | **890,246,138** | **471,374,852** |
| Biological paired donor count | **55** | **55** |
| edgeR model source features | 36,601 | 36,601 |
| After `filterByExpr` testable genes | **16,051** | **14,003** |
| BH FDR<0.05 genes | **13,665** | **11,021** |

All 18 matrices passed full gzip CRC, feature-order and sparse UMI conservation checks. 55/55 matched donors retained ≥30 *actual joined* cells in **every donor×condition×annotation combination**. **These condition-dependent cell totals must not be presented as biological independent replicates or direct cell abundance effects.** Technical/post-exposure label selection and library batch remain possible confounding.

## Matched-universe gene concordance

| Metric | Verified count |
|---|---:|
| Shared Ensembl IDs tested in both models | **13,975** |
| Shared gene signed log2FC concordance | **12,634/13,975 (90.40%)** |
| FDR<0.05 in both | **10,255** |
| Shared gene log2FC Pearson correlation | **0.94990493** |
| Prespecified P9 candidate genes testable in both | **15/16** |
| Among 15 testable, same sign | **15/15** |
| Among 15 testable, FDR<0.05 in both | **15/15** |
| One prespecified candidate excluded by scPred low-expression filter | **ERFE**; no scPred effect/FDR estimate |
| MYOC/TGFBI/MMP2 primary 3 targets | **3/3 same sign and BH-FDR<0.05** |

### Gene-level validation table

| Gene | BeamA-55 log2FC | scPred-55 log2FC | Both FDR<0.05? | Evidence role |
|---|---:|---:|---|---|
| **MYOC** | **+4.00614** | **+3.52539** | yes | GC-responsive assay reference |
| **TGFBI** | **−1.86829** | **−1.94124** | yes | ECM transcript candidate |
| **MMP2** | **−0.59071** | **−0.45883** | yes | ECM-turnover transcript candidate, NOT measured enzyme activity |
| ACTA2 | +0.48093 | +0.61844 | yes | Annotation-stable **within this same cohort**, but opposite in published 7-day TM study |
| SDC1 | −0.55126 | −0.68400 | yes | Annotation-stable within-cohort, but published independent studies disagree |
| EGR3 | +0.75324 | +0.51720 | yes | Context-dependent across external studies |
| LPAR1 | +0.73995 | +0.35458 | yes | Context-dependent across external studies |
| TFPI | +0.61352 | +0.65307 | yes | Context-dependent across external studies |
| LUM | −2.48814 | −2.35374 | yes | Post-hoc ECM gene priority |
| IGFBP3 | −1.86537 | −2.11013 | yes | Post-hoc ECM gene priority |
| VCAM1 | −1.71051 | −2.24726 | yes | Post-hoc ECM gene priority |
| NTM | −2.54239 | −3.46724 | yes | Post-hoc ECM gene priority |
| FKBP5 | +3.97672 | +3.50190 | yes | GC assay control |
| CHAC1 | −1.11301 | −0.49115 | yes | Low-abundance exploratory |
| **ERFE** | −1.21625 | **NOT TESTED** | **not assessable** | **scPred expression filtered**, not a negative FDR |
| PARM1 | +1.15896 | +0.98632 | yes | Low-abundance exploratory |

**Do not collapse the evidence tiers**: the 55+55 annotation models share biological donors, source libraries and many cells, so the exceptionally high sign concordance and large number of double FDR-significant genes are expected to be statistically dependent. The right inference is *robustness to one change of annotation method*, **not** independently replicated transcriptomic response or DEX pharmacological causality.

## Study decision and next specific experiments

- **GO computational robustness**: original BeamA 0.50/0.75 raw UMI threshold robustness passed previously (57 donors, r=0.9991); current *matched-55 alternative annotation* also passes and core MYOC/TGFBI/MMP2 are signed and FDR stable.
- **HOLD complete cellular-bias mechanism**: selection is still post-exposure and the two annotation labels overlap, so treatment-induced cell-state reclassification versus transcriptional change is not identifiable from these comparisons alone. Nonoverlapping cell subsets, negative-control transcripts, and sample/library QC may help characterize the extent of selection effects, but cannot undo perfect condition/library aliasing.
- **HOLD wet-lab mechanism**: no independent randomized matched-donor dose/time design, secreted MYOC, MMP2 gelatinase activity, ECM imaging/stiffness or aqueous-outflow effects.
- **STOP** single-gene steroid glaucoma predictor or mediation through unmatched ocular-surface microbiome, since patient steroid ΔIOP outcome is not linked.
- **Next research-effort priority**: publication-ready figure and statistics comparing original GSE313736 57-donor BeamA, 55-donor BeamA vs 55-donor TM_Fibro, 100nM 48h Mehrotra and 500nM 7-day GSE301517, with explicit effect scales, paired units and causal limitations; original GSE313736 author 63↑/17↓ exact high-effect gene-ID list remains inaccessible in current public Series inventory. New raw-data re-download is not justified without a pre-specified distinct hypothesis.

**Gate B overall**: computational cell-state-selection sensitivity **GO**; mechanistic steroid causality, protein function and clinical prediction **HOLD/STOP**.
