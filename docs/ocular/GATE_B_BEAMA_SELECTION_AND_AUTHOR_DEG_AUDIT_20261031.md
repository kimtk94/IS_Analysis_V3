# Ocular Gate B — BeamA confidence selection bias and original author 80 DEG audit (2026-10-31)

## Gate B1: verified metadata-only sensitivity (PASS, not RNA model)

[GitHub Actions #41198576725](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41198576725) **SUCCESS**; artifact [#13041721835](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41198576725/artifacts/13041721835) downloaded and JSON contents directly inspected.

| GEO final SNG BeamA metadata | Confidence ≥0.75 | Confidence ≥0.50 |
|---|---:|---:|
| Donors with ≥10 cells per condition | 57 | 58 |
| Donors with ≥30 cells per condition | **57** | **57** |
| Donors with ≥50 cells per condition | **57** | **57** |
| Selected Vehicle metadata barcodes among ≥30-cell donors | 40,115 | 48,903 |
| Selected DEX metadata barcodes among ≥30-cell donors | 46,454 | 55,652 |

Under the frozen ≥30 cells per arm threshold, 0.50 retains exactly the same 57 paired donor identities as score 0.75 (set inclusion tested, zero donor additions). A .50 cutoff admits **8,788** additional vehicle and **9,198** additional DEX cell-metadata rows; these are **metadata barcode counts**, not yet matrix-joined cell totals, and are NOT independent biological replicates. The pooled fraction shift does not establish within-donor DEX causal effects.

## Gate B2: source author's 63 up/17 down high-effect gene identities (HOLD)

The [official GSE313736 GEO Series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313736) text reports **63 upregulated / 17 downregulated** high-effect response genes, 74 original donors and 316,799 singlets. Its reported downregulated cutoff `log2 fold change <3` is ambiguously worded and should **not** be assumed equivalent to `log2FC <= -3`.

Actual public GEO Series supplementary inventory contained cell metadata, classification, donor metadata, VCF genotypes and RAW.tar, but **no explicit original author 80-row DEG gene table**. No exact Ensembl-ID set overlap of the original authors' 80 with our 0.75 **83** (67 up/16 down) can be calculated yet. The 83-vs-80 count difference can stem from donor attrition, cell-state filtering, testing universe, fold-change thresholds, and statistical model, not automatically error or novelty.

**Decision:** GO numerical published abstract-vs-our-count comparison, **HOLD original 80-gene identity-level overlap** until the author's actual supplementary gene IDs are obtained. Do not infer an 80-member set from count summary.

## Gate B3: complete raw UMI 0.50 vs 0.75 paired edgeR re-fit (execution underway)

A dedicated on-demand GitHub Actions workflow `.github/workflows/ocular-gate-b-050-raw-edger.yml`:
1. Re-download original 18 GSE313736 raw 10x integer matrices in transient GitHub runner memory; validate barcodes, feature-order, gzip CRC and MatrixMarket NNZ at **score≥0.50**, ≥30 cells per arm.
2. Fit same original `~ pair_id + treatment` edgeR TMM robust NB QL with correctly documented score0.50 in `scripts/ocular/ocular_gate_b_050_edger.R`. Biological unit is donor, not cells.
3. Compare gene-level score0.50 effect/FDR against frozen 0.75 result from [Actions #41085982185](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41085982185); assess MYOC/TGFBI/MMP2, ACTA2, 16-member panel, all-common gene correlations and exact internal high-effect set changes. These are internal reanalysis sensitivity tests, **not independent replication**.
4. Delete runner-private donor-level expression/pair files in an `if:always()` step; upload only aggregate gene-level statistics and anonymized graphics.

**Initial Gate B3 [Actions #41198827259](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41198827259) FAILED** during GEO source download with `BAD_GZIP_MAGIC` before any edgeR model. It was a transport/content-type issue, not proven corruption of 10x count data. We implemented limited retries and an official NCBI GEO download alternative, preserving strict gzip validation. **Retry [Actions #41198957708](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41198957708)** was queued on launch. Its actual biological/statistical conclusion must be appended from a verified artifact; this document does NOT preemptively mark B3 PASS.

## Next paper-facing decisions

- **GO candidate robustness:** only if score0.50 full raw-count re-fit agrees in sign and adjusted inference for frozen three-study candidates **MYOC, TGFBI, MMP2**, without meaningful dependence on classifier cutoff.
- **HOLD technical/drug causality:** each source pool still uses condition-specific 10x library, and BeamA classification is post-DEX. Threshold stability cannot identify condition-specific library effects or remove selection/collider bias; test alternate broad cell-state annotation with a feasible n, or explicitly limit claims.
- **STOP clinical steroid-induced IOP prediction:** GSE313736 cultured cells have no directly linked steroid-exposure ∆IOP outcome. Microbiome 16S/shotgun stays separate until clinical sample mapping.
