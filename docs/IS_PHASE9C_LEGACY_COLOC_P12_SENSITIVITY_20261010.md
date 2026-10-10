# IS P0 — Legacy ABF p12 sensitivity and donor-cell evidence integration

**Analysis date:** 2026-10-10

**Status:** COMPUTATIONAL_COMPLETE_CONDITIONAL / MANUSCRIPT_CAUSAL_CLAIM_BLOCKED

## Data used, without expanding scope retrospectively

- Legacy 646 GTEx v8 gene–tissue ABF PASS tests, limited to four original BBJ IS loci.
- Historical Bayesian hypotheses H0–H4, including 43 genes and 16 tissue datasets. Original outputs were not independently recomputed from SNP-level beta/SE files.
- R3_7 GSE256493 adult healthy temporal-lobe human reference, 80,515 cells across six `Patient` labels, 10 predeclared donor-paired contrasts.
- Broad ancestry-aware discovery remains 80 **provisional distance windows**, 2,425 positional gene–region assignments, 2,225 unique gene IDs, with one separate NEURL1 legacy-only anchor row; **the 646 GTEx tests cannot represent this larger universe**.

## Prior sensitivity design

- Software: installed coloc **5.2.3**. Run the package's own internal `prior.adjust` against each complete H0–H4 ABF summary (not manually multiplying H4 alone).
- Conditional baseline assumes `p1=1e-4`, `p2=1e-4`, `p12=1e-5`. **The original legacy coloc run script and prior provenance have not been independently established**, so these are assumed settings, not audited historical metadata.
- Grid: `1e-6, 3e-6, 1e-5, 3e-5, 1e-4`, holding p1/p2 constant.
- Baseline `p12=1e-5` reproduced all 646 original five-hypothesis posterior vectors within 1e-9; 3,230 conditioned rows generated.
- The focal comparison chooses the tissue with **highest baseline H4 per gene**; this is post hoc exploratory tissue selection, not multiple-testing-adjusted inference.
- Inference is only conditional prior sensitivity of historical **single-signal ABF summaries**. It does not substitute for variant-level coloc, SuSiE, prior-source verification, allele harmonization, multi-signal LD, molecular QTL power or independent replication.

## Selected baseline-best tissue PP.H4

| Gene | p12=1e-6 | 1e-5 | 1e-4 | Interpretative status |
|---|---:|---:|---:|---|
| FGF5 | 0.260 | 0.779 | 0.972 | Sensitive; stronger evidence requires authentic prior and SNP/LD provenance |
| CALHM2 | 0.283 | 0.798 | 0.975 | Sensitive; different locus from FGF5, not directly ranked by H4 |
| SH3PXD2A | 0.051 | 0.351 | 0.844 | Strong prior dependency |
| COL4A2 | 0.047 | 0.331 | 0.832 | Strong prior dependency |
| ALDH2 | 0.009 | 0.079 | 0.462 | Coding variant branch, do not force bulk eQTL model |

The expanded gene universe must not be narrowed to these five based on old ABF posteriors. Even near-one PP.H4 under a generous prior does not establish a causal gene.

## Integration with Phase9F R3_7 donor evidence

The independent six-donor pseudobulk re-audit found 9/10 comparisons with unanimous directional CPM differences and 8/10 with unanimous detection-fraction differences. SH3PXD2A showed descriptively higher expression in Oligodendrocytes/Fibroblasts/Smooth muscle compared with the **combined** Microglia and Macrophages annotation; COL4A1/2 had higher mural/fibroblast CPM than endothelial cells in eligible donor pairs. This does not establish allele-specific regulation.

- FGF5 is **NOT_ASSESSABLE_IN_THIS_REFERENCE** due to no matched Seurat RNA feature, not a biological zero.
- COL4A2 pericyte versus endothelial CPM is 6/6 positive, but detection-rate comparison is only 3/6 positive.
- ALDH2 endothelial versus immune class is 3/6 CPM positive; no uniform donor effect.
- Neurons 15,903/15,942 from donor TL8 — do not infer six-donor neuronal replication.
- All human reference data are from **normal temporal lobe**, not disease cases, and cannot validate stroke-specific transcriptional effects.

Integrated noncausal ledger:
`results/is/audits/is_p0_integrated_genetics_cell_20261010_v1/`

Conditional ABF p12 reweighting:
`results/is/audits/legacy_coloc_p12_sensitivity_20261010_v1/`

Source frozen prior-conditioned analysis and 9-gene evidence ledger outputs are versioned and separate. Inputs are SHA256 logged; raw outputs are intentionally not committed to Git.

## Release gates for an actual paper

1. **G0:** Identify original 646-test execution environment and exact prior settings; ensure complete variant-level shared SNP files including variants that did not reach significance.
2. **G1:** Reassess GWAS region independence and harmonize disease subtype, ancestry and effect-allele representation for each locus.
3. **G2:** Rerun SNP-level ABF + sensitivity and LD-aware multiple-signal coloc on matched GWAS/QTL input; report H3/H4 and harmonization (NOT current summary-only conditional results).
4. **G3:** Independent molecular QTL study, or coding/perturbation evidence with orthogonal confounding assessment; test pQTL, blood pressure, pleiotropy as relevant.
5. **G4:** Disease-state donor cell-type QTL/ATAC for SH3PXD2A and COL4A1/A2, contrasting healthy R3_7 expression patterns without claiming disease-specificity.
6. **G5:** Separate association, tissue detection, colocalization, mediation, and causality conclusions in the manuscript. Preserve all tested gene–tissue combinations and excluded/missing loci.

References to method:
- https://chr1swallace.github.io/coloc/reference/sensitivity.html
- https://chr1swallace.github.io/coloc/reference/coloc.abf.html
- https://rdrr.io/cran/coloc/src/R/sensitivity.R

**No new independent causal mechanism was established by this stage.**
