# IS Phase9F-E R3_7 — Donor-aware pseudobulk protocol

**Status: REAL_DONOR_RUN_COMPLETE; INDEPENDENT_RECALC_PASS (2026-10-10).** The prespecified comparison inventory and QC gates below were frozen before the biological dataset run.

## Purpose

Assess donor-to-donor consistency of expression-localization hypotheses in
GSE256493 adult/control temporal-lobe vascular/perivascular Seurat data.
Avoid interpreting 80,515 cells as 80,515 independent biological replicates.

- Source: GSE256493, adult control temporal lobe, RNA counts.
- Source RDS SHA256: `337c3d82c070c7ea71ad21114af847e5f4c76ab199e70c2ebcdce6b03db104dd`.
- Donor key: `Patient` (TL5, TL6, TL8, TL7, TL9, TL10); cell types: author's `Cellclusters` (11).
- The author class `Microglia and Macrophages` cannot resolve the two cell identities.
- R3_6 feature audit: FGF5 and C4orf22 not matched to RNA feature
  rownames; these genes receive `FEATURE_NOT_IN_MATRIX` and NA expression.
- No stroke disease-state contrast, SNP genotypes, cell eQTL, causal mediation
  or differential-expression hypothesis test is performed.

## Prespecified donor-level workflow

1. Verify Seurat RNA assay **raw counts** and that exactly one counts layer
   maps to author cell annotations. Abort on missing Patient IDs, NA, mismatched
   cell axes, negative libraries or more than one count layer.
2. For each measured target gene and donor x cell-type subset, sum raw counts
   and the cell-level RNA library sizes, then calculate pseudobulk CPM as
   `1,000,000 * sum(target_counts) / sum(library_counts)`.
3. Compute donor x cell-type target detection fraction
   `n_cells_with_target_count_gt_zero / n_cells` and `log2(CPM+1)`.
4. Require at least **20 cells** and nonzero library size per
   donor–celltype observation; observations failing the gate are shown but
   are ineligible for paired comparisons.
5. For each **predeclared** gene/author-celltype pair, compute paired donor
   log2(CPM+1) differences, paired donor count, median difference, number
   of positive/negative/tied differences, and positive fraction.
6. Require at least **4 donors passing both cell-type coverage gates** for
   `DESCRIPTIVE_PAIRED_QC_PASS`. Otherwise label
   `INSUFFICIENT_PAIRED_DONORS`; do not assign an exploratory advantage.
7. **No inferential P values** from these comparisons. Donor compositional
   imbalance, sample processing, age/sex, technical variation and RNA
   coverage remain possible confounding explanations.

## Predeclared comparisons

| Gene | Reference class A | Reference class B |
|---|---|---|
| SH3PXD2A | Oligodendrocytes | Microglia and Macrophages |
| SH3PXD2A | Fibroblasts | Microglia and Macrophages |
| SH3PXD2A | Smooth muscle cells | Microglia and Macrophages |
| COL4A2 | Smooth muscle cells | Endothelial cells |
| COL4A2 | Fibroblasts | Endothelial cells |
| COL4A2 | Pericytes | Endothelial cells |
| COL4A1 | Fibroblasts | Endothelial cells |
| COL4A1 | Smooth muscle cells | Endothelial cells |
| CALHM2 | Endothelial cells | Microglia and Macrophages |
| ALDH2 | Endothelial cells | Microglia and Macrophages |

Comparisons were frozen before running R3_7 on the biological dataset.
No multiple-test-adjusted significance is claimed.

## Outputs (new isolated folder)

Folder:
`MASTER_DEGREE/IS_COLAB/phase9f_e_human_vascular/results/R3_7_DONOR_PSEUDOBULK_20261010/`

- `HUMAN_DONOR_PSEUDOBULK.tsv` — gene x donor x author cell-type
  counts, library size, detection and CPM, with coverage flags.
- `HUMAN_DONOR_CELLTYPE_COVERAGE.tsv` — number of eligible and represented donors.
- `HUMAN_DONOR_FEATURE_STATUS.tsv` — gene coverage including missing features.
- `HUMAN_DONOR_PAIRED_COMPARISONS.tsv` — ten predeclared, donor-paired descriptive comparisons.
- `HUMAN_DONOR_PSEUDOBULK_README.md` — threshold provenance and inference boundary.
- Original R3_6 outputs are regenerated in a **new folder** for comparison.
- `HUMAN_RUN_MANIFEST.json` — added donor audit metadata, independent from genetic results.

## Run and QC

Run the single-cell notebook in Colab with the original Drive input.
The notebook regenerates original R3_6 outputs and donor outputs together;
this is a **real dataset run**, not accomplished by building the notebook.

Build script: `scripts/is/build_phase9f_r3_7_donor_colab.py`.
Donor helper source: `scripts/is/phase9f_donor_pseudobulk.R`.
Synthetic test:
`Rscript tests/fixtures/is_phase9f_donor_pseudobulk_selftest.R scripts/is/phase9f_donor_pseudobulk.R`.

The generated notebook has passed Python `ast.parse` and full embedded R
`parse()` static tests. The entire 209-test existing regression suite
passed in the audited server environment.

**Release gate:** All five new donor artifacts must exist and be nonempty,
analysis must reach `PHASE9F_E_HUMAN=COMPLETE`, and sample sizes / donor
coverage must be examined before any biology claims or web promotion.

## Post-run verified observations (after pre-registration)

- Drive original R3_7 results: https://drive.google.com/drive/folders/1JJJm7xkze2ene_vdDcm6nx_Ni-sU84Yi
- Server immutable source import: results/is/stage5_functional/phase9f_e_human_r3_7_donor_import_20261010
- Server independent analysis: results/is/audits/phase9f_r3_7_donor_reaudit_20261010_v1
- Verified Drive audit folder: https://drive.google.com/drive/folders/1rj44uRREGFwpX6_KI5NaR7MTf5ICPXCb
- Independently checked 80,515 raw cells, six donor labels, 585 gene-donor-celltype rows, 99 gene-celltype coverage rows and ten predeclared paired comparisons.
- **10/10** comparisons meet minimum 4 paired donors; **9/10** show a unanimous direction for CPM, versus **8/10** for detection fraction.
- SH3PXD2A (Oligodendrocytes, Fibroblasts and Smooth muscle cells vs combined Microglia and Macrophages) is CPM-higher in all eligible 5-6 donors.
- COL4A2 Pericytes vs Endothelial cells: CPM is higher 6/6, but detection fraction is higher 3/6; CPM magnitude and fraction of expressing cells answer different questions.
- ALDH2 endothelial vs combined immune class shows 3/6 positive CPM and 3/6 negative CPM paired differences; the earlier pooled expression contrast was not consistently reproduced at donor level.
- **15,903 of 15,942 Neurons (99.76%) belong to donor TL8**; donor-level claims for this class are not supported by robust representation. Neuron progenitor, Stem cell and Astrocyte categories are also donor-skewed.
- Five-source input provenance is SHA256-frozen in the separate output JSON. The published data remains **normal brain descriptive localization**, not stroke-specific differential expression, cell eQTL or genetic causal evidence.
- Independently generated paired-delta and donor-imbalance SVG figures are outside Git; source data and full results are not committed.
