# IS GSE331114 8-sample cell QC / reference-annotation gate — 2026-11-11

**Scope:** real public human ischemic-stroke brain CD45+ enriched single-cell data; **not** validated cell-type or genotype/eQTL analysis. This internal research note is not an AI-authored submission (course form prohibits AI submissions).

## Evidence lineage (true original data)

- Official GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE331114
- [Eight full MEX source QC #43010007530](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43010007530): 4 AIS GSM9739251–4 and 4 control GSM9739255–8; eight raw 10x MEX matrix/feature/barcode trio QC PASS.
- [8-sample gene CPM synthesis #43012648250](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43012648250): 26 genes (15 historical IS candidates, 11 immune markers), per-sample CPM; descriptive n=4+4 results; no P/FDR. CD45+ cell composition was not adjusted.
- [Per-cell re-download and QC #43013961145](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43013961145): **8/8 jobs SUCCESS**; all three files for each sample independently re-downloaded and SHA-256 matched previous audit; per-cell detected genes, UMIs and mitochondrial fraction measured; exploratory immune-marker co-detection.
- [Eight-sample marker synthesis #43014520677](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/43014520677): **SUCCESS, offline unit tests 6/6 PASS**. Downloaded exactly the original eight sample marker-QC artifacts. The GitHub job log first shows *synthetic* unit-test fixture data; only the second synthesis summary after `actions/download-artifact` represents **real GEO data**.

## Observed cell QC (all eight libraries)

- Per-GSM QC passing barcodes **98.23–99.87%**, with prespecified heuristic filters: 200–7,500 detected genes and mitochondrial UMI fraction strictly below 20%. Such high pass rates are expected to be evaluated against possible upstream Cell Ranger filtering; this does NOT establish doublet-free or batch-free cells.
- Per-library median mitochondrial fraction **0.945–2.024%**. Some samples have major differences in library size and recovered barcodes; group composition is not balanced or randomized by sample.
- All eight processed MEX source SHA-256s and aggregate UMI totals matched prior recorded values. No variant or donor-specific eQTL data in this transcriptome series.

## Real GEO data: descriptive per-GSM immune marker co-detection

**Marker support** = at least two distinct detected genes in one of five prespecified panels among QC-passing cells, with *exclusive* single-panel support; cells showing multiple panels are marked MULTIPLE_PANELS and cells without an exclusive qualified panel UNASSIGNED. This is NOT cell-type annotation or immune-cell abundance estimation. Marker detection varies with UMI depth, tissue injury, ambient RNA, doublets, and marker selection.

| Exploratory exclusive marker panel | AIS median fraction of QC barcodes | Control median fraction | Difference |
|---|---:|---:|---:|
| T markers (CD3D, CD3E, TRAC) | 5.4% | 3.9% | +1.6 percentage points |
| B markers (CD79A, CD79B, MS4A1) | 1.3% | 0.4% | +0.9 pp |
| NK markers (NKG7, GNLY, KLRD1) | 3.1% | 1.1% | +2.0 pp |
| Myeloid markers (LST1, LYZ, FCN1, CD14) | 37.2% | 32.9% | +4.3 pp |
| Microglia markers (P2RY12, TMEM119, SALL1) | 0.4% | 1.9% | −1.5 pp |
| Multiple qualified marker panels | 9.1% | 7.1% | +2.0 pp |
| Unassigned by threshold | 39.2% | 53.2% | −13.9 pp |

**NOT statistically significant tests**: n=4 original AIS versus n=4 control GEO library samples; biological donor independence/clinical covariate comparability not confirmed. These fractions are not actual proportions of validated T, B, NK, myeloid, or resident microglia cell types.

The pan-microglial markers P2RY12 and TMEM119 may be downregulated in activated microglia; using a 2-marker rule can generate false negatives. Conversely, multi-panel co-detection does not alone imply doublets. Need independent reference annotation and disease-state sensitivity.

## Why formal celltype pseudobulk remains HOLD

- Only 10x MEX count matrices + GEO sample-level metadata validated. **No original author-exported Seurat cell-level type labels** confirmed in available MEX; the submitter stated Seurat was used, but that does not mean its per-barcode annotation files were deposited.
- Biological replicate claim in GEO series description does not give unique donor IDs, sex/age, time from stroke, lesion site, tissue-processing covariates, or batch linkage. Metadata field audit in the source artifacts should be inspected without inferring information from location.
- To perform **n=4+4 donor × cell-type** pseudobulk: independently verify donor linkage, annotate each cell with externally validated reference transfer (e.g. Azimuth/CellTypist plus manual marker review, verify reference tissue/disease compatibility), exclude QC/doublets, establish minimum cells in at least 3 biological replicates per group, predefine covariate model within tiny n, sum raw counts over cell type × donor, then perform exploratory edgeR/DESeq2 with proper multiple testing and uncertainty.
- A marker threshold or entire-sample CPM does not substitute for the above. **Do not** treat cell barcodes as n. **Do not** claim post-stroke transcriptional changes mediate genetic susceptibility.
- Methodological grounding: [Seurat pseudobulk guidance](https://satijalab.org/seurat/articles/de_vignette), [Harvard Chan single-cell pseudobulk guide](https://bcbio.github.io/singlecell-reports/03_differential_expression/scRNA_pseudobulk.html); [microglial marker limitations](https://pmc.ncbi.nlm.nih.gov/articles/PMC8498770/); [Azimuth reference markers](https://azimuth.hubmapconsortium.org/references/human_pbmc/).

## Study-level decision

- **GSE331114 immune-response context:** **GO for exploratory phenotype-associated expression QC, and broad marker-support characterization**. True differential cell-type abundance and cell-type-specific DEG **HOLD**.
- **ROSMAP:** retain genotype-linked brain eQTL role after independent source/data access checks; human stroke immune scRNA does NOT replace variant-level eQTL for colocalization.
- **GTEx:** retain normal tissue expression and possible tissue eQTL; distinguish normal tissue expression from disease reaction.
- **EUR/EAS IS GWAS:** original GIGASTROKE files passed CRC, structural QC, publisher MD5. Cross-population eQTL matched locus and LD availability are separate gates. MAGMA gene analysis is secondary; LDSC-SEG enrichment conditional on h2 and LD score reference quality.
- **Newness:** human stroke response annotation is not itself novel, and observing gene expression after stroke is distinct from genetic susceptibility or mechanistic mediation. Preserve previous 2026 literature overlap matrix.

## Reproducibility files

- `scripts/is/audit_gse331114_cell_marker_qc.py` and `tests/is/test_gse331114_cell_marker_qc.py`
- `.github/workflows/is_gse331114_cell_marker_qc.yml`
- `scripts/is/aggregate_gse331114_cell_marker_qc.py` and `tests/is/test_aggregate_gse331114_cell_marker_qc.py`
- `.github/workflows/is_gse331114_marker_synthesis.yml`

CI SUCCESS means the real files/marker metrics were processed. **It does not establish full biological validity**.
