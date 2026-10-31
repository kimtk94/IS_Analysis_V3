# HCC Stage26 — Original 5,481-gene membership and public RDS transformation identity

Date: 2026-10-31. State: source provenance audit **in progress** until live CI result is inspected. Research PR #251, governance PR #252 remain Draft.

## Original paper evidence hierarchy

Ng et al. Nature Communications 2022, DOI https://doi.org/10.1038/s41467-022-29960-8 explicitly report **5,490** cross-comparable genes/proteins for differential dysregulation quadrants and **5,481** genes for per-gene mRNA–protein correlation (CNA–mRNA 15,272). These are different denominators; neither is an exact source gene-ID list by itself. The published Figure 2 legend says downloadable **Source Data** contain values for panels. The Nature article also makes **Supplementary Data 4** available to support Figure 2 / correlation pathways; the direct official supplement URLs resolved via the paper:
- Figure source XLSX: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-022-29960-8/MediaObjects/41467_2022_29960_MOESM19_ESM.xlsx
- Supplementary Data 4: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-022-29960-8/MediaObjects/41467_2022_29960_MOESM7_ESM.xlsx
- Supplementary Data 2: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-022-29960-8/MediaObjects/41467_2022_29960_MOESM5_ESM.xlsx

Stage26 source probe retrieves the bounded XLSX original through publisher or PMC, independently verifies ZIP structure, scans all sheet rows for explicit ENSG-like identifiers, and only exports aggregate sheet counts and salted-free SHA256 of an **ordered deterministic Ensembl gene set** if it is literally present. A 5,481 row count or number of ENSG-like IDs is **not** direct proof of exact original 2022 set identity. If only gene symbols are available, full symbols→2022 GRCh37 Ensembl mapping and exact gene identity remain HOLD. The primary source may not expose the complete original ID universe; do not invent it.

## Empirical HCCfull vs HCCfullz

Public Suter repository `cbg-ethz/HCC`, source commit `438f52271d0a86a6701baa7a96e80c516ba575c2`, original `HCCfull.rds` SHA256 `461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a`. The original 50 biopsy source matrices have five modalities M/CN/T/P/PP.

Stage26 computes per modality:
- original/fullz file SHA256 and size; exact row feature ID/column sample ID order; numeric/finite/missing coverage;
- exact same numeric cell counts;
- direct residual if `fullz` is base-R per-feature (row) z-score `(full-rowMean)/rowSD` versus standardized **by sample column**, with median/max residuals and proportion below 1e-5;
- fullz empirical per-gene centered row mean and unit-variance SD checks.

A nearly exact row z-score relationship is evidence of an observable numerical transformation, not confirmation that the source original Ng counts / normal-reference signal or patient/batch imputation pipeline was reproduced. If relationship differs, investigate inter-stage source processing and do not auto-override. No feature or patient IDs, individual expression matrices or gene-wise correlations can be uploaded.

## Verification gates

- Paper source file provenance: only source workbook retrieval/format/aggregate counts GO; original 5,481 gene identities HOLD until matching explicit IDs and independently reviewed mapping policy.
- Full/fullz technical relation: evaluate actual regression residual/ordering; still V03 PARTIAL until Suter source computation and Ng normal-control lineage verified.
- Stage24 48-patient median rho 0.284086 and Stage25 LOPO 0.273127–0.292900, bootstrap empirical 0.243922–0.329928 remain descriptive (no paper exact replication, clinical effect or etiology inference).
- Separate Actions workflow jobs can fail independently. Action #41180971741: https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41180971741 .
- After successful actual run, append next Decision entry with observed results, CI job conclusions, artifact SHA and any verified per-sheet candidate IDs. Do not mark any biological significance PASS based on numeric similar medians.

## Stage26 executed evidence and source sheet-name correction (2026-10-31)

**Correction to exploratory premise:** The original Source Data workbook has distinct workbook sheet metadata and physical worksheet XML indexes. `xl/workbook.xml` sheetId is **not** the same as `sheetN.xml`; only `xl/_rels/workbook.xml.rels` resolves this relationship. The physically numbered `sheet27.xml` (15,272 distinct ENSG) is the **SuppFig13f** worksheet; `sheet28.xml` (5,481 distinct ENSG) is **SuppFig13g**. **Neither is Fig2a/b**. Count matching to the paper's main Figure 2 denominators cannot establish exact original Figure2 gene membership. The earlier D-035 interpretation is superseded by D-036.

Corrected source-to-current exact Ensembl string overlap (original *SuppFig13g* 5,481 vs current Suter strict reviewed-UniProt 5,343): 5,253 shared, 228 publisher-only, 90 current-only; publisher match fraction 95.8402%, current match fraction 98.3156%, Jaccard 0.942919. This is a reproducible published-supplement gene-set comparison, not exact reproduction of the original Fig2b per-gene correlation data set. Old 41181240430 arithmetic passed but mislabeled its source as Figure2b; schema v2 and corrected job https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41181715524 explicitly fix the source labels and keep Fig2b verification HOLD.

Empirical Suter RDS lineage CI #41181240430: RNA full=20,356, fullz=19,768, overlapping=19,768; protein full/fullz=5,585, PP=7,889, CN=20,339, M=3,519. Fullz T/P/PP individual-feature row SD median is 1 but the median absolute residual for direct row-z transform of HCCfull is T=0.780626, P=0.826316, PP=0.795543 (thus **NOT** direct `scale(HCCfull)`). HCCfullz SHA256 `9c1e4e05472a0aa477b81aa7b4c4e1b27809df9a2b71c1f57e247bb031263007`; full SHA256 `461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a`.

**Status:** technical source and transform screening GO, exact original Figure2b IDs and Ng original analytic processing HOLD, V03 PARTIAL, V07 PARTIAL, manuscript biology PASS NONE. Next identify Figure2b underlying numeric data or upstream script/2022 annotations and inspect HCCfullz generation and feature selection independently.
