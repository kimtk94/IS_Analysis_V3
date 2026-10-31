# HCC Stage27 — Suter published code provenance, fullz rank preservation, and Ng Figure 2 source data identity

Date: 2026-10-31. No changes to biological/clinical significance gates without actual CI evidence. Research PR #251, governance PR #252 remain Draft.

## Confirmed upstream repository provenance

Immutable public repository: `cbg-ethz/HCC` commit `438f52271d0a86a6701baa7a96e80c516ba575c2`. Public tree has **63 blobs** and **14 tracked R files**, including `HCCfullclustering.R`, `HCCclutsering.R`, `HCC_prior.R`, `helpresfns.R`, and `README.txt`. The pinned `HCCfullclustering.R` explicitly reads `HCCinputs/HCCfullz.rds` for clustering algorithms (PCA, CIMLR, iClusterPlus); it also reads `HCCfull.rds` for MOFA analysis with in-memory variable-selection. `HCCinputs` contains pre-built `HCCfull.rds` and `HCCfullz.rds`. Published Suter manuscript identifies these as processed clustering objects. **Do not infer** from an analysis script's `readRDS` that it produced the .rds file, and do not assert that the source preprocessing code does not exist anywhere.

Stage27 source CI #41185391498 repository-code job **SUCCESS**: source 6 scripts verified by hash and scanned; named HCCfullz references 22, source loading present, zero `saveRDS(HCCfullz, ...)` writer calls within these scripts. This is a bounded source-code discovery, not a full preprocessing computation reproduction.

## Rank preservation versus raw numerical differences

Stage26 established:
- Full RNA 20,356 genes, fullz RNA 19,768 genes; common 19,768, other matrices same feature cardinality.
- Each layer has 50 source tumor biopsy columns with matching sample IDs.
- Fullz T/P/PP has per-feature SD median 1, but `scale(HCCfull)` differs substantially (T median absolute residual 0.780626; P 0.826316; PP 0.795543).

Stage27 tests, for each *exact matched row gene/protein ID* across all 50 biopsies, **per-feature Spearman** and Pearson correlation between full and fullz. Report distribution medians, 10/90 percentiles, fractions Spearman ≥0.99 / <0.9 / negative, and finite-feature counts. Values are aggregate-only; NO patient/feature names and NO per-gene rho values to workflow artifacts. High rank preservation would indicate a nearly monotonic within-feature transformation; **it still cannot establish source computational identity, normal reference, imputation, batch corrections or biological replication.**

## Figure 2 main manuscript vs Supplementary Figure 13

The publisher Source Data XLSX is SHA256 `16376438aa16a09a13db417a4ae1f4cc16b87ec2c47630aeae651004501be420` (4,446,471 bytes). Use `xl/workbook.xml` **plus** `xl/_rels/workbook.xml.rels` to map human sheet names to physical XML source sheets. Public sheet `SuppFig13g` has 5,481 ENSG IDs; this is not `Fig2b`. Exact overlap of `SuppFig13g` with current strict reviewed UniProt 2026 one-to-one 5,343 genes: 5,253 common / 228 published-only / 90 current-only. This result is a *supplementary-specific set overlap*, not original Figure 2b gene-list replication.

A separate Stage27 source profiler examines publisher-labeled `Fig2a`, `Fig2b`, `SuppFig13f`, `SuppFig13g` using only sanitized header categories and aggregate per-column numeric distributions. Original Ng text says Figure 2a CNA–mRNA per gene rho n=15,272 and Figure 2b mRNA–protein per gene rho n=5,481 (median rho=0.287). This descriptive expected median is not sufficient to infer matching 2022 historical gene IDs from a sheet with 5,481 ENSG IDs unless the source Figure 2b column or provenance explicitly connects them.

## Gate decisions

- Code source and Suter full vs fullz numeric structural audit: GO for technical provenance only.
- Whether original 2022 Figure 2b 5,481 **exact gene IDs** are present and belong to the source sheet: HOLD pending source field/header interpretation and original annotation.
- HCCfullz source generation / normal controls / missingness/feature selection full reproduction: HOLD. V03 PARTIAL.
- Original paper Spearman rho and current 48-patient 5,343-gene rho are descriptive only, V07 PARTIAL.
- Clinical, prognostic, etiology-specific mechanisms, per-gene q-values and replication: no new PASS.

Actual workflows:
- Original Stage26 (three corrected jobs PASS): https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41181715524
- Stage27 code and rank runs: https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41185391498 ; updated three-job workflow https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41185474022

New source/quantitative evidence to be appended to D-037 after CI inspection and repeat successful execution, not anticipated here.

## Definitive Figure 2 source correction — supersedes the previous premise (D-038)

The original source XLSX was parsed twice with different notions of "sheet index." The earlier D-036 corrective report incorrectly identified the *lexicographically sorted 28th worksheet XML* (`sheet8.xml`) with the real physical `sheet28.xml`. These are **not the same file**. Verification must map `sheetN.xml` to workbook relationship IDs, **and numeric-sort the file paths** instead of comparing a loop index to physical sheet number.

GitHub Stage27 #41185474022 successful independent original Figure 2 workbook inspection proves:
- Real publisher **`Fig2a` = `sheet7.xml`**, 15,272 explicit Ensembl IDs and numerical Spearman median **0.203142**.
- Real publisher **`Fig2b` = `sheet8.xml`**, 5,481 explicit Ensembl IDs and numerical Spearman median **0.287059**.
- Real `SuppFig13f` = `sheet27.xml`, 122 data rows; real `SuppFig13g` = `sheet28.xml`, 117 data rows. They are NOT the original 15,272/5,481 gene universe sheets.
- The older Stage26 ID overlap of **5,253** with current 5,343 strict one-to-one Ensembl IDs actually came from the *physical `sheet8.xml`*, so its true source is **Fig2b**. Numeric operation is valid; earlier source label in one intermediate output and D-036 is superseded.
- Stage26 parser fixed to sort by physical numeric sheet value and affirm source panel labels; new original `Fig2b` membership v3 CI https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41185931131 . The actual Figure 2b ID set is technically verified, but the original patient's per-gene rho and 2022 annotation process remain NOT reproduced.

**Updated gates:** Original 5,481 Figure2b gene ID *source identification* GO_TECHNICAL; V07 remains PARTIAL_DESCRIPTIVE_ONLY, original Nguyen sample-level correlation effect replication HOLD. Full/fullz rank discrepancy still V03 PARTIAL.

A parallel Stage27b 50×50 sample-axis identity/label-permutation test remains planned under https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41185636861 and does NOT assume source numeric equivalence.
