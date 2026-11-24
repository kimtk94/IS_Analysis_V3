# XING 2023 — PMC 2026 실제 정량 데이터 검증 (2026-10-30)

**Source:** Xing et al., *Cell Reports Medicine* 2023, DOI https://doi.org/10.1016/j.xcrm.2023.101315
**Authoritative official source:** https://pmc-oa-opendata.s3.amazonaws.com/README.txt ; PMCID `PMC10783603.1`
**Actual successful GitHub Actions:** https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41075092093 ; https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41075680502 ; https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41075928463

## PROVEN biological matrix evidence

- **mmc3.xlsx** (supplementary Table S2) successfully downloaded from metadata-validated anonymous NCBI PMC Cloud HTTPS: **16,978,396 bytes**. SHA-256 **`7da732a933eba7a84936b8033d37564fd27125e29d059f9d7e306cf456965e17`**.
- Workbook OOXML verified; 5 sheets: `Description`, `1. 6512 proteins`, `2. Clinicopathological charater`, `3.Pathway characteristics`, `4. 761 subtype-specfic proteins`.
- On actual protein sheet: **6,512 data rows**, **307 populated Excel cell positions per row**, including a working hypothesis of 3 feature annotations + **304 quantitative specimen columns**. The sheet's format extent `A1:XFD6514` overstates actual populated columns; never treat XFD as biological measurements.
- Actual second-row **304 sample header identifiers** were inspected in memory. All **304** matched an explicit anonymous `numeric_patient_key + T/N` coding, yielding **152 paired tumor/normal sample-header keys**, no header-only unmatched samples. **This is a header-level pair**, not independently matched clinical patient records or RNA expression.
- Full protein matrix **6,512 × 304 = 1,979,648** potential finite numeric cells. **1,311,220** cells are finite numeric measurements (**66.2350%**, missing/other **33.7650%**). Per-sample quantified protein feature counts: min **2,951**, median **4,349**, max **5,459**. **4,383 proteins** were measurable in at least 152 of 304 samples; **1,102 proteins** measurable in all 304.
- The study reports **160 included patients, proteomics QC 152, phosphoproteomics 132, RNA-seq 57, WES 58**. These different populations must not be conflated. The current confirmed independent protein sample-header pairs are 152, while **verified paired patient-level RNA–Protein N is still unknown**.
- Original raw workbook was used transiently in CI and **not committed**. Aggregated sanitized audit JSON is uploaded as a GitHub Actions artifact. Publisher **CC BY-NC-ND** licensing: do not redistribute the full workbook without permission.
- Another attempt on legacy PMC HTML endpoint yielded HTTP 200 with HTML, not XLSX. Proper 2026 PMC Cloud metadata references `s3://pmc-oa-opendata/...` objects: convert using the **official documented HTTPS mapping** and validate content, never guess a DOI-specific path.

## Scientific gate interpretation

**Protein-only quantitative matrix: source + header + numeric count validated.** This is a legitimate **P0 protein expression feasibility PASS** for Xing 2023, but not V02 paired multiomics PASS. Next verify official clinical metadata Table S1 `mmc2.xlsx` and RNA differential Table S6 `mmc7.xlsx`, then perform in-memory source-to-clinical sample ID join with audited denominators. A differential-expression output table is **not** the raw RNA gene-by-sample matrix. Do not treat paper-reported RNA=57 or proteome=152 as verified RNA–protein intersection.

**Novelty warning:** Xing et al. already demonstrated three HCC proteomic subtypes and external validation across Jiang/Gao/Ng. Original scientific contribution must involve predeclared etiology-specific RNA–Protein pathway interaction, independent replication, and robust patient-level cell context; subtype clustering alone is insufficient.
