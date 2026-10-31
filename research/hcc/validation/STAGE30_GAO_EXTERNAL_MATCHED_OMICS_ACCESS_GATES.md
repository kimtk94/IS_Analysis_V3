# HCC Stage30 — independent Gao HBV-HCC matched RNA / proteome verification gates

Created: 2026-10-31. **This is an access and scientific validation protocol, not confirmation that numeric matrices are currently downloaded or sample joins verified.**

## Why

- The Suter/Ng Figure2b 5,253 ENSG same-source per-gene RNA/protein effect-vector Spearman **0.972556921** (48 patient mean, job #124499343862) proves strong computational concordance with the **original Ng source**, but is NOT independent patient-cohort replication. This is now GO_TECHNICAL_REANALYSIS / V07 PARTIAL, not a novel biologic claim.
- A separate HBV-related HCC cohort from Gao et al., *Cell* 2019, **159 paired tumor and adjacent liver cases**, DOI `10.1016/j.cell.2019.08.052` is a plausible validation source. Source publication: https://pubmed.ncbi.nlm.nih.gov/31585088/ .
- Exact independent public dataset roles:
  - Protein abundance: NCI PDC `PDC000198`, official study https://pdc.cancer.gov/pdc/study/PDC000198
  - Phosphoproteome: `PDC000199`, https://pdc.cancer.gov/pdc/study/PDC000199
  - RNA/transcriptomic NODE project candidate: `OEP000321`, https://www.biosino.org/node/project/detail/OEP000321
- Published independent downstream 2022 study specifically says it obtained Gao HCC proteomics from CPTAC/PDC and transcriptomics from `OEP000321`: https://pmc.ncbi.nlm.nih.gov/articles/PMC9006883/ ; **source-level Gene × Sample matrices plus matching patient IDs NOT yet verified in this project**.
- PDC API docs: https://pdc.cancer.gov/pdc/api-documentation (public endpoint POST JSON `https://pdc.cancer.gov/graphql`; latest version by PDC ID). Public study and alias record requires study UUID/version tracking, files/quantitation and case→sample→aliquot join.
- PDC download and join guide: https://pdc.cancer.gov/pdc-docs/usage . Official PDC specifically cautions that gene/protein quantitative data is at the aliquot level while clinical data is at case level and recommends the biospecimen manifest for bridging them. Exporting file manifests generates URLs expiring in **7 days**; never hardcode an expiring download URL as permanent provenance.

## Required real-data gate order (fail closed)

| Gate | Evidence (actual downloaded or API output) | GO | HOLD or STOP |
|---|---|---|---|
| G30-A PDC study | Verify PDC000198 and PDC000199 studies (PDC study UUID, version, file types, cases, aliquots); source endpoint response SHA, fetched UTC | Both study identities actually present | API refusal/unstable => ACCESS_HOLD, not zero cases |
| G30-B quantitative matrices | Actual numerical gene/protein abundance data with gene symbol/UniProt name, **aliquot sample columns**, tumor/normal sample type, missingness QC, processing units and source SHA | Values finite and enough cases for direct per-gene correlation | If only heatmaps or no numeric export => HOLD_QUANT |
| G30-C RNA | NODE OEP000321 actual expression matrix and study accession/manifest; determine FPKM/TPM/counts, reference genome, T/N coding, sample identifiers, original file SHA | Reproducible patient-linked HCC RNA file in source | NODE URL landing is not sufficient; HTTP failure => HOLD_SOURCE, never invent RNA sample access |
| G30-D join | PDC case/sample/aliquot mapping ↔ RNA exact patient labels, ensure no crossover, exclude normal references from tumor-only correlation, select 1 quantitative tumor RNA + 1 proteome measurement per independent patient (or explicitly predefine multi-biopsy aggregation) | Positive >=30 **independent** complete matched cases with metadata-documented same tissue/timepoint and quantified features (pilot threshold; target near published N159 if feasible) | Many-to-many or unknown labels => STOP cross-omics effect estimation |
| G30-E external hypothesis | Freeze Ng per-gene pathway panel PRIOR to Gao estimate; compare patient-wise gene-wise Spearman and pathway pattern, with assay/platform variability and FDR/uncertainty, no cherry-picking | Replicated sign/pathway effect with uncertainty and explicitly novel nuance | HBV-only Gao cannot validate ALD-vs-HCV differences; absence of matched RNA/P => restrict to triangulation |
| G30-F outcome/scRNA | Only after sample-identity, metadata and clinical-event codebook verified; scRNA donor-level pseudobulk for tissue-cell context only | Independent donors, covariate and temporal codebook defensible | No survival/radiomics/cell claim without source joint IDs |

## CI implementation targets

Create independent, small **read-only** GitHub Actions jobs that report public source study/version/aliquot counts and access failures WITHOUT exporting patient IDs or quantitation. Do not download large raw mass spectra until a valid quantitation path, test matrices and biospecimen manifest are confirmed. If import/public API errors happen, distinguish API schema/network failures from true biological missingness; only conditional successful access may transition G30-A/B/C. Pin source metadata artifact digest and record retrieval timestamp. Once manifests/materialized quantitation confirmed, download to controlled Drive or research server with encrypted per-patient material (not public GitHub artifacts), and publish only aggregate QC.

## Current status — updated from actual source CI 2026-10-31

- `G30-A`: **GO_CATALOG_TECHNICAL**. Official PDC GraphQL verified both study IDs PDC000198/199, each with **171 cases and 331 aliquots** in current catalog metadata. The Gao original article's **159 enrolled patients** is NOT the same denominator; technical QC needed. CI #41191160746 **3/3 SUCCESS**. PDC Protein Assembly actual listings 198=**8 files (3 TSV/5 TXT)**, 199=**9 (4 TSV/5 TXT)**, each with source MD5 metadata; CI #41191530087 SUCCESS.
- `G30-B`: **PARTIAL_FILE_CONTENT_ACCESS; HOLD_FULL_MATRIX**. Official PDC signed download URL HTTP 206 Range only, with no permanent signed link: actual 111-column TSV source headers, 520/550 preview numeric cells in each study, and phospho 334-column tables; CI #41191748823 SUCCESS. Cross-assay file header comparison CI #41192003264 SUCCESS shows up to **331 shared source header tokens** across two studies and an identically hashed 111-column header. These are numeric **source-prefix** observations, not validated full gene×aliquot quantitative matrices; no phenotype/normalisation proof.
- `G30-C`: **HOLD_NODE_RNA**. NODE OEP000321 returned HTTP 200 landing shell but the response did not contain the study ID/HCC actual project data; no RNA expression matrix, sample identifiers or permanent source digest verified.
- `G30-D`: **STOP_INFERENCE** pending true same-case RNA + protein file/aliquot join. Independent PDC official `biospecimenPerStudy` protein/phospho **metadata-only** case/sample/aliquot overlap check in CI #41192179017 (pending as of this revision). Even source PDC Protein↔Phospho matched IDs are NOT independent RNA↔Protein case data.
- `G30-E/F`: **NOT_STARTED**, no independent scientific, etiologic, prognostic, or clinical PASS.
- This project’s V02/V03/V07 remain PARTIAL; Stage29 source Ng Figure2b 5253-gene rho-vector Spearman 0.972556921 within the same Ng/Suter cohort is a *technical positive control*, not an independent validation. PR #251 / #252 remain Draft. See decisions D-045 through D-050.

## Verified PDC biospecimen metadata reconciliation (real CI 2026-10-31)

- [CI #41192179017](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41192179017) **SUCCESS**, source GraphQL `biospecimenPerStudy` (PDC000198 and PDC000199), artifact #13039740706. Both studies returned **exactly the same source response SHA256** `b09aeed5ad9fc55d885e85de2a4f226787bcb1141f9a9032639666398a5340a6`.
- EACH: **171 distinct case submitter IDs**, **331 sample IDs**, **331 aliquot IDs**. Type rows: **165 Tumor / 165 Normal / 1 OtherOrUnknown**. **160 cases with both tumor and normal metadata**. All source case/sample/aliquot identifiers nonempty. Cross-PDC Proteome/PP intersection: **171 common case IDs / 331 common sample IDs / 331 common aliquot IDs**, all metadata-only.
- **Scientific interpretation:** PDC Protein↔PP pairing **GO_BIOSPECIMEN_METADATA**, NOT quantitative gene-by-aliquot pairing and NOT RNA↔Protein same patient overlap. In particular 331 aliquots ≠ 331 independent patients; 160 tumor/normal-paired metadata cases ≠ 159 published assay-eligible patients, until exact original cohort QC and exclusion reasons verified. `G30-D` **remains STOP** for the RNA–Protein cohort because OEP000321 expression file and sample names are still unverified.
- [Stage30 exact assembly-header↔aliquot audit CI #41192598558](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41192598558) registered: read actual PDC source 331 aliquot submitter IDs privately, compare with up to six Protein Assembly file header labels per study using bounded 16KiB source reads, report only aggregate exact-match case/sample/aliquot/tumor/normal counts. Until this result exists, **do not** claim an actual full gene-by-aliquot quantitative abundance matrix or validated source aliquot columns.
