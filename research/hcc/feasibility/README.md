# HCC Multi-omics Feasibility — Evidence Gate (2026-10-30)

**Status: HOLD for definitive mechanism/recurrence study design; GO for a metadata-only pilot.** This is a research feasibility audit, not a claim that all omics data have been downloaded or joined.

## Existing related work
- Radiogenomics is **a separate project**: Google Drive `HCC_Radiogenomics_Public` and [Codex PR #228](https://github.com/kimtk94/Codex/pull/228) (still open at audit).
- This branch adds **RNA/protein + etiology + clinical + single-cell/spatial** dataset qualification to the research repo, avoiding changes to the existing IS/CKD pipelines.
- The 2026-10-22 HCC Research Proposal template reportedly prohibits AI-written submission text. This document is **internal technical evidence**, not a submission.

## Key findings and public-source cautions
1. **Ng et al 2022**: 114 patients / 122 tumors; RNA across 122 tumor samples, but proteome/phosphoproteome on **51 tumors** (paper-level intersection). Its RNA sequencing is **EGA controlled access** (EGAS00001005074); public PRIDE datasets PXD025705/PXD025836 are proteome/phosphoproteome, respectively. Supplementary/source data may enable a restricted pilot without raw RNA, but actual sample-level accessible matrices and mapping are **unverified**. https://www.nature.com/articles/s41467-022-29960-8
2. **Gao et al Cell 2019**: 159 HBV-HCC patients paired tumor/adjacent tissue in paper; public processed matrix access and RNA–protein linking not yet verified. https://pubmed.ncbi.nlm.nih.gov/31585088/
3. **TCGA-LIHC**: potential RNA/clinical, plus RPPA which is a *targeted antibody-based assay*, not whole MS proteome. Need to inspect etiology and clinical completeness. https://portal.gdc.cancer.gov/projects/TCGA-LIHC
4. **GSE14520**: 488 *microarray samples* across two platforms, not 488 unique patients; RFS/time/event columns must be verified from the actual supplement. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14520
5. **GSE149614**: 10 patients / 21 GEO samples; processed scRNA matrix exists, raw data EGA controlled. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614
6. **GSE242889**: 5 patients (3 MVI+, 2 MVI−), 10 tumor/adjacent samples. MVI is a proxy biological phenotype, **not observed recurrence**. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242889
7. **GSE319709**: 18 scRNA sample records from PBMC/tumor/peritumor, HBV/MASH, CD45+ only. Cannot infer malignant hepatocyte origin. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE319709
8. **GSE238264**: 7 processed spatial samples **after neoadjuvant cabozantinib/nivolumab**; not a treatment-naive external recurrence cohort. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE238264
9. **PXD050969**: iProX proteome/lipidome/metabolome atlas; public raw matrix access and matching to RNA are unverified. https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD050969

## Empirical public-metadata probe (2026-10-30)
[GitHub Actions #41052039696](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41052039696) successfully fetched **two small public GEO supplements** and passed **10 regression tests**.

- `GSE14520_Extra_Supplement.txt.gz`: **488 sample rows** (247 Tumor, 239 Non-Tumor, 2 Normal); `Recurr status` and `Recurr months` each populated on **242 rows**. Raw status values `1=136`, `0=106`, and median recorded recurrence months **30.4**. Do not label these as 136 events or 242 unique patients until codebook / patient/sample joins are checked.
- `GSE149614_HCC.metadata.updated.txt.gz`: **71,915 cell metadata rows**, **21 unique sample labels**, **10 unique patient labels**.
- `HBV viral status` in GSE14520 is categorical (`CC`, `AVR-CC`, `N`, `.`), **not automatically a clean HBV+ vs HBV− coding**.
- Raw EGA genomic files were **not downloaded**. Source details are preserved as non-sensitive aggregates in `datasets.json`.

Evidence of a patient-level proteotranscriptomic intersection and independent recurrence cohort remains **UNVERIFIED**.

## Stage-2 actual public Ng/Gao sample metadata (2026-10-30)

Public-source audit [Actions #41053421959](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41053421959) fetched the official Ng Supplementary Data 1 XLSX (SHA-256 `edf1d2afd3baaedd56edee8536c16f10e4a62cb6cdac929863a062671ee73d25`). The sample-level metadata **confirmed** 122 tumor biopsies / 114 distinct patients; 51 biopsies / **49 patients** are flagged RNAseq + proteome + phosphoproteome `Y` (51 samples are NOT 51 patients). The file alone does not provide or establish a cross-matrix quantitative sample join, so the hard gate stays HOLD.

Among the 49 distinct patients with all three modality flags: ALD only 13; HCV only 11; NAFLD only 7; HBV only 4; HCV+ALD 7; ALD+HBV 3; ALD+NAFLD 2; no known liver disease 2. Mixed etiologies and patient multiplicity require prespecified treatment and adequate power.

Public PDC GraphQL [Actions #41053276768](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41053276768) verified **PDC000198 and PDC000199**, each reporting `cases_count=171`, `aliquots_count=331`. These are repository counts, **not** evidence of 171 Gao paper cases (the Cell publication reports 159 paired patients), nor a matched proteotranscriptomic patient N. PDC processed quantitative matrix versus BioSino OEP000321 RNA join remains open.

**Next unlock**: obtain legally usable *processed expression matrices* with compatible patient/specimen IDs (Ng public supplement/PRIDE or authorized EGA; Gao public PDC and BioSino). Confirm batch, tissue status, unique patient count, feature/gene mapping, then only promote `verified_pair_n`. The candidate ALD-only versus HCV-only comparison is exploratory, not a powered inference.

## Stage-3 quantitative file location (2026-10-30)

CI evidence: [Actions #41054812247](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41054812247), 17 regression tests PASS; public file-metadata audit available as Actions artifact. This stage **did not** access gene×patient expression values or join RNA/protein matrices.

- **Ng**: Four publisher-hosted XLSX workbooks (Supplementary Data 2, 13, 14 and Source Data) fetched and inspected for sheet shapes. They may contain analytic summaries; sheet dimensions alone do not establish an expression matrix. PRIDE PXD025705/PXD025836 first 100 file records each fetched; three `P-PROT_HCC_SET{1,2,3}_DDA_Progenesis_Mascot_SQ_Results.xlsx` found in **PXD025836** file list. Do not assume these are complete gene-by-patient quantitative matrices.
- **Gao**: PDC GraphQL total file counts: **PDC000198=6,347** (including 8 Protein Assembly/Text files); **PDC000199=3,896** (including 9 Protein Assembly/Text files). Most are spectra/PSMs; counts do not measure usable individuals. RNA expression from BioSino has not been paired to these files.
- **Next targeted inventory**: enumerate PDC Protein Assembly filenames and the remaining PRIDE result/other files; classify results as (1) quantitative gene×sample, (2) peptide/protein identification, (3) instrument raw, (4) secondary figure/table, then compare sample identifiers with the existing Ng 49-patient measured subset.
- `scripts/hcc/probe_quant_candidates.py` contains the next candidate-filename inventory utility **but has not been executed by CI**; its output must not be treated as findings.

**Decision stays HOLD** for quantitative RNA×Protein pairing and independent clinical endpoint validation. Do not promote `verified_pair_n` in the JSON registry until a reviewed non-sensitive patient/biopsy join summary exists.

## Stage-4 actual public PDC quantitative schema (2026-10-30)

- [GitHub Actions #41056348071](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41056348071) confirmed HTTP **206 Range** retrieval of up to 32 KiB from each of four public PDC files, without downloading the full matrices or persisting PDC signed URLs. The **Proteome TMT11** TSV is **128,690,348 bytes / 667 header columns**, the **Phosphosite TMT11** TSV **110,179,920 bytes / 334 columns**. Their second data rows have matching column widths. Both files have a confirmed public quantitative TSV structure; this is NOT a proof of full quantitative matrix integrity.
- The corresponding **Proteome.sample.txt** and **Phosphoproteome.sample.txt** files (6,931 and 7,328 bytes) were retrieved and parsed; both comprise **33 data rows and 15 columns**.
- [GitHub Actions #41056621082](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41056621082) screened metadata fields against matrix header names. **Zero literal normalized identifier matches** were found. This does NOT imply no mapping: TMT isobaric label experiments use multiplex plex/channel to sample/aliquot mappings.
- **PDC official guidance**: `*.tmt11.tsv` columns represent relative log ratios (Log ratio vs Unshared Log ratio); do not treat as absolute protein abundance. `*.sample.txt` is analytical-run context; `aliquot_submitter_id` must be mapped to `case_submitter_id` using PDC Biospecimens and Study Run Metadata. See https://pdc.cancer.gov/pdc-docs/data-analysis-guides and https://pdc.cancer.gov/pdc/faq/Files_Download.
- **Next unlock**: parse PDC study-run plex/channel metadata, resolve aliquot/case mapping without persisting patient identifiers publicly, and independently verify BioSino RNA processed expression/sample keys. Only then compare matched quantitative protein and RNA values.

**Gates unchanged:** full quantitative matrix/data join = HOLD; paired RNA/protein patient count = UNVERIFIED; recurrence validation = HOLD. Successful CI validates data access and schema, not biological associations.

## Stage-5 PDC tumor/normal and 11-plex lineage (2026-10-30)

**Evidence:** [Actions #41059030994](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41059030994) ran the documented public PDC `paginatedCasesSamplesAliquots` API for PDC000198 / PDC000199, aggregated **without persisting case, sample or aliquot IDs**.

- Both studies have the exact **same set** of **171 case IDs / 331 sample IDs / 331 aliquot IDs**. These are storage-level lineages, not 171 publication-eligible cases or quantitation matches.
- Each study contains **165 Primary Tumor samples, 165 Solid Tissue Normal samples, and 1 Not Reported sample**; when grouped by case ID, there are **160 cases with BOTH tumor and normal**, **5 tumor-only**, **5 normal-only**, and **1 Not Reported** case. This is distinct from the **159 analyzed tumor–normal patient pairs in Gao Cell 2019**, indicating paper eligibility/QC exclusions still need adjudication.
- `aliquot_run_metadata` contains **33 distinct TMT experiment numbers, 363 distinct experiment×channel combinations (33×11), and 726 run metadata rows**, reflecting repeated/aliquot-specific records. Counts of run records or channels are **not independent patients**.
- Protein table has 667 columns, phosphosite table 334. **Zero direct case/sample/aliquot full-ID token matches** against the table headers, so linking quantitative columns to PDC case IDs requires the run metadata plex/channel mapping (not blind header equality).
- **V02 remains BLOCKED**: this work has not verified which exact quantitative columns correspond to patient cases, the full protein matrix values, or *any* matched public RNA gene-by-sample values.
- **Next:** inspect standardized `(experiment_number, TMT label) → aliquot` mapping for uniqueness and reference/pool channels, then determine quantitative header-to-plex naming, reconcile paper N=159, and locate BioSino OEP000321 public processed RNA expression matrix or document access limitation.

## Stage-6 Gao TMT analyte mapping vs NODE RNA access (2026-10-30)

**Executed evidence:** [CI #41062185538](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41062185538), 42 offline tests PASS. Earlier direct `studyExperimentalDesign` requests timed out; the final audit instead used working PDC `paginatedCasesSamplesAliquots` data, including `aliquot_run_metadata` with label, experiment, fraction, replicate and analyte fields. The public response was aggregated in memory without committing patient or aliquot IDs.

- **Across each PDC study**, there are 726 run-metadata rows, 33 experiment numbers, 363 `experiment×channel` slots; 330 slots contain two distinct aliquots when **analyte is omitted**. This is **not** 330 duplicate patient measurements.
- The apparent collision is explained by **two analyte classes**: `Proteome` and `Phosphoproteome`, **363 run records, 363 experiment/channel slots, 331 unique aliquots, and 0 same-slot multi-aliquot conflicts per analyte**. One must include analyte when mapping TMT channels.
- Fraction grouping alone reduces 330 collisions to 60; replicate field is **missing on all 726 records** and cannot disambiguate.
- Protein quantification table header: **667 columns**, including **330 Log Ratio and 330 Unshared Log Ratio** fields. Phosphosite header: **334 columns**, with **330 Log Ratio** fields. The extra metadata fields and quantitative column naming/clinical-alignment rules must still be audited.
- **Critical not yet achieved:** full protein quantitative column → aliquot/case identity join remains unverified; **RNA matrix is not downloaded**. The public NODE `OEP000321` page returned HTTP 200 but only a 906-byte HTML shell that **did not contain the accession or a processed RNA matrix download URL**. A webpage loading does not grant dataset access. Secondary literature reports use of an RSEM / upper-quartile-normalized processed RNA matrix, suggesting investigation of NODE's documented download workflow or legal source-article supplement, but is not direct file evidence.
- **Research decision:** public PDC multiplex lineage analysis **PARTIAL/PASS for metadata mapping only**, RNA–Protein quantitative paired **V02 remains BLOCKED**. Do not imply n=159, 160, 165, 171 or 330 paired RNA–Protein cases until an actual in-data transcript/protein matrix join and tissue/time mapping is demonstrated.

Sources: [official PDC public API](https://pdc.cancer.gov/pdc/publicapi-documentation) and [PDC TMT relative quantification FAQ](https://pdc.cancer.gov/pdc/faq/Files_Download); [NODE OEP000321](https://www.biosino.org/node/project/detail/OEP000321).

## Stage-7 PDC native gene×aliquot matrix API / Gao RNA archive fallback (2026-10-30)

**Evidence:** [GitHub Actions #41064173332](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41064173332), **49 local tests PASS; biological access is separately assessed**.

- [Official NCI PDC `quantDataMatrix` endpoint](https://pdc.cancer.gov/pdc/publicapi-documentation) documents direct gene×aliquot normalized relative quantitation (requested `data_type=log2_ratio`); this is preferred to reconstructing a protein-to-case join by guessing TMT column header formats. However, PDC000198 and PDC000199 **both timed out** after 25 seconds in the actual CI. No matrix response payload or sample identifier intersection was obtained. A future run should try PDC direct GCT/download portal, smaller queries if documented, or a separately staged bounded fetch; do not interpret 504/timeout as a missing dataset.
- [HCCDB v2.0](https://academic.oup.com/gpb/article/22/1/qzae011/7600437) archives Gao OEP000321 transcriptome as **HCCDB25**, corroborated by independent studies. It offers a plausible alternate legal processed-RNA source. In actual CI, `http://lifeome.net:809/` served **HTTP 206, HTML 6,305 bytes**, and `http://lifeome.net/database/hccdb2` the same shell. HTTPS routes failed (connection/SSL). The HTML did **not** contain literal HCCDB25/OEP000321 or a verified matrix file URL; downloading a full RNA matrix and joining patient IDs is still outstanding.
- Source-sensitive note: publication data-availability statements alone cannot establish current public download availability or patient/sample ID fidelity. **Current research decision: V02 BLOCKED**, not FAIL, due to unverified data access / mapping.

## Stage-8 HCCDB25 public frontend route inventory (2026-10-30)

A **focused** [GitHub Actions #41066051263](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41066051263) successfully tested the bounded public frontend probe (3 unit tests PASS), avoiding a full GEO/PDC rerun. HCCDB HTTP portal at `http://lifeome.net:809/` returned 6,305 bytes of HTML and referenced **6 public JavaScript assets**. Five were inspected. The small `/static/js/app.3f5b30f0.js` bundle contained the SPA route literals `/download` and `/repository`, with additional UI strings "Download" and "Repository".

**Interpretation:** these are client-side navigation routes, **not** demonstrated HTTP download endpoints or verified HCCDB25 RNA sample files. The route audit could not yet recover an expression matrix URL, a gene×sample matrix or matched case IDs. V02 remains BLOCKED. To proceed, navigate documented HCCDB25 public repository route through an authorized browser/workflow or obtain an official exported HCCDB25 expression matrix and sample manifest; alternatively request a direct access link from the data custodian. Never guess a dataset download URL from a route name.

## Criteria (pilot engineering thresholds, NOT universally accepted clinical criteria)
- At least one observed, documented **patient/specimen-level RNA↔protein join with N≥30** (paper abstract is insufficient).
- At least two verified etiology categories **each N≥10 within those joined data**. Need inspect overlap and confounders (stage, treatment, batch, platform).
- For optional recurrence prediction: one eligible event-defined dataset with ≥30 events, **independent** external cohort with ≥20 events. These are *minimum exploration triggers*, not statistical power guarantees.
- RFS, DFI, PFI, OS, MVI, and treated response are different endpoints. TCGA-CDR recommendations: https://pmc.ncbi.nlm.nih.gov/articles/PMC6066282/

No `verified_*` flag may be changed without a non-sensitive supporting audit artifact with counts and join rules. Never commit raw restricted files or identifiable patient keys. `HOLD` = evidence missing, not `NO GO`.

## CI usage
```bash
python -m unittest discover -s tests -p 'test_hcc_feasibility.py' -v
python scripts/hcc/feasibility_audit.py --out-dir results/hcc/feasibility
# Optional metadata URL reachability only (NOT a file download):
python scripts/hcc/feasibility_audit.py --probe-links --out-dir results/hcc/feasibility
```
GitHub Actions: `hcc-multiomics-feasibility.yml` on PR / branch main change / manual dispatch; writes `audit.json` and `audit.md` as artifacts. Success of CI means the guardrails and generator ran; **not** biological/clinical data availability.

## Next concrete tests (in order)
1. Review Ng 2022 paper Source Data and PRIDE file listing; establish whether 51 tumor RNA & protein values can be paired **without** EGA controlled data. If no, request authorized data and consider alternative public pairs.
2. Download **metadata only** from Gao/TCGA/GSE14520; calculate distinct patient counts, eligible HBV/HCV/MASH vs alcohol categories, missingness, outcome definitions.
3. Check outcome follow-up/event counts for eligible discovery/validation cohorts, separated by training/validation and primary therapy.
4. Only after stage 1–3, design low-cost RNA/protein pathway concordance pilot; scRNA/spatial are contextual validation, not artificially patient-matched modalities.

Repository registry: `research/hcc/feasibility/datasets.json`.
