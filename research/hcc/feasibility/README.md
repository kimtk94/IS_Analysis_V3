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

## Stage-9 source-paper calibration: prior **158-patient RNA/protein integration** (2026-10-30)

**New high-priority paper:** Cheng et al., *Briefings in Bioinformatics* 24(2), 2023, [DOI 10.1093/bib/bbac601](https://doi.org/10.1093/bib/bbac601) ([publisher methods](https://academic.oup.com/bib/article/24/2/bbac601/7025463)). **This study already integrated the Gao HBV-HCC cohort at 158 patients after omitting samples with insufficient multi-omics data**. It states RNA quantification from BioSino `OEP000321` was **UQ_FPKM**; proteome and phosphosite came from Gao original supplementary files. Published differential-expression methods convert FPKM to TPM, apply log2(TPM+1), and use limma, BH correction, plus relative protein abundance normalization.

**IMPORTANT N-MATCH CORRECTION (publisher Methods):** The paper selects **158 pairs** for overall multi-omic analyses, but its **integrative clustering actually uses 156 tumor samples** after **2 unpaired samples were excluded because the RNA and proteomic sample identifiers did not match**. Treat `158` and `156` as distinct published populations. Our V02 should report both a pre-join eligible N and post-join sample-matched N. Publisher methods: https://academic.oup.com/bib/article/24/2/bbac601/7025463 .

- **Published, not independently replicated here**: 1,977 differential mRNAs, 1,273 differential proteins; 408 RNA/protein concordant deregulated genes (224 jointly upregulated, 184 jointly downregulated), 3 integrated molecular subgroups, and external validation. The 408-gene signature is a **required literature positive control**, not our thesis's novel discovery.
- **Crucial distinction**: Gao original **159 paired patients**, PDC **160 tumor+normal lineage cases** among 171 cases, and Cheng analytic **158 paired patients** are different cohort selection/count concepts. The Cheng 158 figure is a publication-reported result, **NOT an intersection of files we have downloaded/joined**.
- **HCCDB25 RNA external archive corroboration**: a separate study describes `OEP000321` as **158 tumor samples** in the HCCDB v2 RNA reanalysis ([Frontiers 2024](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1436131/full)). Gene expression file download and patient IDs remain unverified here.
- **Novelty guardrail**: repeating proteotranscriptomic subtype clustering in the Gao HBV cohort, or the 408-gene RNA/protein concordance, is prior art. Keep RQ1 targeted at cross-etiology transferability + independent replication + patient-level cellular interpretation.
- **V02 status still BLOCKED**: independent paired **quantitative matrix** join is mandatory regardless of earlier publications claiming paired data. The reproducible target is **158 paper-eligible pairs**, subject to filtering and sample ID verification, *not an assumed output count*.

**Additional first-party literature portal route**: a 2024 study explicitly reports [`http://lifeome.net:809/#/repository?type=buke`](https://www.frontiersin.org/journals/cell-and-developmental-biology/articles/10.3389/fcell.2024.1407428/full) for HCCDB datasets including HCCDB25. This is a web-app view, **not a direct RNA file link**.

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

## Stage-12 Xing 2023 clinical completeness / candidate RFS-OS schema (2026-10-31)

**Executed GitHub CI:** [source-matched column completeness #41092590299](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41092590299), [candidate event/duration QC #41092972318](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41092972318). All steps and synthetic regression tests passed.

- Verified the official PMC Cloud mmc2.xlsx and mmc3.xlsx SHA256. Exactly **152 distinct proteome T/N sample header pairs** correspond to **152 unique same-row T/N sample labels** in mmc2 Table S1 detailed clinicopathological worksheet. **31 non-ID clinical columns (C:AG)** were audited. Patient IDs, row-level clinical attributes and complete XLSX workbooks were not published.
- Header-based values present among matched rows (denominator 152): Age 152, Sex 152, AFP 152 (143 numeric), HBV-associated H 143, HCV-associated M 138, MVI 152, tumor size 152, TNM V/W 152 and BCLC X 152. Presence does not certify coding or valid category levels.
- Candidate RFS-status/duration **AC/AE:** 152 numeric pairs, observed codes 0=26 and 1=126, all positive durations, median uninterpreted duration 7.1333. Candidate survival **AF/AG:** 152 pairs, observed codes 0=82 and 1=70, all positive durations, median uninterpreted duration 30.4833. Neither code direction nor time unit/time zero/censoring has been verified; **these are NOT proven recurrence/death event counts**.
- **Gates:** V02 BLOCKED (no RNA expression-by-patient matrix join); V12 PARTIAL (endpoint codebook missing); V13 BLOCKED (no independent locked survival score or external validation). No survival model has been fit.
- Scripts: scripts/hcc/probe_xing2023_clinical_layout.py and scripts/hcc/probe_xing2023_clinical_completeness.py. Tests: tests/test_hcc_xing2023_clinical_layout.py and tests/test_hcc_xing2023_clinical_completeness.py. Artifact: xing2023_clinical_completeness.json. Evidence registry: research/hcc/feasibility/XING2023_EXTERNAL_QC.json.
- **Next:** verify Table S1 source codebook; then calculate genuine event N and complete-case adjustment variables without publishing patient records. Independently procure lawful RNA gene-by-sample values for V02.


## Stage-13 Xing 2023 clinical complete-case prescreen — 2026-10-31

- **Actual GitHub Actions PASS:** https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41155499877 (official mmc2/mmc3 source hash checks, 152/152 exact clinical-row sample-label mapping, synthetic tests, aggregate-only source analysis). Implemented in scripts/hcc/probe_xing2023_clinical_complete_case.py and tests/test_hcc_xing2023_clinical_complete_case.py.
- Predeclared baseline variables (age, sex, tumor size, BCLC, candidate event/status + time): **152/152 numeric/presence-complete** for each candidate endpoint, RFS and survival.
- Extended clinical covariates (adds numeric AFP, MVI, TNM stage): **143/152 candidate complete**, 9 lost. AFP itself has 143 finite numeric values though 152 have a nonblank value.
- HBV/HCV-related fields together have **132/152** finite numeric paired coverage. HCV-related field alone 138, HBV-related field alone 143; **132 is the actual intersection**, not 143+138-152 inferred from marginal counts.
- Extended covariate set including both viral marker fields and candidate endpoint: **125/152 complete**, 27 lost; these are not 125 patients of any verified HCC etiology or a matched RNA/Protein count. This is a sample-key-matched *protein/clinical-only* cohort.
- Conditional candidate raw endpoint code frequencies (NOT recurrence/death events): extended N=143 RFS 0=25, 1=118, survival 0=78, 1=65; extended + viral N=125 RFS 0=20, 1=105, survival 0=71, 1=54. They must not be labelled event rates or used for survival inference before source-specific codebook adjudication.
- **Decision:** V12 PARTIAL, V13 BLOCKED, V02 BLOCKED. This completes an engineering data-completeness preflight, NOT clinical event semantics, confounder adjustment, Cox/KM survival claims, external recurrence validation, or multiomics interactions.
- **Next P0:** validate true RFS/OS event definitions/0–1 direction, unit, time-zero, censoring and eligibility directly in the publication/supplemental source; verify HBV/HCV marker clinical meaning and subgroup counts; compute reviewable biologically eligible complete cases; separately procure actual patient-level RNA quantification for V02.


## Stage-14 Xing 2023 HBV/HCV code distribution and etiology feasibility (2026-10-31)

**Actual source-audited CI PASS:** https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41155985179 ; in-memory, SHA256-pinned official `mmc2.xlsx`/`mmc3.xlsx`, 152 unique protein-to-clinical T/N row alignments, 4 synthetic regression tests. Code: `scripts/hcc/probe_xing2023_etiology_codebook_preflight.py`. No patient IDs, original spreadsheets, raw clinical attributes, or unallowlisted categorical values have been published.

- HBV-associated worksheet column H: **143/152** finite numeric, exclusively binary 0/1, observed code **1=131, 0=12**, 9 missing/non-numeric.
- HCV-associated worksheet column M: **138/152** finite numeric, exclusively binary 0/1, observed code **1=1, 0=137**, 14 missing/non-numeric.
- Both markers numeric in **132/152** source-matched clinical rows, but this is NOT the number of RNA/protein pairs or cases assigned to an etiology.
- Source liver disease text column P contained **0 directly recognizable literal HBV/HCV/ALD/MASH etiology labels** among 152 values under conservative allowlisted matching; other source field content remains unclassified and is not printed.
- **Scientific caution:** the original worksheet marker *codebook* and positivity meaning of numeric 0/1 are unconfirmed. In particular, 131 HBV code-1 cases and **only 1 HCV code-1** must not be labelled 131 infected or 1 infected patients. If code 1 is positive, the HCV comparison has gross under-representation. This is a **strong conditional data-sufficiency warning**, not a biological or etiologic result.
- **Decision:** Xing is a compelling protein + clinical QC and possible *within-cohort* outcome reference; it is **HOLD for cross-etiology HBV vs HCV interaction** without official codebook and independent cohort support. V04 PARTIAL; V02 BLOCKED; V12 PARTIAL; no interaction, RFS or etiologic population inference approved.
- **Follow-up priority:** confirm binary code interpretation and actual etiology categories with official Table S1 headers/source codebook; count patient subsets after licensed access; prefer alternative etiology-diverse Ng2022 (49 distinct patients flagged for three assays, but quantitative join unverified) and other appropriately powered cohorts. Don't use a source-paper's published patient counts as a replacement for real quantitative RNA-protein joins.

## Stage-16: publicly available Ng/Suter multilayer quantitative RNA + protein RDS (2026-10-31)

**New source route (not EGA bypass):** Suter et al., *PLOS Comput Biol* 2022 DOI [10.1371/journal.pcbi.1009767](https://doi.org/10.1371/journal.pcbi.1009767) releases preprocessed reduced and non-reduced Ng HCC multiomics in [cbg-ethz/HCC](https://github.com/cbg-ethz/HCC). Suter Methods §2.1 confirm 51 biopsies/49 patients originally, one hypermutated biopsy excluded, yielding **50 biopsies/48 patients**. The initial paper count n=48 was subsequently reproduced by the exact source-identity audit described in Stage-17.

**Real source-binary CI evidence:** [GitHub Actions #41160494201](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41160494201) successfully downloaded and parsed RDS objects, checked source SHA256, five numeric modalities, sample label alignment, and source format orientation; immutable source commit `438f52271d0a86a6701baa7a96e80c516ba575c2`.

- `HCCfull.rds`: M **3,519×50**, CN **20,339×50**, T **20,356×50**, P **5,585×50**, PP **7,889×50** features × samples; T and P have all **50 identical sample column labels** in identical order.
- `HCCfullz.rds`: T **19,768×50**, P **5,585×50**, other layers as above, **50 T/P common sample labels**.
- `HCCdata.rds`: **50 samples in rows** × reduced features: M 24, CN 292, T 188, P 116, PP 158. Source regression now checks the proper axis; previous column-width mismatch was an audit orientation error, not absence of matched sample axes.
- Reproducible evidence: `research/hcc/feasibility/NG_SUTER2022_PUBLIC_RDS_EVIDENCE.json`, `scripts/hcc/audit_ng_public_multiomics_rds.R`, `scripts/hcc/verify_ng_public_rds_evidence.py` and related synthetic tests.

**Research decision:** GO for *public processed multiomics matrix access / exploratory feasibility*, not clinical/mechanistic inference. **V02 is PARTIAL as of Stage-17** after exact biopsy/patient source mapping; source expression/protein units, preprocessing, imputation and batch documentation remain pending; V04 PARTIAL, V12 PARTIAL, V13 BLOCKED. The source cohort and its clustering networks are previously published prior art, **not independent replication**. Do not claim 50 distinct patients (Suter reports 48).
**Next:** map the 50 RDS source sample labels to Ng Supplementary Data 1 tumor biopsy and patient keys inside a protected CI step, export aggregate matched/unique-patient counts only, then evaluate etiology and missingness as power constraints.

## Stage-17: protected exact Ng 2022 biopsy-to-patient joining on public quantitative matrices (2026-10-31)

**Empirical strict-source PASS:** [GitHub Actions #41160770924](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41160770924). A private, ephemeral 50-sample axis extracted from the immutable cbg-ethz/HCC `HCCfull.rds` was compared in-memory with Ng 2022's original SHA256-pinned Supplementary Data 1 (`edf1d2afd3baaedd56edee8536c16f10e4a62cb6cdac929863a062671ee73d25`). **All 50/50 source biopsy labels matched exactly under strict normalized identity**, **zero unmatched**, all 50 annotated as RNA-seq=Y and proteomics=Y; **48 unique nonmissing patient IDs**. Source patient duplication is exactly two excess biopsy rows, consistent with Suter's report of two patients with two HCC tumors. No raw IDs, clinical row values, sample-label hashes, RDS or private key files are published; the runner's ephemeral labels were deleted.

**Data audit implication:** Processed numeric T/P matrices have 50 identical sample labels *and* those labels now map to real source biopsies/48 unique patients. This removes the previous uncertainty about file-level same-biopsy identity and enables **exploratory paired RNA–protein feasibility using 50 biopsies/48 patients**, while clinical inference still needs careful modelling for duplicated patients.

**V02 scientific gate upgraded from BLOCKED to PARTIAL, not PASS.** Missing verification: original RNA/protein normalization units and directions, Suter's feature filtering/imputation/batch correction (S5 Appendix), clinical etiology definition and distribution of groups among 48 patients, outcome event/censoring codebook, independent external replication. Source-identity verification alone does not establish mechanistic effects or clinical significance.

Code: `scripts/hcc/probe_ng_suter_clinical_sample_join.py`, `scripts/hcc/audit_ng_public_multiomics_rds.R`; evidence: `research/hcc/feasibility/NG_SUTER2022_PUBLIC_RDS_EVIDENCE.json`.

## Stage-18: 48 matched-patient source etiology distribution — corrected and CI validated (2026-10-31)

**Source-level direct result**: [Ng/Suter patient-join CI #41161332568](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41161332568) PASS on original SHA256-validated Ng Supplemental Data 1 and public quantitative RDS, 50/50 matched biopsies, 48 distinct patient keys. All exposed values are aggregate source-defined label categories; patient keys, individual phenotypes and small cells (<5) are not published.

| Etiology label within the 48 distinct patients | Source-category candidate n |
| --- | ---: |
| ALD only | 12 |
| HCV only | 11 |
| NAFLD only | 7 |
| ALD + HCV | 7 |
| HBV only / ALD + HBV / ALD + NAFLD / source-other | suppressed (<5 each) |

**Non-exclusive exposures across mixed categories**: source-labelled ALD **24**, HCV **18**, HBV **7**, NAFLD **9** of 48; these exposures overlap and cannot be summed as independent patients. No per-patient conflicting etiology labels among the two patients with multiple biopsies.

**Bug found and corrected**: `Non-alcoholic liver disease (NAFLD)` must not be categorized as ALD merely because it contains the substring `alcoholic`. Source classification was fixed to require explicit `ALD` or alcoholic liver disease not preceded by `non-`, with a synthetic regression reproducing the source wording. The preliminary ALD=31 inclusive tally was an **invalid overcount**; the validated corrected ALD any-exposure value is **24**.

**Analysis and power decisions:** 12 ALD-only vs 11 HCV-only are the most immediate independent-patient etiology contrast. This supports descriptive and effect-size feasibility but **not** a priori confirmatory transcriptomic/proteomic significance, given the small n and potentially thousands of tested features. Mixed etiologies cannot be naively assigned to either exclusive group. HBV-only remains a small suppressed stratum. Evaluate tissue stage, age, sex, batch and processed data transformations before choosing statistical tests.

**Gate decisions:** V02 PARTIAL (50 exact biopsies/48 source patients and actual quantitative T/P matrices), V04 PARTIAL (source categories exact but clinical etiology codebook unreviewed), V12 PARTIAL, V13 BLOCKED. PASS claims for mechanism or recurrence **remain 0**. Evidence details in `NG_SUTER2022_PUBLIC_RDS_EVIDENCE.json`.

## Stage-19: 48-patient ALD-only versus HCV-only confounder + RNA/protein numeric QC (2026-10-31)

**Actual source-based CI PASS**: [clinical source balance #41161928470](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41161928470) and [numeric/feature QC #41162144728](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41162144728). Source-specified Ng clinical workbook SHA256 and pinned Suter RDS revision preserved; no patient rows/IDs, sensitive raw expression values, or clinical event assumptions were published.

| Covariate (deduplicated patient) | ALD only | HCV only |
|---|---:|---:|
| N patients | 12 | 11 |
| Age median [Q1, Q3], years | 67 [62, 74.75] | 69 [56, 76] |
| BCLC stage 0/A | 7 | 5 |
| BCLC B/C/D | 5 | 6 |
| Male | 12 | 9 |

Clinical workbook allowlist headers retrieved: age, sex, BCLC, Edmondson grade, MELD and cirrhosis. Treatment history and sequencing/proteomics batch columns were **not identified** under conservative header rules; do not assume these variables were measured or confounding removed. A within-patient discordant Edmondson annotation exists in the multi-biopsy cases (small count suppressed) and must not be collapsed arbitrarily.

**Raw numeric QC on public preprocessed HCCfull.rds**: RNA T 20,356 features×50 biopsies, all Ensembl-format gene IDs, zero NA/nonfinite entries, 416,949 negative observations, 158,871 exact zeros, overall median 0. Protein P 5,585 features×50 biopsies, zero NA/nonfinite entries, 142,614 negative observations, 0 exact zeros, overall median -0.0211. **Exact RNA–protein feature label intersection = 0** because identifiers differ, *not because no corresponding biological genes exist*. A reviewed mapping between Ensembl gene IDs and proteome accessions/symbols is a hard prerequisite for paired gene-level concordance.

**Scientific interpretation**: The processed matrices are centered/transformed in appearance and cannot be treated as raw RNA counts or raw mass-spectrometric intensities. Exact transformation, normalization and missing-value imputation cannot be proven from the numeric distribution and require direct review of Suter 2022 S5 Appendix and original Ng methods. Avoid DESeq2/count-based modelling on these processed RDS values. For now GO *descriptive sample pairing and covariate/omics QC only*, HOLD *ALD-versus-HCV high-dimensional differential effects*, BLOCKED *novel and independent prognostic claims*. V02 and V04 remain PARTIAL, V12 PARTIAL, V13 BLOCKED.

Evidence record: `research/hcc/feasibility/NG_SUTER2022_STAGE19_CLINICAL_NUMERIC_QC.json`. New code: `scripts/hcc/probe_ng_suter_covariate_balance.py` and `scripts/hcc/assess_ng_rds_quantitative_qc.py`, with synthetic regression tests integrated in `.github/workflows/hcc-ng-public-rds.yml`.

## Stage-20: original bnInfo source gene crosswalk for selected RNA/protein features (2026-10-31)

**Actual Actions PASS:** [Ng RDS source mapping CI #41162706947](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41162706947). The original public `HCCinfo.rds` and reduced `HCCdata.rds` were downloaded with immutable Git source pin `438f52271d0a86a6701baa7a96e80c516ba575c2` and SHA256 checks; `bnInfo`-like `allnames` tables provide source-defined gene symbols for the selection. Exactly 188 of 188 T and 116 of 116 P reduced feature keys map; their unique gene symbols have **114 shared genes, all 114 unambiguous one-to-one anchors**. No raw gene names, patient keys or expression values were exported.

**Important boundary:** This is a **source-based selected-feature crosswalk**, NOT a whole-genome 20,356-transcript vs 5,585-protein mapped dataset. Original 188/116 features were selected by Suter et al. to study HCC clustering and networks. This carries strong selection/analysis reuse risk; an inferential concordance p-value calculated from those 114 features alone would be post-selection and cannot be interpreted as an independent discovery. Review source S5 preprocessing, feature selection and full Ensembl/UniProt mapping with independent annotations before high-dimensional RNA–protein significance.

**Gates:** technical reduced crosswalk GO; V02 PARTIAL, V03 PARTIAL, V04 PARTIAL, V05 PARTIAL, V07 concordance NOT_STARTED, V12 PARTIAL, V13 BLOCKED. **Biological PASS remains 0**. Evidence detail: `research/hcc/feasibility/NG_SUTER2022_STAGE19_CLINICAL_NUMERIC_QC.json` (Stage20 embedded) and scripts `audit_ng_suter_reduced_feature_mapping.R` / `verify_ng_reduced_feature_mapping.py`.

## Stage-21: source-verified ALD12/HCV11 confounder imbalance and clinical viability HOLD (2026-10-31)

Original 50 Suter processed biopsies were strictly reconciled to Ng Supplementary Data 1, 48 unique patients. The [actual pinned-source full CI #41163101027](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41163101027) verified ALD-only 12 patients and HCV-only 11 patients with the following **absolute standardized mean differences**: age **0.077**, male sex **0.667**, BCLC 0/A vs B/C/D **0.260**. Source age medians 67 (IQR 62–74.75) vs 69 (IQR 56–76); early BCLC 7/12 vs 5/11; male 12/12 vs 9/11. These are descriptive imbalance metrics, not association significance tests.

The official source `child` column was newly recognized as a *candidate* Child–Pugh category in a subsequent source-verified CI [#41163221941](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41163221941). Class A/B/C was parseable for only **8 ALD** and **3 HCV** patients, so the initially calculated Child–Pugh SMD of 0.343 is **not reliable**. A fail-closed safeguard now suppresses Child–Pugh SMD where either group has fewer than five interpretable values; its final regression CI is separately tracked. Source MELD medians were **11** (12 numeric ALD patients) and **8** (10 numeric HCV patients); yes-labelled cirrhosis counts were **11** and **8**, respectively. Interpret these only within the original source metadata, not as adjusted patient-risk predictions.

The clinical source has no directly verified technical batch variable in the protected allowlist, and treatment line fields cannot be used as pre-biopsy baseline predictors. Suter reports normalized Ng input followed by imputation and batch correction where applicable (Methods 2.1 and S5 Appendix; https://doi.org/10.1371/journal.pcbi.1009767). Ng used upper-quartile normalized RSEM expected counts and protein intensities expressed on log2 relative scale (https://doi.org/10.1038/s41467-022-29960-8); the released RDS contains negative values and must **not** be treated as unmodified counts for DESeq2.

**Decision:** GO descriptive cohort and omics QC, HOLD direct high-dimensional ALD-versus-HCV differential expression, pathway discovery p-values, and mechanism/clinical claims. Source-sex and BCLC confounding, small independent patient n=23, partial Child–Pugh availability, residual batch confounding and previously selected 114 reduced gene anchors preclude confirmatory GO. V02 PARTIAL, V04 PARTIAL, V12 PARTIAL, V13 BLOCKED. Next: audit whole-profile feature mapping and preprocessing/technical batch provenance and design independently replicated cohorts; avoid overfitted covariate models in 23 patients.

## Stage-22: full-resolution RNA/protein identifier namespaces and original Ng benchmark (2026-10-31)

- **Empirical pinned-source CI** [#41170474747](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41170474747) PASS. Ng/Suter `HCCfull.rds` transcriptome T has **20,356/20,356 ENSG-like unique Ensembl gene labels**, proteome P **5,530/5,585 canonical UniProt-accession-like labels** (5,585 unique protein labels). Among the previously selected reduced features, **188 T and 116 P identifiers** match their own full source feature identifiers. Raw ENSG vs UniProt name overlap is zero by namespace, not evidence of biological nonoverlap.
- **Primary published benchmark**: Ng et al., *Nat Commun* (2022), DOI [10.1038/s41467-022-29960-8](https://doi.org/10.1038/s41467-022-29960-8), methods mapped protein UniProt accessions to Ensembl IDs for **5,490 comparable genes/proteins** (quadrant analyses) and **5,481 genes** for mRNA–protein correlation. Do not conflate these universes or equate 5,530 UniProt-shaped source labels with 5,530 valid mapped gene pairs.
- **Live data source screening**: UniProt reviewed human release currently has **5,515** of the 5,530 canonical-looking source protein accessions; however its `xref_ensembl` TSV field did **not** contain ENSG gene IDs (`ENST` transcripts are not `ENSG` gene IDs). First TSV-only [CI](https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41170858103) therefore failed-closed at zero *validated gene mappings*. Public JSON crossreference `GeneId` lookup is being validated in a separately isolated workflow, with 40-accession batches, retries for transient HTTP 500/429 and aggregate-only outputs. **Do not promote a full mapping count until a successful complete job**.
- **Processing constraints**: Ng original RNA differential tests required raw RSEM expected counts as input to edgeR TMM/quasi-likelihood; protein used limma on log2 relative expression, BH FDR. The existing Suter full RDS is processed/imputed/batch-corrected data, not raw gene counts. Rerunning edgeR/DESeq2 count models on the published processed matrices would be statistically unsupported. Ng source already reported dysregulated pathways and 5,481 gene-level correlations, so merely reproducing correlation is not a new biological result.
- **Gate decisions:** technical source namespace GO; V02/V03/V04 PARTIAL, V07 NOT_STARTED, V12 PARTIAL, V13 BLOCKED. Full stable UniProt–ENSG mapping, batch provenance, inferential etiology interaction and independent validation **HOLD**.
