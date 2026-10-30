# HCC Research Decision Log (immutable decisions, append-only)

**Purpose:** preserve decisions after each analysis, especially why a checkpoint moved or a manuscript claim changed. Never rewrite older decisions; append a dated correction.

| UTC/KST date | Decision ID | Gate(s) | Evidence | Decision | Meaning for thesis |
|---|---|---|---|---|---|
| 2026-10-30 | D-001 | V00–V02 | https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41053421959 | Ng 51 triple-omic assay-flag biopsies = 49 distinct patients. Actual gene×sample matrix join is **not verified**. | Matched quantitative mechanistic claims are blocked; original paper's N must not substitute for a quantified join. |
| 2026-10-30 | D-002 | V04–V05 | https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41053421959 | Within 49 flagged patients: ALD-only 13, HCV-only 11, NAFLD-only 7, HBV-only 4; mixed etiologies separate. | Etiology effect modification is exploratory. Predeclare interaction contrast and assess BCLC, treatment, sex/age and platform confounding. |
| 2026-10-30 | D-003 | V01–V02 | https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41055576655 | PDC Protein Assembly proteome/phosphoproteome files and PRIDE SWATH/Progenesis candidates identified. | File presence is not proven RNA–protein join; avoid promoting P02 and V02. |
| 2026-10-30 | D-004 | V12–V13 | https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41052039696 | GSE14520 has 242 samples with raw Recurr status/time columns; unique patient rows and event code not audited. | No independent recurrence model/clinical significance claim; preserve RQ3 as conditional. |
| 2026-10-30 | D-005 | V06–V09 | https://doi.org/10.1038/s41586-019-0987-8 ; https://doi.org/10.1016/j.cell.2019.08.052 ; https://doi.org/10.1038/s41467-022-29960-8 | SOAT1/cholesterol, PYCR2/ADH1A, and HCC metabolic RNA/protein differences are **published prior findings**. | These are positive controls; novelty requires independent cross-etiology portability or biological source evidence. |
| 2026-10-30 | D-006 | All | research/hcc/validation/HCC_VALIDATION_MASTER.md | Separate method feasibility, effect significance, biological meaning, novelty and clinical utility. | Never claim high research significance because q<0.05 or CI is green. |

## Entry template

Date / ID / repo commit / checkpoint / dataset / source SHA256 / patient N / event N (if relevant) / outcome definition / analysis and plot artifact path / effect estimate / 95% CI / multiple-test correction / independent validation / evidence sufficiency / remaining bias / status change / exact allowed claim / next step.

## Automatic continuation rule

For a future request such as “다음”, first inspect the protocol, checkpoint registry and **latest CI**. Complete the first actionable public-data blocker. Save real quantitative evidence and update the ledger. If the hypothesis turns out null, negative, confounded, or restricted by access, log the limitation rather than manufacturing a PASS.

Study-governance PR is independent from the ongoing code/data-feasibility PR #251; merge or review both before treating them as available on main.

## 2026-10-30 — D-007: Evidence-gate hardening and Drive handoff

- Protocol governance PR: https://github.com/kimtk94/IS_Analysis_V3/pull/252
- Automated verification: https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41057719470 (12 offline tests PASS)
- Google Drive reference: https://docs.google.com/document/d/1Ipk7kS6ukGYqwgwhSWHICs42EBqpu3eRFd82mBOkaAE/edit
- Rule update: `V12`/`V13` clinical event/count gates alone do **not** authorize a biology-guided clinical significance claim when the biological validation prerequisites are incomplete. Full biological prerequisite set includes V00, V01, V02, V03, V04, V05, V06, V07, V08, V09, V14 and V15.
- To mark any gate PASS, record explicit reviewer approval and reproducible analysis artifact path under `results/hcc/`, SHA256 of source (64 hex), exact code commit (40 hex), and checked date. V02/V04/V10 additionally require independent patient N; V12/V13 require verified event N.
- The static evidence checker validates *format and declared evidence*, **not** remote artifact contents, statistical effects, new discovery, or prospective clinical usefulness. A scientific reviewer still needs to inspect outputs. Even every workflow gate being PASS only makes a claim eligible for review.
- The user-facing handoff document resides with the original feasibility audit; the Markdown protocol and machine-readable status remain canonical in GitHub. Existing conclusions: molecular and clinical findings unproven.

## 2026-10-30 — D-008: PDC tumor-normal lineage and 11-plex analysis, RNA–Protein gate still BLOCKED

- Evidence: https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41059030994 and code on feasibility PR #251.
- Public PDC `paginatedCasesSamplesAliquots` on PDC000198 (Proteome) and PDC000199 (Phosphoproteome) yielded an identical set of **171 case, 331 sample, 331 aliquot identifiers** in both studies. Output saves only aggregates, not case/aliquot identifiers.
- Per study: **165 Primary Tumor samples, 165 Solid Tissue Normal samples, 1 Not Reported sample**. Distinct cases with **both tumor and normal: 160**; 5 tumor-only, 5 normal-only, 1 unclassified case. Publication's 159 paired analysis subjects should not be replaced with 160 without adjudicating inclusion/QC exclusions.
- Experimental metadata includes **33 TMT experiment numbers × 11 labels = 363 unique channel slots**, with 726 run metadata records; not 363 patients or 726 quantification observations. Header 667 (Protein) / 334 (Phosphosite) columns cannot be joined by direct full case/sample/aliquot string matching (0 matches).
- Scientific implication: this resolves which PDC cases are linked at the biospecimen layer, not the quantitative-matrix value layer, and **does not contain a BioSino RNA expression join**. This strengthens V01's source-level evidence but leaves **V02 BLOCKED**, all biological/clinical primary claims ineligible.
- Next required test: aliquot-to-(experiment, label) uniqueness and bridge/reference channels; replicate vs biological sample distinction; quantitative TSV header parsing; BioSino OEP000321 RNA gene-by-sample mapping and compatible patient keys.

## 2026-10-30 — D-009: Gao TMT analyte collision resolved in public metadata, V02 NOT promoted

- Source evidence: https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41062185538 ; official PDC `paginatedCasesSamplesAliquots` and HCC feasibility PR #251.
- In each PDC study, 726 `aliquot_run_metadata` rows span 33 experiment numbers. If keyed only by experiment+TMT label, 330 of 363 channel slots are ambiguously assigned to two aliquots. **The ambiguity is due to mixing two assay analytes:** `Proteome` and `Phosphoproteome`.
- Split by analyte: each class has **363 records, 363 experiment×label slots, 331 distinct aliquots, and 0 ambiguous slots**. `replicate_number` is blank in all 726 records; `fraction` alone leaves 60 ambiguous slots. The study field, analyte, experiment and label must not be conflated.
- Protein quantitative TSV has 330 Log Ratio + 330 Unshared Log Ratio columns (667 total including metadata); phosphosite has 330 Log Ratio columns (334 total including metadata). PDC values are **relative TMT ratios**, not absolute abundance. No public assay-specific TSV column→case mapping has yet been verified.
- NODE `OEP000321` was probed; HTTP 200 returned only a 906-byte public HTML shell, no source-accession text and no processed RNA matrix URL. Secondary publications report RSEM-quantified RNA from OEP000321, but actual file/permission/IDs were not verified.
- **Scientific gate implication:** PDC Proteome/Phosphoproteome experimental lineage is materially better resolved. Nonetheless cross-modal quantitative **V02=BLOCKED**, RQ1 inferential signal not testable, recurrence clinical claim ineligible. Next: match protein TSV normalized columns to `(analyte,experiment,TMT label)` mapping; locate and legally obtain processed OEP000321 RNA and verify same-patient overlap.

## 2026-10-30 — D-010: PDC native matrix API and HCCDB25 access (V01 partial / V02 blocked)

- Source evidence: https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41064173332 (49 regression tests passed). Official API documentation: https://pdc.cancer.gov/pdc/publicapi-documentation
- A **better PDC integration route** is documented: `quantDataMatrix` returns a gene-by-aliquot quantitative matrix and the PDC GCT heatmap embeds aliquot-level metadata. Therefore prefer direct aliquot-level normalized log2-ratio values over speculative experiment/channel-based TSV column mapping. However, real PDC000198/PDC000199 direct API requests both **timed out after 25 seconds** in CI. Actual quantitative matrix headers and patient joins were not recovered by that operation.
- **Alternate RNA source:** Gao `OEP000321` is archived as **HCCDB25** in HCCDB v2.0, documented by HCCDB and independent publications. Public `http://lifeome.net:809/` responded HTTP 206, **6,305 bytes of generic HTML**, with no HCCDB25 literal accession or verified matrix file URL. HTTPS routes failed (SSL or connection).
- Data access decision: PDC native matrix **UNVERIFIED**, HCCDB25 processed RNA file **UNVERIFIED**, V02 **BLOCKED**. Do not label timeout as absent data, HTML as an expression matrix, or the paper's 159 individuals as verified RNA–protein patient pairs.
- Follow-up order: inspect HCCDB front-end public script manifests / documented download API, probe PDC portal GCT availability, then verify same biological patient/biopsy ID across *actual quantitative values*. Avoid recurring high-load `quantDataMatrix` calls on every PR push; manually trigger bounded calls only when needed.

## 2026-10-30 — D-011: HCCDB25 frontend route discovered, no RNA file verified

- Evidence: https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41066051263 (focused public endpoint CI; 3 tests PASS).
- HCCDB v2.0 public site has **six referenced JS assets** in its main HTML shell (HTTP 206, 6,305 bytes). Of five inspected, the small `/static/js/app.3f5b30f0.js` has `/download` and `/repository` literals; they are **client-side navigation paths**, not confirmed quantitative RNA download URLs.
- No direct HCCDB25 sample-by-gene expression file, sample metadata, or RNA→PDC case join was recovered. Literature provenance HCCDB25=Gao OEP000321 is maintained as a source relationship, not a verified current binary/file transfer.
- **Decision:** V02 BLOCKED. Next: authorized HCCDB25 repository navigation or official export; if access remains blocked, request source data according to original publication/institution policy and invoke independent-cohort RQ4 only when cross-omics pairing is truly impossible.
- Engineering adjustment: PDC `quantDataMatrix` heavy/time-out-prone probes are manual-only in the standard CI; focused HCCDB25 route workflow allows quick low-risk verification.

## 2026-10-30 — D-012: Cheng 2023 RNA–Protein integration is prior art; PDC TMT-126 candidate reference identified

- New source: https://doi.org/10.1093/bib/bbac601 (Cheng et al., *Briefings in Bioinformatics* 2023). Authors reanalyzed Gao CHCC-HBV by excluding incomplete multi-omics samples; **158 paired patients**, RNA `OEP000321` UQ_FPKM plus Gao protein/phosphoprotein supplements. Reported **408 concordant deregulated mRNA/protein genes** (224 up, 184 down), 3 subtypes, independent validation. These are **published findings**, not results of our pipeline.
- The verified PDC **metadata-only** candidate control structure: `TMT_126` appears in 33 experiments and references one distinct aliquot per Proteome/Phosphoproteome analyte; 330 other channel slots correspond to 165 tumor/165 adjacent normal aliquots. Source: https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41068016825 (56 tests PASS). This does not prove a pooled reference in authors' protocols or map full quantitative TSV values.
- Sample denominator reconciliation remains **159 original paper / 171 PDC repository / 160 PDC both tumor-and-normal / 158 Cheng analytic**, distinct concepts requiring cohort selection verification.
- HCCDB25 repository SPA route `http://lifeome.net:809/#/repository?type=buke` corroborated by published research (https://doi.org/10.3389/fcell.2024.1407428); successful public Vue router investigation: https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41068016603. No RNA expression file downloaded.
- **V06** paper-replication control expanded: attempt reproduction of Cheng 408 genes and 3 clusters only after source identity and assay normalization verification. **V15** novelty gate: do not market these results as new; favor independent etiologic transferability and patient-level cell context.
- **V02 stays BLOCKED**, since neither direct `OEP000321` normalized matrix nor a full PDC case→quantitative-column and RNA–protein common-ID join has passed. Do not confuse green CI or published 158 with our joined patient N.

## 2026-10-30 — D-013: 158 eligible but 156 matched tumor samples in Cheng integrative clustering

- Primary paper Methods: https://academic.oup.com/bib/article/24/2/bbac601/7025463 . The study initially selected **158 paired multi-omics HBV-HCC** patients after data completeness filtering; **2 RNA/protein sample-ID-unmatched samples were removed** before the integrative clustering, yielding **156 tumor samples**. The paired 158 eligibility denominator must not be presented as the 156-sample final integrated cluster denominator.
- Required V02 evidence when datasets are obtained: `N_paper_eligible=158` (reproduce or explain), `N_rna_protein_matched_for_cluster=156` (reproduce or explain), explicit identifiers/filters audited in a nonpublic workspace, patient/sample counts retained as sanitized summaries in Git; do not assert counts as results before direct reanalysis.
- Required V06 positive control and V15 novelty restrictions remain as in D-012: Cheng's 408 genes and 3 subtypes are prior art.
- Current V02 **BLOCKED** because neither exact UQ_FPKM raw processed matrix nor real matched patient manifest has been downloaded and independently audited.

## 2026-10-30 — D-014: external 2023/2026 HCC prior art + alternative proteome

- 2023 *Cell Reports Medicine*, DOI https://doi.org/10.1016/j.xcrm.2023.101315: 160 recruited, 152-patient quantitative proteomic supplement `mmc3.xlsx` (~16.2 MB), 132-patient phosphoproteomic supplement `mmc6.xlsx` (~66.2 MB), RNA differential-result supplement `mmc7.xlsx`. Actual public workbook download is being separately verified; mere file listing never constitutes independent quantitative patient matching. The study already cross-validates molecular subtypes in Jiang, Gao, and Ng cohorts.
- 2026 *Nature Communications*, DOI https://doi.org/10.1038/s41467-026-73212-y: 529 integrated HCC samples classified into **nine genomic/immune subtypes**, external validation n=807. Generic new HCC subtype discovery is highly crowded.
- **V15 novelty decision:** do not market three cross-cohort HCC proteomic subclasses or a nine-subtype genetic classification as newly discovered. Focus RQ1 on etiologic interaction and robust external transfer of RNA–Protein pathway concordance, with proper clinical and cell-level validation. Do not promote V02 from publication counts: only actual value+sample join is acceptable.

## 2026-10-30 — D-015: Xing2023 official raw XLSX source and patient-level QA

- **Official PMC distribution switch:** NCBI discontinued the legacy PMC FTP/metadata distribution in August 2026. Modern official anonymous Cloud metadata at `https://pmc-oa-opendata.s3.amazonaws.com/metadata/PMC10783603.1.json` returned HTTP 200 with 16 media objects (6 XLSX). Its `s3://pmc-oa-opendata/` media object references must be resolved to the documented anonymous HTTPS S3 location; a deprecated PMC URL returning HTTP 200 HTML is NOT an XLSX.
- **Real public source files obtained and validated by GitHub Actions**, all from the same Xing et al. *Cell Reports Medicine* 2023 article DOI 10.1016/j.xcrm.2023.101315: `mmc3.xlsx` **16,978,396 bytes SHA-256 7da732a933eba7a84936b8033d37564fd27125e29d059f9d7e306cf456965e17** (152-patient proteome); `mmc2.xlsx` **50,781 bytes** (clinical Table S1); `mmc7.xlsx` **226,184 bytes** (RNA DE and phosphosite Table S6). All OOXML verified; source hashes and evidence in `research/hcc/feasibility/XING2023_EXTERNAL_QC.json`. Original copyrighted workbook bytes were not committed to public GitHub (publisher CC BY-NC-ND).
- **Actual protein matrix inspected:** `1. 6512 proteins` 6,512 feature rows × 304 biological sample columns. Using in-memory T/N-labeled header IDs, **304 sample columns form 152 Tumor/Normal sample-key pairs**; these are *not* independent clinical-ID joins. Of 1,979,648 protein cells, 1,311,220 were finite numerical values (**66.235%**). Median protein measurements per sample **4,349** (min **2,951**, max **5,459**); **4,383 proteins** measured in at least 152 samples, **1,102** in all 304. Actual evidence: https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41075928463
- **Clinical corroboration is NOT yet achieved:** in `mmc2.xlsx` detailed clinicopathological worksheet `A1:AG163`, 160 data rows were read; strict numeric ID matching with the 152 protein sample-key pairs matched at most **2**, potentially coincidental clinical values; stable prefixed-ID strategy matched **0**. Thus DO NOT claim the 152 protein pairs independently join to 152 patient phenotypes. Evidence https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41076963271 and https://github.com/kimtk94/IS_Analysis_V3/actions/runs/41077195020.
- The `mmc7.xlsx` `1. RNA-seq` sheet `A1:G1643` is a **differential-expression/summary table**, not a patient-level RNA quantitative matrix; the paper reports **57 RNA-seq assays** across the 160 enrolled patients, while proteomics QC yielded **152 patients**, phosphoproteome 132 and WES 58. The RNA–Protein intersection is **NOT inferred to be 57**.
- **Decisions:** external protein raw processed quantitative **source + numeric QC PASS**, T/N header pairing **PARTIAL** pending author/clinic ID confirmation, patient-level clinical join **HOLD**, cross-modal RNA×Protein **V02 remains BLOCKED**. Prior 152-proteomic subtype discovery and external validation are not new thesis results; V15 novelty constraint preserved.
- **Next data action:** inspect the detailed clinical workbook's ID-field definitions and non-sensitive format fingerprints; obtain patient-level RNA gene×sample matrix from study-specific authorized repository or custodian and perform a reproducible independent same-biopsy join. Do not treat Excel supplementary RNA differential results as RNA expression values.
